classdef BayesOptRandom
    % quick implementation of BayesOpt relying on GPStuff for the GP and
    % NL_OPT for acquisition function optim. Most options are fixed as
    % this version is only for experimentation purposes
    
    % Random features to approximate the GP with linear regression
    
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['BayesOptRandom' '_'  num2str(optimizerInput.nbInitSamplesBO) '_' num2str(optimizerInput.nbEvalsBO) '_' num2str(optimizerInput.maxNbSamplesBO) '_' optimizerInput.gpHyperOption '_' ...
                optimizerInput.kernelType '_' optimizerInput.featureName '_' optimizerInput.yCenteringType '_' optimizerInput.hypOptimAfter];
        end
        
        %%
        function [allVals] = optimizeStruct(optimizerInput, func)
            if(~isfield(optimizerInput, 'videoFile'))
                optimizerInput.videoFile = [];
            end
            
            [allVals] = static_optimization_algs.BayesOptRandom.optimize(optimizerInput.functionBounds, optimizerInput.nbInitSamplesBO, optimizerInput.nbEvalsBO, optimizerInput.maxNbSamplesBO, ...
                func, optimizerInput.yCenteringType, optimizerInput.videoFile);
        end
        
        %% main function
        function [allVals] = optimize(functionBounds, nbInitEvals, nbTotalEvals, maxNbSamples, fun, yCenteringType, videoFile)
            usedSamples = [];
            usedVals = [];
            allVals = [];
            dim = size(functionBounds, 2);
            lb = functionBounds(1, :);
            ub = functionBounds(2, :);
            %%% initializing GP
            lik = lik_gaussian('sigma2', 1e-6);
            pn = prior_fixed();
            lik = lik_gaussian(lik,'sigma2_prior', pn);
            gpcf = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1);
            % hyper-param priors of GP
            pl = prior_unif();
            pm = prior_sqrtunif();
%             pl = prior_fixed();
%             pm = prior_fixed();
            gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
            gp = gp_set('lik', lik, 'cf', gpcf);
            
            %%% animation
            if(~isempty(videoFile))
                open(videoFile);
            end
            
            funGp = [];
            for iter = 1:nbTotalEvals
                x = usedSamples;
                y = usedVals;
                %%% select next point
                if(iter <= nbInitEvals)
                    % generate a point from the unif distrib
                    newSample = rand(1, dim) .* (ub - lb) + lb;
                else
                    % optimize acquisition
                    %%% GP hyper-param
                    if(strcmp(yCenteringType, 'min'))
                        y = usedVals - min(usedVals);
                    elseif(strcmp(yCenteringType, 'mean'))
                        y = usedVals - mean(usedVals);
                    elseif(strcmp(yCenteringType, 'max'))
                        y = usedVals - max(usedVals);
                    else
                        error('Unrecognized y centering type in LocalBayes');
                    end
                    [newSample, gp, funGp] = static_optimization_algs.BayesOptRandom.optimizeAcquisition...
                        (x, y, functionBounds, gp);
                end
                
                newVal = fun.eval(newSample);
                fprintf('iteration %d\t val %f\n', iter, newVal);
                allVals = [allVals; newVal];
                
                usedSamples = [usedSamples; newSample];
                usedVals = [usedVals; newVal];
                if(length(usedVals) > maxNbSamples)
                    usedSamples = usedSamples(end-maxNbSamples+1:end, :);
                    usedVals = usedVals(end-maxNbSamples+1:end);
                end
                
                %%% current distrib plotting
                if(~isempty(videoFile))
                    frame = static_optimization_algs.BayesOptRandom.getFrame(fun, x, funGp, iter);
                    writeVideo(videoFile, frame);
                end
            end
            
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        function features = getFeatures(x, W, b, magn)
            features = sqrt(2 * magn / size(x, 1)) * cos(bsxfun(@plus, W * x, b));
        end
        
        function g = acquisition(theta, x, W, b, magn)
            g = theta' * static_optimization_algs.BayesOptRandom.getFeatures(x', W, b, magn);
        end
        
        function [newSample, gp, fun] = optimizeAcquisition(x, y, functionBounds, gp)
            opt = optimset('TolFun', 1e-3, 'TolX', 1e-3);
            gp = gp_optim(gp, x, y, 'opt', opt);
            fprintf('noise %f, lengthscale %f, magnitude %f\n', gp.lik.sigma2, gp.cf{1}.lengthScale, gp.cf{1}.magnSigma2);
            
            % building random function from GP
            dim = size(x, 2);
            nbFeatures = 500;
            W = randn(nbFeatures, dim) / gp.cf{1}.lengthScale;
            b = rand(nbFeatures, 1) * 2 * pi;
            features = static_optimization_algs.BayesOptRandom.getFeatures(x', W, b, gp.cf{1}.magnSigma2);
            A = features * features' + gp.lik.sigma2^2 * eye(nbFeatures);
            cholA = chol(A);
            invCholTheta = gp.lik.sigma2 * (cholA \ eye(nbFeatures));
            meanTheta = cholA \ ((cholA' \ features) * y);
            theta = meanTheta + invCholTheta * randn(nbFeatures, 1);
            fun = @(x) static_optimization_algs.BayesOptRandom.acquisition(theta, x, W, b, gp.cf{1}.magnSigma2);
            
            % Calling NL_OPT
            optOption.max_objective = fun;
            optOption.algorithm = NLOPT_GN_DIRECT_L;
            %             optOption.algorithm = NLOPT_LN_BOBYQA;
            %             optOption.algorithm = NLOPT_LD_LBFGS;
            %             obj.optOption.algorithm = NLOPT_LD_MMA;
            %             obj.optOption.algorithm = NLOPT_LD_TNEWTON;
            %             obj.optOption.algorithm = NLOPT_LD_SLSQP;
            %             obj.optOption.algorithm = NLOPT_LD_VAR1;
            optOption.verbose = 0;
            optOption.lower_bounds = functionBounds(1, :);
            optOption.upper_bounds = functionBounds(2, :);
            optOption.maxeval = 2000;
            %             [~, bestIndex] = max(y);
            %             params0 = x(bestIndex, :);
            params0 = (functionBounds(1, :) + functionBounds(2, :)) * .5;
            [newSample, ~, retCode] = nlopt_optimize(optOption, params0);
            if(retCode < 0)
                warning('Error in NL_OPT optimization');
            end
            optOption.maxeval = 200;
            optOption.algorithm = NLOPT_LN_BOBYQA;
            [newSample, ~, retCode] = nlopt_optimize(optOption, newSample);
            if(retCode < 0)
                warning('Error in NL_OPT optimization');
            end
        end
        
        function frame = getFrame(fun, samples, funGp, iter)
            persistent currPlot;
            if(isempty(currPlot))
                currPlot = figure;
            end
            
            %function and distrib
            subplot(1, 2, 1);
            fun.plot();
            hold on;
            if(~isempty(samples))
                plot(samples(:, 1), samples(:, 2), '*r');
            end
            hold off;
            xLimits = xlim;
            yLimits = ylim;
            
            subplot(1, 2, 2);
            if(~isempty(funGp))
                nbPoint = 100;                
                line1 = xLimits(1):(xLimits(2)-xLimits(1))/nbPoint:xLimits(2);
                line2 = yLimits(1):(yLimits(2)-yLimits(1))/nbPoint:yLimits(2);
                [x1, x2] = meshgrid(line1, line2);
                x = [x1(:) x2(:)];                                
                y = funGp(x);
                contour(line1, line2, reshape(y, length(line2), length(line1)));

            end
            text(-5, -5, ['iteration = ' num2str(iter)]);
            frame = getframe(gcf,[0 0 560 420]);
        end
        
    end
end
