classdef BayesOpt
    % quick implementation of BayesOpt relying on GPStuff for the GP and
    % NL_OPT for acquisition function optim. Most options are fixed as
    % this version is only for experimentation purposes
    
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['BayesOpt' '_'  num2str(optimizerInput.nbInitSamplesBO) '_' num2str(optimizerInput.nbEvalsBO) '_' num2str(optimizerInput.maxNbSamplesBO) '_' optimizerInput.gpHyperOption '_' ...
                optimizerInput.kernelType '_' optimizerInput.featureName '_' optimizerInput.yCenteringType '_' optimizerInput.hypOptimAfter];
        end
        
        %%
        function [allVals] = optimizeStruct(optimizerInput, func)
            if(~isfield(optimizerInput, 'videoFile'))
                optimizerInput.videoFile = [];
            end
            
            [allVals] = static_optimization_algs.BayesOpt.optimize(optimizerInput.functionBounds, optimizerInput.nbInitSamplesBO, optimizerInput.nbEvalsBO, optimizerInput.maxNbSamplesBO, ...
                func, optimizerInput.kernelType, optimizerInput.featureFunction, optimizerInput.yCenteringType, optimizerInput.gpHyperOption, optimizerInput.hypOptimAfter, optimizerInput.videoFile);
        end
        
        %% main function
        function [allVals] = optimize(functionBounds, nbInitEvals, nbTotalEvals, maxNbSamples, fun, kernelType, featureFunction, yCenteringType, gpHyperOption, hypOptimAfter, videoFile)
            usedSamples = [];
            usedVals = [];
            allVals = [];
            dim = size(functionBounds, 2);
            lb = functionBounds(1, :);
            ub = functionBounds(2, :);
            %%% initializing GP
            lik = lik_gaussian('sigma2', 1e-6);
            %             pn = prior_logunif();
            pn = prior_fixed();
            lik = lik_gaussian(lik,'sigma2_prior', pn);
            if(strcmp(kernelType, 'sexp'))
                gpcf = gpcf_sexp('lengthScale', 1, 'magnSigma2', 0.01);
                % hyper-param priors of GP
                pl = prior_unif();
                pm = prior_sqrtunif();
                gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
            elseif(strcmp(kernelType, 'linear'))
                gpcf = gpcf_linear('coeffSigma2', 1);
                % hyper-param priors of GP
                pl = prior_logunif();
                gpcf = gpcf_linear(gpcf, 'coeffSigma2_prior', pl);
            else
                error('unrecognized kernel in LocalBayes');
            end
            gp = gp_set('lik', lik, 'cf', gpcf);
            
            gp_rec = gp;
            
            %%% animation
            if(~isempty(videoFile))
                open(videoFile);
            end
            
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
                    if(~isempty(featureFunction))
                        x = featureFunction(x);
                    end
                    doHypOptim = false;
                    if(mod(iter - nbInitEvals - 1, hypOptimAfter) == 0)
                        doHypOptim = true;
                    end
                    [newSample, gp, gp_rec] = static_optimization_algs.BayesOpt.optimizeAcquisition...
                        (x, y, functionBounds, gp, gp_rec, doHypOptim, hypOptimAfter, gpHyperOption, featureFunction);
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
                    frame = static_optimization_algs.BayesOpt.getFrame(fun, x, y, gp, iter);
                    writeVideo(videoFile, frame);
                end
            end
            
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        function [newSample, gp, gp_rec] = optimizeAcquisition(x, y, functionBounds, gp, gp_rec, doHypOptim, hypOptimAfter, gpHyperOption, featureFunction)
            persistent indexHyper;
            if(isempty(indexHyper))
                indexHyper = 1;
            end
            % Setting the function
            if(doHypOptim)
                if(strcmp(gpHyperOption, 'MCMC'))
                    indexHyper = 1;
                    try
                        burnin = 100;
                        gp_rec = gp_mc(gp, x, y, 'nsamples', 2 * hypOptimAfter + burnin, 'display', 0);
                        gp_rec = thin(gp_rec, burnin, 2); %delete first burnin samples
                    catch
                        disp('error in mcmc sampling of gp parameters');
                    end
                    gp.lik.sigma2 = gp_rec.lik.sigma2(indexHyper);
                    try %params of sexp covariance func
                        gp.cf{1}.lengthScale = gp_rec.cf{1}.lengthScale(indexHyper, :);
                        gp.cf{1}.magnSigma2 = gp_rec.cf{1}.magnSigma2(indexHyper);
                    catch %params of linear covariance func
                        gp.cf{1}.coeffSigma2 = gp_rec.cf{1}.coeffSigma2(indexHyper);
                    end
%                     gp = gp_rec;
                    indexHyper = indexHyper + 1;
                elseif(strcmp(gpHyperOption, 'MAP'))
                    opt = optimset('TolFun', 1e-3, 'TolX', 1e-3);
                    gp = gp_optim(gp, x, y, 'opt', opt);
                else
                    error('wrong hyperparm option in LocalBayes');
                end
            elseif(strcmp(gpHyperOption, 'MCMC'))
                gp.lik.sigma2 = gp_rec.lik.sigma2(indexHyper);
                try %params of sexp covariance func
                    gp.cf{1}.lengthScale = gp_rec.cf{1}.lengthScale(indexHyper, :);
                    gp.cf{1}.magnSigma2 = gp_rec.cf{1}.magnSigma2(indexHyper);
                catch %params of linear covariance func
                    gp.cf{1}.coeffSigma2 = gp_rec.cf{1}.coeffSigma2(indexHyper);
                end
                indexHyper = indexHyper + 1;            
            end
            
            [gp.lik.sigma2, gp.cf{1}.lengthScale, gp.cf{1}.magnSigma2]
            
            % building acquisition function (thomson sampling)
            [~, C] = gp_trcov(gp, x);
            Cy = C \ y;
            cholC = chol(C);
            invCholC = cholC \ eye(size(C));
            
            %             func = @(xt) static_optimization_algs.GP.randFromCov(gp, xt, x, Cy, invCholC);            
%             func = @(xt) static_optimization_algs.GP.EILib(gp, xt, x, y, max(gp_pred(gp, x, y, x)));
            func = @(xt) static_optimization_algs.GP.EI(gp, xt, x, Cy, invCholC, max(static_optimization_algs.GP.predFromCov(gp, x, x, Cy, invCholC)));
%             static_optimization_algs.GP.ThompsonSampling(gp, []);
%             func = @(xt) static_optimization_algs.GP.ThompsonSampling(gp, xt, x, y);

            % Calling NL_OPT
            optOption.max_objective = func;
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
        
        function frame = getFrame(fun, samples, vals, gp, iter)
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
            center(1) = (xLimits(2) + xLimits(1)) / 2;
            length(1) = xLimits(2) - xLimits(1);
            yLimits = ylim;
            center(2) = (yLimits(2) + yLimits(1)) / 2;
            length(2) = yLimits(2) - yLimits(1);
            
            %gp
            subplot(1, 2, 2);
            if(~isempty(samples))
                static_optimization_algs.GP.plotGP(gp, samples, vals, center, length);
            end
            text(-5, -5, ['iteration = ' num2str(iter)]);
            frame = getframe(gcf,[0 0 560 420]);
        end
    end
end
