classdef BayesOptGPML
    % quick implementation of BayesOpt relying on GPStuff for the GP and
    % NL_OPT for acquisition function optim. Most options are hand defined.
    % This version is only for experimentation purposes
    
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['BayesOptGPML' '_'  num2str(optimizerInput.nbInitSamplesBO) '_' num2str(optimizerInput.nbEvalsBO) '_' num2str(optimizerInput.maxNbSamplesBO) '_' optimizerInput.gpHyperOption '_' ...
                optimizerInput.kernelType '_' optimizerInput.featureName '_' optimizerInput.yCenteringType];
        end
        
        %%
        function [allVals] = optimizeStruct(optimizerInput, func)
            if(~isfield(optimizerInput, 'videoFile'))
                optimizerInput.videoFile = [];
            end
            
            [allVals] = static_optimization_algs.BayesOptGPML.optimize(optimizerInput.functionBounds, optimizerInput.nbInitSamplesBO, optimizerInput.nbEvalsBO, optimizerInput.maxNbSamplesBO, ...
                func, optimizerInput.kernelType, optimizerInput.featureFunction, optimizerInput.yCenteringType, optimizerInput.gpHyperOption, optimizerInput.videoFile);
        end
        
        %% main function
        function [allVals] = optimize(functionBounds, nbInitEvals, nbTotalEvals, maxNbSamples, fun, kernelType, featureFunction, yCenteringType, gpHyperOption, videoFile)
            usedSamples = [];
            usedVals = [];
            allVals = [];
            dim = size(functionBounds, 2);
            lb = functionBounds(1, :);
            ub = functionBounds(2, :);
            %%% initializing GP
            meanfunc = [];
            covfunc = @covSEiso;
            likfunc = @likGauss;
            hyp = struct('mean', [], 'cov', [0 0], 'lik', -1);
            
            %%% animation
            if(~isempty(videoFile))
                open(videoFile);
            end
            
            for iter = 1:nbTotalEvals
                x = usedSamples;
                y = usedVals;
                hyp_param = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, y);
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
                    [newSample] = static_optimization_algs.BayesOptGPML.optimizeAcquisition...
                        (x, y, functionBounds, hyp_param, @infGaussLik, meanfunc, covfunc, likfunc);
                end
                
                newVal = fun.eval(newSample);
                allVals = [allVals; newVal];
                
                usedSamples = [usedSamples; newSample];
                usedVals = [usedVals; newVal];
                if(length(usedVals) > maxNbSamples)
                    usedSamples = usedSamples(end-maxNbSamples+1:end, :);
                    usedVals = usedVals(end-maxNbSamples+1:end);
                end
                
                %%% current distrib plotting
                if(~isempty(videoFile))
                    frame = static_optimization_algs.BayesOptGPML.getFrame(fun, x, y, hyp_param, @infGaussLik, meanfunc, covfunc, likfunc, iter);
                    writeVideo(videoFile, frame);
                end
            end
            
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        function [newSample] = optimizeAcquisition(x, y, functionBounds, hyp_param, infGaussLik, meanfunc, covfunc, likfunc)
            % Setting the function
            func = @(xt) gp(hyp_param, infGaussLik, meanfunc, covfunc, likfunc, x, y, xt);
            
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
            optOption.maxeval = 200;
%             [~, bestIndex] = max(y);
%             params0 = x(bestIndex, :);
            params0 = (functionBounds(1, :) + functionBounds(2, :)) * .5;
            [newSample, ~, retCode] = nlopt_optimize(optOption, params0);
            if(retCode < 0)
                warning('Error in NL_OPT optimization');
            end
            optOption.algorithm = NLOPT_LN_BOBYQA;
            [newSample, ~, retCode] = nlopt_optimize(optOption, newSample);
            if(retCode < 0)
                warning('Error in NL_OPT optimization');
            end
        end
        
        function frame = getFrame(fun, samples, vals, hyp_param, infGaussLik, meanfunc, covfunc, likfunc, iter)
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
                static_optimization_algs.GP.plotGPML(hyp_param, infGaussLik, meanfunc, covfunc, likfunc, samples, vals, center, length);
            end
            text(-5, -5, ['iteration = ' num2str(iter)]);
            frame = getframe(gcf,[0 0 560 420]);
        end
    end
end
