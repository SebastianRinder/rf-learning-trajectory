classdef RepsBanditsIProjectionChol
    methods (Static)
        function sign = getSignature(optimizerInput)
            sign = ['RepsBanditsIProjectionChol' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.initVar)];
        end
        
        function [perf] = optimizeStruct(optimizerInput)
            perf = static_optimization_algs.RepsBanditsIProjectionChol.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbIter, optimizerInput.fun, optimizerInput.videoFile);
        end
        
        function [perf] = optimize(initDistrib, epsiKL, nbSamplesPerIter, nbIter, fun, videoFile)
            perf = zeros(nbIter, 2);
            nbSamplesForEval = 100;
            currentDistrib = initDistrib;
            dim = currentDistrib.getDim;
            
            %%% IProjection optimization option
            optOptim =  optimset('GradObj','on', 'DerivativeCheck', 'off', 'Display', 'off');
%             optOptim.TolFun = .01;
            
            %% animation
            if(~isempty(videoFile))
                open(videoFile);
            end
            
            for iter = 1:nbIter
                %% evaluation
                newSamples = currentDistrib.getSamples(max(nbSamplesPerIter, nbSamplesForEval));
                vals = fun.eval(newSamples);
                newPerf = mean(vals);
                perf(iter, :) = [((iter-1) * nbSamplesPerIter) newPerf];
                
                %% samples
                newSamples = newSamples(1:nbSamplesPerIter, :);
                vals = vals(1:nbSamplesPerIter);
                
                %% current distrib plotting
                if(~isempty(videoFile))
                    frame = static_optimization_algs.RepsBanditsIProjection.getFrame(fun, currentDistrib, newSamples);
                    writeVideo(videoFile, frame);
                end
                %% update
                [eta, sampleWeights] = static_optimization_algs.RepsBandits.optimizeDual(vals, epsiKL);
                
                %%% optimizing the I-Projection
                initPrecision = currentDistrib.getPrecision;
                initMu = currentDistrib.getMu;
                initChol = currentDistrib.getCholP;
                params0 = [initMu'; initChol(triu(ones(dim)) == 1)];
                lbParams = -inf * ones(size(params0));
                lbParams(dim+1:end) = 0;
                logImportanceProbas = currentDistrib.getLogProbas(newSamples);
                params = fmincon(@(params) static_optimization_algs.IProjection.eval_chol(params, initPrecision, initMu, newSamples, vals, logImportanceProbas, eta), params0, [], [], [], [], lbParams, [], [], optOptim);
                
                muOptim = params(1:dim);
                cholOptim = zeros(dim);
                cholOptim(triu(ones(dim)) == 1) = params(dim+1:end);
                currentDistrib = static_optimization_algs.Normal(muOptim, cholOptim, true);
            end
            
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        function frame = getFrame(fun, policy, samples)
            persistent currPlot;
            if(isempty(currPlot))
                currPlot = figure;
            end
            hold off;
            fun.plot();
            hold on;
            policy.plot();
            plot(samples(:, 1), samples(:, 2), '*r');
            frame = getframe(gcf,[0 0 560 420]);
        end
        
    end
end
