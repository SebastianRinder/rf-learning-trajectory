classdef RepsBanditsIProjection
    methods (Static)
        function sign = getSignature(optimizerInput)
            sign = ['RepsBanditsIProjection' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.initVar)];
        end
        
        function [perf] = optimizeStruct(optimizerInput)
            perf = static_optimization_algs.RepsBanditsIProjection.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbIter, optimizerInput.fun, optimizerInput.videoFile);
        end
        
        function [perf] = optimize(initDistrib, epsiKL, nbSamplesPerIter, nbIter, fun, videoFile)
            perf = zeros(nbIter, 2);
            nbSamplesForEval = 100;
            currentDistrib = initDistrib;
            dim = currentDistrib.getDim;
            
            %%% IProjection optimization option
            optOptim =  optimset('GradObj', 'on', 'DerivativeCheck', 'off', 'Display', 'off');
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
                params0 = [initMu'; initPrecision(tril(ones(dim)) == 1)];
                logImportanceProbas = currentDistrib.getLogProbas(newSamples);
%                 params = fmincon(@(params) static_optimization_algs.IProjection.eval(params, dim, initPrecision, initMu, newSamples, vals, importanceProbas, eta), params0, [], [], [], [], [], [], [], optOptim);
                params = fminunc(@(params) static_optimization_algs.IProjection.eval(params, initPrecision, initMu, newSamples, vals, logImportanceProbas, eta), params0, optOptim);
                
                muOptim = params(1:dim);
                covOptim = zeros(dim);
                covOptim(tril(ones(dim)) == 1) = params(dim+1:end);
                covOptim = covOptim + covOptim' - diag(diag(covOptim));
                covOptim = inv(covOptim);
                currentDistrib = static_optimization_algs.Normal(muOptim, covOptim);
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
            plot(samples(:, 1), samples(:, 2), '*r');
            policy.plot();
            frame = getframe(gcf,[0 0 560 420]);
        end
        
    end
end
