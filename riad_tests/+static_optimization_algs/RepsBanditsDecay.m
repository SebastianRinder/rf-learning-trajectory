classdef RepsBanditsDecay
    methods (Static)
        function sign = getSignature(optimizerInput)
            sign = ['RepsGaussSampleDecay' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.sampleDecay) '_' ...
                num2str(optimizerInput.initVar)];
        end
        
        function [perf] = optimizeStruct(optimizerInput)
            perf = static_optimization_algs.RepsBanditsDecay.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbIter, optimizerInput.sampleDecay,...
                optimizerInput.fun, optimizerInput.videoFile);
        end
        
        function [perf] = optimize(initDistrib, epsiKL, nbSamplesPerIter, nbIter, sampleDecay, fun, videoFile)
            perf = zeros(nbIter, 2);
            nbSamplesForEval = 100;
            currentDistrib = initDistrib;
            sampleWeights = [];
            samples = [];
            vals = [];
            
            %% animation
            open(videoFile);
            
            for iter = 1:nbIter
                %% new samples
                newSamples = currentDistrib.getSamples(nbSamplesPerIter);
                vals = [vals; fun.eval(newSamples)];
                sampleWeights = [sampleWeights; ones(nbSamplesPerIter, 1)];
                samples = [samples; newSamples];
                
                %% evaluation
                newPerf = mean(fun.eval(currentDistrib.getSamples(nbSamplesForEval)));
                perf(iter, :) = [((iter-1) * nbSamplesPerIter) newPerf];
                
                %% current distrib plotting
                frame = static_optimization_algs.RepsBanditsDecay.getFrame(fun, currentDistrib, newSamples);
                writeVideo(videoFile, frame);
                
                %% update
                [eta, sampleWeights] = static_optimization_algs.RepsBandits.optimizeDual(vals, epsiKL, sampleWeights);
                currentDistrib = currentDistrib.wmle(samples, exp(vals / eta) .* sampleWeights);
                sampleWeights = sampleWeights * sampleDecay;
            end            
            close(videoFile);
        end
        
        function frame = getFrame(fun, policy, samples)
            persistent currPlot;
            if(isempty(currPlot))
                currPlot = figure;
            end
            hold off;
            fun.plot();
            hold on;
            policy.plot(10);
            plot(samples(:, 1), samples(:, 2), '*r');
            frame = getframe(gcf,[0 0 560 420]);
        end
    end
end
