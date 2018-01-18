classdef RepsBandits
    methods (Static)
        function sign = getSignature(optimizerInput)
            sign = ['RepsBandits' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.initVar)];
        end
        function [perf, kls] = optimizeStruct(optimizerInput)
            [perf, kls] = static_optimization_algs.RepsBandits.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbIter, optimizerInput.maxIterReuse, optimizerInput.fun, optimizerInput.videoFile);
        end
        
        function [perf, kls] = optimize(initDistrib, epsiKL, nbSamplesPerIter, nbIter, maxIterReuse, fun, videoFile)
            perf = zeros(nbIter, 2);
            nbSamplesForEval = 100;
            currentDistrib = initDistrib;
            kls = zeros(nbIter, 1);
            allSamples = [];
            allVals = [];

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
                allSamples = [allSamples; newSamples];
                allVals = [allVals; vals];
                
                %% current distrib plotting
                if(~isempty(videoFile))
                    currPlot = figure;
                    hold on;
                    fun.plot();
                    plot(allSamples(:, 1), allSamples(:, 2), '*r');
                    currentDistrib.plot(10);
                    writeVideo(videoFile, getframe(gcf,[0 0 560 420]));
                    close(currPlot);
                end
                
                %% update
                importanceSamplingW = ones(size(allVals));
                if(length(allVals) ~= nbSamplesPerIter)
                    %importance weights
                    logInitProbas = currentDistrib.getLogProbas(allSamples);
                    initPrecision = currentDistrib.getPrecision();
                    logImportanceProbas = static_optimization_algs.NonParamImportanceSampling.getLogProbas(allSamples, initPrecision, logInitProbas);
                    importanceSamplingW = exp(logInitProbas - logImportanceProbas);
                end
                
                [eta, sampleWeights] = static_optimization_algs.RepsBandits.optimizeDual(allVals, epsiKL, importanceSamplingW);
                oldDistrib = currentDistrib;
                currentDistrib = currentDistrib.wmle(allSamples, exp(allVals / eta) .* sampleWeights);
                kls(iter) = currentDistrib.getKL(oldDistrib);
                
                if(iter > maxIterReuse)
                    allSamples = allSamples(nbSamplesPerIter+1:end, :);
                    allVals = allVals(nbSamplesPerIter+1:end, :);
                end
            end
            
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        function [eta, sampleWeights] = optimizeDual(rewards, epsiKL, sampleWeights)
            if(~exist('sampleWeights', 'var'))
                sampleWeights = ones(size(rewards));
            end
            sampleWeights = sampleWeights / sum(sampleWeights);
            
            infEta = 1e-16;
            %             persistent initEta;
            %             if(isempty(initEta))
            %                 initEta = 1;
            %             end
            options =  optimset('GradObj','on', 'DerivativeCheck', 'off', 'Display', 'off');
            eta = fmincon(@(eta) static_optimization_algs.RepsBandits.dual(rewards, eta, epsiKL, sampleWeights), 1, [], [], [], [], infEta, inf,[],options);
            %             initEta = eta;
        end
        
        function [g, dg] = dual(rewards, eta, epsiKL, sampleWeights)
            sumExpRewards = sum(exp(rewards / eta) .* sampleWeights);
            sumExpRR = sum(exp(rewards / eta) .* rewards .* sampleWeights);
            g = eta * epsiKL + eta * log (sumExpRewards);
            dg = epsiKL + log (sumExpRewards) - sumExpRR / (eta * sumExpRewards);
        end
    end
end
