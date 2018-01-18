classdef AdaptiveMoreWithModelSwitch
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['AdaptiveMoreWithModelSwitch' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbInitSamples) '_' ...
                num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.initVar)];
        end
        %%
        function [perf, kls] = optimizeStruct(optimizerInput)
            [perf, kls] = static_optimization_algs.AdaptiveMoreWithModelSwitch.optimize(optimizerInput.initDistrib, optimizerInput.minEpsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbInitSamples, optimizerInput.nbIter, optimizerInput.maxIterReuse, optimizerInput.confidenceMW, optimizerInput.regularization, optimizerInput.fun, optimizerInput.videoFile);
        end
        
        %% main function
        function [perf, kls] = optimize(initDistrib, minEpsiKL, nbSamplesPerIter, nbInitSamples, nbIter, maxIterReuse, confidence, regularization, fun, videoFile)
            perf = zeros(nbIter, 2);
            kls = zeros(nbIter, 1);
            nbSamplesForEval = 100;
            dim = initDistrib.getDim;
            currentDistrib = initDistrib;
            allSamples = [];
            allVals = [];
            
            optimizerConstants.minEpsiKL = minEpsiKL;
            optimizerConstants.maxEpsiKL = minEpsiKL * 100;
            optimizerConstants.nbLineSearch = 100;
            optimizerConstants.probaError = confidence;
            optimizerConstants.regularization = regularization;
            
            %%% animation
            if(~isempty(videoFile))
                open(videoFile);
            end
            
            for iter = 1:nbIter
                %%% evaluation
                numberOfSamplesEval = max(nbSamplesPerIter, nbSamplesForEval);
                if(iter == 1)
                    numberOfSamplesEval = max(numberOfSamplesEval, nbInitSamples);
                end
                newSamples = currentDistrib.getSamples(numberOfSamplesEval);
                vals = fun.eval(newSamples);
                newPerf = mean(vals);
                
                %%% samples
                if(iter == 1)
                    numberOfSamplesUpdate = nbInitSamples;
                    perf(iter, :) = [0 newPerf];
                else
                    numberOfSamplesUpdate = nbSamplesPerIter;
                    perf(iter, :) = [((iter-2) * nbSamplesPerIter + nbInitSamples) newPerf];
                end
                
                newSamples = newSamples(1:numberOfSamplesUpdate, :);
                vals = vals(1:numberOfSamplesUpdate);
                allSamples = [allSamples; newSamples];
                allVals = [allVals; vals];
                
                %%%plotting
                if(~isempty(videoFile))
                    frame = static_optimization_algs.AdaptiveMore.getFrame(fun, currentDistrib, allSamples, iter);
                    writeVideo(videoFile, frame);
                end
                
                %%% update
                [currentDistrib, kl] = static_optimization_algs.AdaptiveMoreWithModelSwitch.updateSearchDistribution(currentDistrib, allSamples, allVals, optimizerConstants);
                
                %%% delete old samples
                if(iter > maxIterReuse)
                    if(iter == maxIterReuse)
                        allSamples = allSamples(nbInitSamples+1:end, :);
                        allVals = allVals(nbInitSamples+1:end, :);
                    else
                        allSamples = allSamples(nbSamplesPerIter+1:end, :);
                        allVals = allVals(nbSamplesPerIter+1:end, :);
                    end
                end
                kls(iter) = kl;
                
                %!debug                
                if(currentDistrib.getEntropy < -2)
                    perf(iter+1:end, :) = [((iter-1:nbIter-2) * nbSamplesPerIter + nbInitSamples)' (perf(iter, 2) * ones(nbIter-iter, 1))];
                    break;
                end
            end
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        function [nextDistrib, kl, concaveQuadModel] = updateSearchDistribution(currentDistrib, allSamples, allVals, optimizerConstants, guessedConcaveQuadModel)
            if(~exist('guessedConcaveQuadModel', 'var'))
                guessedConcaveQuadModel = [];
            end
            
            minEpsiKL = optimizerConstants.minEpsiKL;
            maxEpsiKL = optimizerConstants.maxEpsiKL;
            nbLineSearch = optimizerConstants.nbLineSearch;
            regularization = optimizerConstants.regularization;
            probaError = optimizerConstants.probaError;
            
            %%% Importance weights
            logInitProbas = currentDistrib.getLogProbas(allSamples);
            logImportanceProbas = static_optimization_algs.NonParamImportanceSampling.getLogProbas(allSamples, currentDistrib.getPrecision, logInitProbas);
            
            %%% Learn Concave Quadratic Reward Model
            allVals = allVals / max(abs(allVals));
            sampleWeights = exp(logInitProbas - logImportanceProbas - max(logInitProbas - logImportanceProbas));
            sampleWeights = static_optimization_algs.Quadratic.normalizeWeights(sampleWeights);
            [~, concaveQuadModel] = static_optimization_algs.Quadratic.learnQuadModel(allSamples, allVals, regularization, sampleWeights, .999);
            
            
            %%% I- Update policy, move mean
            [etaStar, ~] = static_optimization_algs.AdaptiveMore.lineSearch(allSamples, allVals, logImportanceProbas, currentDistrib, concaveQuadModel, minEpsiKL, maxEpsiKL, nbLineSearch, probaError, 0);            
            [nextDistrib, ~, kl] = static_optimization_algs.More.updatePolicy(currentDistrib, concaveQuadModel, etaStar, 1);            
        end
        
    end
end
