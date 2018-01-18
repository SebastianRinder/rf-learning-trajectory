classdef AdaptiveMoreStudentEntropy
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['AdaptiveMoreStudentEntropy' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbInitSamples) '_' ...
                num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.initVar)];
        end
        %%
        function [perf, kls] = optimizeStruct(optimizerInput)
            avrRewardLB = @(vals, logImportanceProbas, logProbas) static_optimization_algs.ConcentrationInequalities.studentLowerBound(vals, logImportanceProbas, logProbas, optimizerInput.confidenceStudent);
            [perf, kls] = static_optimization_algs.AdaptiveMoreStudentEntropy.optimize(optimizerInput.initDistrib, optimizerInput.minEpsiKL, optimizerInput.entropyReduction, ...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbInitSamples, optimizerInput.nbIter, optimizerInput.maxIterReuse, avrRewardLB, optimizerInput.regularization, optimizerInput.fun, optimizerInput.videoFile);
        end
        
        %% main function
        function [perf, kls] = optimize(initDistrib, minEpsiKL, entropyReduction, nbSamplesPerIter, nbInitSamples, nbIter, maxIterReuse, avrRewardLB, regularization, fun, videoFile)
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
            optimizerConstants.avrRewardLB = avrRewardLB;
            optimizerConstants.entropyReduction = entropyReduction;
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
                
                %%% update
                %%% Trimming very low weight samples
                logInitProbas = currentDistrib.getLogProbas(allSamples);
                [sWeights, is] = sort(exp(logInitProbas), 'descend');
                firstId = find(cumsum(sWeights / sum(sWeights)) > .999, 1, 'first');
                if(isempty(firstId))
                    firstId = length(sWeights) + 1;
                end
                is = is(1:firstId-1);
                allSamplesTrim = allSamples(is, :);
                allValsTrim = allVals(is, :);
                
                %%% current distrib plotting
                if(~isempty(videoFile))
                    frame = static_optimization_algs.AdaptiveMore.getFrame(fun, currentDistrib, allSamplesTrim, iter);
                    writeVideo(videoFile, frame);
                end

                %%% get new distrib
                [currentDistrib, kl] = static_optimization_algs.AdaptiveMoreStudentEntropy.updateSearchDistribution...
                    (currentDistrib, allSamplesTrim, allValsTrim, optimizerConstants);
                
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
            %%% unpacking constants required by this optimizer 
            minEpsiKL = optimizerConstants.minEpsiKL;
            maxEpsiKL = optimizerConstants.maxEpsiKL;
            nbLineSearch = optimizerConstants.nbLineSearch;
            avrRewardLB = optimizerConstants.avrRewardLB;
            entropyReduction = optimizerConstants.entropyReduction;
            regularization = optimizerConstants.regularization;

            %%% Importance weights
            logInitProbas = currentDistrib.getLogProbas(allSamples);
            logImportanceProbas = static_optimization_algs.NonParamImportanceSampling.getLogProbas(allSamples, currentDistrib.getPrecision, logInitProbas);
            
            %%% Learn Concave Quadratic Reward Model
            allVals = allVals / max(abs(allVals));
            sampleWeights = exp(logInitProbas - logImportanceProbas - max(logInitProbas - logImportanceProbas));
            sampleWeights = static_optimization_algs.Quadratic.normalizeWeights(sampleWeights);
            [~, concaveQuadModel] = static_optimization_algs.Quadratic.learnQuadModel(allSamples, allVals, regularization, sampleWeights, .999);

            %%% Update policy
            desiredEntropy = currentDistrib.getEntropy - entropyReduction;
            [etaStar, ~] = static_optimization_algs.AdaptiveMoreStudent.lineSearchStudent(allSamples, allVals, logImportanceProbas, currentDistrib,...
                concaveQuadModel, minEpsiKL, maxEpsiKL, nbLineSearch, avrRewardLB, 0, desiredEntropy);
            [nextDistrib, ~, kl] = static_optimization_algs.More.updatePolicy(currentDistrib, concaveQuadModel, etaStar, 1, desiredEntropy);
        end
        
    end
end
