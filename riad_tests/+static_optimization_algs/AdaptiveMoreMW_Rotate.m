classdef AdaptiveMoreMW_Rotate
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['AdaptiveMoreMW_Rotate' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbInitSamples) '_' ...
                num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.initVar)];
        end
        %%
        function [perf, kls] = optimizeStruct(optimizerInput)
            [perf, kls] = static_optimization_algs.AdaptiveMoreMW_Rotate.optimize(optimizerInput.initDistrib, optimizerInput.minEpsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbInitSamples, optimizerInput.nbIter, optimizerInput.maxIterReuse,...
                optimizerInput.confidenceMW, optimizerInput.regularization, optimizerInput.fun, optimizerInput.videoFile);
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
                %%% Importance weights
                logInitProbas = currentDistrib.getLogProbas(allSamples);
                logImportanceProbas = static_optimization_algs.NonParamImportanceSampling.getLogProbas(allSamples, currentDistrib.getPrecision, logInitProbas);
                
                %%% Learn Concave Quadratic Reward Model
                sampleWeights = exp(logInitProbas-logImportanceProbas);
                [~, quadModel] = static_optimization_algs.Quadratic.learnConcaveQuadModel(allSamples, allVals, regularization, sampleWeights); 


                %%% I- Update policy, move mean
                maxEpsiKL = 1;
                [etaStar, isBadModel] = static_optimization_algs.AdaptiveMore.lineSearch(allSamples, allVals, logImportanceProbas, currentDistrib, quadModel, minEpsiKL, maxEpsiKL, 100, confidence, 0);
%                 if(isBadModel)
%                     disp('bad quadratic model, switching to IProjection for this time around');
%                     quadModel = static_optimization_algs.AdaptiveMoreIProj.learnConcaveQuadModelWithIProj(currentDistrib, maxEpsiKL, allSamples, allVals, logImportanceProbas, 1);
%                     [etaStar, isBadModel] = static_optimization_algs.AdaptiveMore.lineSearch(allSamples, allVals, logImportanceProbas, currentDistrib, quadModel, minEpsiKL, maxEpsiKL, 100, confidence, 0);
%                     if(isBadModel)
%                         disp('neither the regression nor the iprojection worked.');
%                     end
%                 end
                movedMeanDistrib = static_optimization_algs.More.updatePolicy(currentDistrib, quadModel, etaStar, 1);
                                
                %%% II- Update policy, rotate covariance
                desiredEntropy = movedMeanDistrib.getEntropy;
                logInitProbas = movedMeanDistrib.getLogProbas(allSamples);
                sampleWeights = exp(logInitProbas - logImportanceProbas);
                [~, quadModel] = static_optimization_algs.Quadratic.learnCentredConcaveQuadModel(allSamples, allVals, movedMeanDistrib.getMu, regularization, sampleWeights);
                etaStar = static_optimization_algs.AdaptiveMore.lineSearch(allSamples, allVals, logImportanceProbas, movedMeanDistrib, quadModel, 1e-6, maxEpsiKL, 50, confidence, 0, desiredEntropy);
                nextPol = static_optimization_algs.More.updatePolicy(movedMeanDistrib, quadModel, etaStar, 1, desiredEntropy);
                
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
                kls(iter) = currentDistrib.getKL(nextPol);
                currentDistrib = nextPol;
                
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
    end
end
