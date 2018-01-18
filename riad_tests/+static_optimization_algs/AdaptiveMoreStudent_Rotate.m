classdef AdaptiveMoreStudent_Rotate
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['AdaptiveMoreStudent_Rotate' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbInitSamples) '_' ...
                num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.initVar)];
        end
        %%
        function [perf, kls] = optimizeStruct(optimizerInput)
            avrRewardLB = @(vals, logImportanceProbas, logProbas) static_optimization_algs.ConcentrationInequalities.studentLowerBound(vals, logImportanceProbas, logProbas, optimizerInput.confidenceStudent);
            [perf, kls] = static_optimization_algs.AdaptiveMoreStudent_Rotate.optimize(optimizerInput.initDistrib, optimizerInput.minEpsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbInitSamples, optimizerInput.nbIter, optimizerInput.maxIterReuse, avrRewardLB, optimizerInput.regularization, optimizerInput.fun, optimizerInput.videoFile);
        end
        
        %% main function
        function [perf, kls] = optimize(initDistrib, minEpsiKL, nbSamplesPerIter, nbInitSamples, nbIter, maxIterReuse, avrRewardLB, regularization, fun, videoFile)
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
                
                
                %%% update
                %%% Importance weights
                logInitProbas = currentDistrib.getLogProbas(allSamples);
                [sWeights, is] = sort(exp(logInitProbas), 'descend');
                firstId = find(cumsum(sWeights / sum(sWeights)) > 2, 1, 'first');
                if(isempty(firstId))
                    firstId = length(sWeights) + 1;
                end
                is = is(1:firstId-1);
                logInitProbas = logInitProbas(is);
                allSamplesTrim = allSamples(is, :);
                allValsTrim = allVals(is, :);
                logImportanceProbas = static_optimization_algs.NonParamImportanceSampling.getLogProbas(allSamplesTrim, currentDistrib.getPrecision, logInitProbas);
                
                
                
                %%% Learn Concave Quadratic Reward Model
                sampleWeights = exp(logInitProbas-logImportanceProbas);
                [~, quadModel] = static_optimization_algs.Quadratic.learnConcaveQuadModel(allSamplesTrim, allValsTrim, regularization, sampleWeights);
               
                %%% current distrib plotting
                if(~isempty(videoFile))
                    frame = static_optimization_algs.AdaptiveMore.getFrame(fun, currentDistrib, allSamplesTrim, iter);
                    writeVideo(videoFile, frame);
                end
               
                
                %%% I- Update policy: mean
                maxEpsiKL = 1;

                etaStar = static_optimization_algs.AdaptiveMoreStudent.lineSearchStudent(allSamplesTrim, allValsTrim, logImportanceProbas, currentDistrib, quadModel, minEpsiKL, maxEpsiKL, 100, avrRewardLB, 0);
%                 if(isBadModel)
%                     disp('bad quadratic model, switching to IProjection for this time around');
%                     quadModel = static_optimization_algs.AdaptiveMoreIProj.learnConcaveQuadModelWithIProj(currentDistrib, maxEpsiKL, allSamplesTrim, allValsTrim, logImportanceProbas, 1);
%                     [etaStar, isBadModel] = static_optimization_algs.AdaptiveMoreStudentWithModelSwitch.lineSearchStudent(allSamplesTrim, allValsTrim, logImportanceProbas, currentDistrib, quadModel, minEpsiKL, maxEpsiKL, 100, avrRewardLB, 0);
%                     if(isBadModel)
%                         disp('neither the regression nor the iprojection worked.');
%                     end
%                 end
                movedMeanDistrib = static_optimization_algs.More.updatePolicy(currentDistrib, quadModel, etaStar, 1);

                %%% II- Update policy, rotate covariance   
                desiredEntrop = movedMeanDistrib.getEntropy;
                logInitProbas = movedMeanDistrib.getLogProbas(allSamplesTrim);
                sampleWeights = exp(logInitProbas - logImportanceProbas);
                [~, quadModel] = static_optimization_algs.Quadratic.learnCentredConcaveQuadModel(allSamplesTrim, allValsTrim, movedMeanDistrib.getMu, regularization, sampleWeights);
                [etaStar, ~] = static_optimization_algs.AdaptiveMoreStudent.lineSearchStudent(allSamplesTrim, allValsTrim, logImportanceProbas, movedMeanDistrib, quadModel, minEpsiKL, maxEpsiKL, 100, avrRewardLB, 0, desiredEntrop);
                nextPol = static_optimization_algs.More.updatePolicy(movedMeanDistrib, quadModel, etaStar, 1, desiredEntrop);
                
                
                %%% delete old samples
                if(iter > maxIterReuse)
                    if(iter == maxIterReuse) % first time, delete nbInitSamples
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
        
        function [etaStar, isBadModel] = lineSearchStudent(samples, vals, logImportanceProbas, initPolicy, quadFunction, minKL, maxKL, nbSearch, avrRewardLBFun, verbose)
            %%%logging for verbose mode
            if(~exist('verbose', 'var'))
                verbose = 0;
            end
            if(verbose)
                loggingLabel = {'etas', 'kls', 'entrops', 'modelRewardImprovement', 'avrReward', 'avrRewardLb'};
                logging = zeros(3, nbSearch);
            end            
            isBadModel = false;
            
            %%% linesearch eta -> rewardLb
            rewardLb = zeros(1, nbSearch);
            etas = zeros(1, nbSearch);
            avrRewards = zeros(1, nbSearch);            
            logEtaSup = log(static_optimization_algs.More.optimizeDual(initPolicy, quadFunction, minKL));
            logEtaInf = log(static_optimization_algs.More.optimizeDual(initPolicy, quadFunction, maxKL));
            
            for k = 1:nbSearch
                etas(k) = exp(logEtaInf + (k-1) * (logEtaSup - logEtaInf) / (nbSearch - 1));
                [newPol, rImprove, kl, entrop]  = static_optimization_algs.More.updatePolicy(initPolicy, quadFunction, etas(k), 0);
                logProbas = newPol.getLogProbasP(samples);
                [rewardLb(k), avrRewards(k)] = avrRewardLBFun(vals, logImportanceProbas, logProbas);
                
                %%% logging
                if(verbose)
                    logging(:, k) = [kl; entrop; rImprove];
                end
            end
            
            logInitProbas = initPolicy.getLogProbasP(samples);
            avrInitReward = mean(vals .* exp(logInitProbas - logImportanceProbas));
            %all avr rewards are smaller than the init one. Maybe the quad
            % model is bad?
            if(~any(avrRewards > avrInitReward))
                isBadModel = true;
            end
            
            
            [~, i] = max(rewardLb);
            etaStar = etas(i);
            
            %%% plotting
            if(verbose)
                logging = [etas; logging; avrRewards; rewardLb];
                for k = 1:length(loggingLabel)
                    figure;
                    plot(logging(k, :));
                    title(loggingLabel{k});
                end
            end
        end
    end
end
