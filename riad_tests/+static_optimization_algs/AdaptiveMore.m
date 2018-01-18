classdef AdaptiveMore
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['AdaptiveMore' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbInitSamples) '_' ...
                num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.initVar)];
        end
        %%
        function [perf, kls] = optimizeStruct(optimizerInput)
            [perf, kls] = static_optimization_algs.AdaptiveMore.optimize(optimizerInput.initDistrib, optimizerInput.minEpsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbInitSamples, optimizerInput.nbIter, optimizerInput.maxIterReuse, optimizerInput.confidence, optimizerInput.regularization, optimizerInput.fun, optimizerInput.videoFile);
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
                
                %%% current distrib plotting
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
%                 if(~exist('quadModel', 'var'))
%                     [~, quadModel, phiSamples] = static_optimization_algs.Quadratic.learnConcaveQuadModel(allSamples, allVals, regularization, sampleWeights);
%                 else 
%                     phiSamples = [phiSamples; static_optimization_algs.Quadratic.getQuadFeatures(newSamples, 1)];
%                     [~, quadModel] = static_optimization_algs.Quadratic.optimizeThetas(quadModel, dim, phiSamples, allVals, regularization, sampleWeights);
%                 end
                
                [~, quadModel] = static_optimization_algs.Quadratic.learnConcaveQuadModel(allSamples, allVals, regularization, sampleWeights); 


                %%% Update policy
                maxEpsiKL = 1;
                etaStar = static_optimization_algs.AdaptiveMore.lineSearch(allSamples, allVals, logImportanceProbas, currentDistrib, quadModel, minEpsiKL, maxEpsiKL, 100, confidence, 0);
                [currentDistrib, ~, kl] = static_optimization_algs.More.updatePolicy(currentDistrib, quadModel, etaStar, 1);
                
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
        
        %%
        function [etaStar, isBadModel] = lineSearch(samples, vals, logImportanceProbas, initPolicy, quadFunction, minKL, maxKL, nbSearch, conf, verbose, desiredEntropy)
            %%%logging for verbose mode
            if(~exist('verbose', 'var'))
                verbose = 0;
            end
            if(~exist('desiredEntropy', 'var'))
                desiredEntropy = [];
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
            trimedRatio = .9;
            minSamples = 8;
            
            logEtaSup = log(static_optimization_algs.More.optimizeDual(initPolicy, quadFunction, minKL, desiredEntropy));
            logEtaInf = log(static_optimization_algs.More.optimizeDual(initPolicy, quadFunction, maxKL, desiredEntropy));

            logInitProbas = initPolicy.getLogProbasP(samples);
            [trimedInitRewards, tooLessInit] = static_optimization_algs.ConcentrationInequalities.trimImprobable(logInitProbas, logImportanceProbas, vals, trimedRatio, minSamples);

            for k = 1:nbSearch
                etas(k) = exp(logEtaInf + (k-1) * (logEtaSup - logEtaInf) / (nbSearch - 1));
                [newPol, rImprove, kl, entrop]  = static_optimization_algs.More.updatePolicy(initPolicy, quadFunction, etas(k), 0, desiredEntropy);
                logProbas = newPol.getLogProbasP(samples);
                %                 rewardLb(k) = static_optimization_algs.ConcentrationInequalities.andersonLowerBound(vals, logImportanceProbas, logProbas, conf);
                [trimedRewards, tooLess] = static_optimization_algs.ConcentrationInequalities.trimImprobable(logProbas, logImportanceProbas, vals, trimedRatio, minSamples);
                if(~tooLess && ~tooLessInit)
                    rewardLb(k) = ranksum(trimedRewards, trimedInitRewards);
                else
                    rewardLb(k) = 1;
                end
                avrRewards(k) = mean(vals .* exp(logProbas - logImportanceProbas));
                
                %%% logging
                if(verbose)
                    logging(:, k) = [kl; entrop; rImprove];
                end
            end
            avrInitReward = mean(vals .* exp(logInitProbas - logImportanceProbas));
            %all avr rewards are smaller than the init one. Maybe the quad 
            % model is bad?
            if(~any(avrRewards > avrInitReward))
                isBadModel = true;
            end
                
            acceptedMoves = find(rewardLb < conf);
            if(~isempty(acceptedMoves))
                [maxVal, i] = max(avrRewards(acceptedMoves));
                if(maxVal > avrInitReward)
                    etaStar = etas(acceptedMoves(i));
                    if(verbose)
                        fprintf('Selected eta %g with avrReward %g and p-value %g\n', etaStar, maxVal, rewardLb(i));
                    end
                else
                    if(verbose)
                        disp('Only significantly different move decreases the average reward');
                    end
                    etaStar = exp(logEtaSup);
                end
            else
                if(verbose)
                    disp('No significantly different move');                    
                end
                etaStar = exp(logEtaSup);
            end
            
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
        
        function frame = getFrame(fun, policy, samples, iter, quadModel)
            persistent currPlot;
            if(isempty(currPlot))
                currPlot = figure;
            end
            hold off;
            fun.plot();
            hold on;
            plot(samples(:, 1), samples(:, 2), '*r');
            policy.plot();
            if(exist('quadModel', 'var'))
                static_optimization_algs.Quadratic.plotQuad(quadModel, [0 0], [20 20]);
            end
            text(-5, -5, ['iteration = ' num2str(iter)]);
            frame = getframe(gcf,[0 0 560 420]);
        end
    
    end
end
