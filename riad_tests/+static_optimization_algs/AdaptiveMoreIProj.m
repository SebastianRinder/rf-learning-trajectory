classdef AdaptiveMoreIProj
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['AdaptiveMoreIProj' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbInitSamples) '_' ...
                num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.initVar)];
        end
        %%
        function [perf, kls] = optimizeStruct(optimizerInput)
            [perf, kls] = static_optimization_algs.AdaptiveMoreIProj.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL,...
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
                    frame = static_optimization_algs.RepsBanditsIProjection.getFrame(fun, currentDistrib, allSamples);
                    writeVideo(videoFile, frame);
                end
                
                %%% update
                %%% Importance weights
                logInitProbas = currentDistrib.getLogProbas(allSamples);
                logImportanceProbas = static_optimization_algs.NonParamImportanceSampling.getLogProbas(allSamples, currentDistrib.getPrecision, logInitProbas);
                
                %%% Learn Concave Quadratic Reward Model (fitted with IProjection)
                maxEpsiKL = 1;
                if(~exist('etaModelLearning', 'var'))
                    etaModelLearning = 1;
                end
                [quadModel, etaModelLearning] = static_optimization_algs.AdaptiveMoreIProj.learnConcaveQuadModelWithIProj(currentDistrib, maxEpsiKL, allSamples, allVals, logImportanceProbas, etaModelLearning);

                %%% Update policy
                etaStar = static_optimization_algs.AdaptiveMore.lineSearch(allSamples, allVals, logImportanceProbas, currentDistrib, quadModel, minEpsiKL, maxEpsiKL, 100, confidence, 1);
                [currentDistrib, ~, kl] = static_optimization_algs.More.updatePolicy(currentDistrib, quadModel, etaStar, 1);             
                
                %%% delete old samples
                if(iter > maxIterReuse)
                    if(length(allVals) ~= maxIterReuse * nbSamplesPerIter)
                        allSamples = allSamples(nbInitSamples+1:end, :);
                        allVals = allVals(nbInitSamples+1:end, :);
                    else
                        allSamples = allSamples(nbSamplesPerIter+1:end, :);
                        allVals = allVals(nbSamplesPerIter+1:end, :);
                    end
                end
                kls(iter) = kl;
                
                %!debug                
                if(currentDistrib.getEntropy < -5)
                    perf(iter+1:end, :) = [((iter-1:nbIter-2) * nbSamplesPerIter + nbInitSamples)' (perf(iter, 2) * ones(nbIter-iter, 1))];
                    break;
                end
            end
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        function [quadModel, etaInit, kl] = learnConcaveQuadModelWithIProj(initPolicy, maxKL, allSamples, allVals, logImportanceProbas, etaInit)
            initPrecision = initPolicy.getPrecision;
            initLogDetPrecision = initPolicy.getLogDetPrecision;
            initMu = initPolicy.getMu;
            dim = length(initMu);
            params0 = [initMu'; initPrecision(tril(ones(dim)) == 1)];
            
            optFun = @(eta)@(params) static_optimization_algs.IProjection.eval(params, initPrecision, initMu, allSamples, allVals, logImportanceProbas, eta);
            klFun = @(covOptim, detCovOptim, muOptim) static_optimization_algs.AdaptiveMoreIProj.KL(initPrecision, initLogDetPrecision, initMu, covOptim, detCovOptim, muOptim);            
            optOptim =  optimset('GradObj','on', 'DerivativeCheck', 'off', 'Display', 'off');
            for iter = 1:20
                params = fminunc(optFun(etaInit), params0, optOptim);
                muOptim = params(1:dim)';
                precisionOptim = zeros(dim);
                precisionOptim(tril(ones(dim)) == 1) = params(dim+1:end);
                precisionOptim = precisionOptim + precisionOptim' - diag(diag(precisionOptim));
                cholP = chol(precisionOptim);
                invCholP = cholP \ eye(dim);
                covOptim = invCholP * invCholP';
                logDetCovOptim = - 2 * sum(log(diag(cholP)));
                kl = klFun(covOptim, logDetCovOptim, muOptim);
                if(kl > maxKL && kl < 3 * maxKL)
                    break;
                elseif(kl < maxKL)
                    etaInit = etaInit / 2;
                else
                    etaInit = etaInit * 2;
                end
            end
%             iter
            quadModel.A = -2 * precisionOptim;
            quadModel.w = 2 * precisionOptim * muOptim';
            quadModel.b = 0;
        end
        
        function y = KL(precision, logDetP, mu, cov0, logDetCov0, mu0)
            y = .5 * (trace(precision * cov0) + (mu-mu0) * precision * (mu-mu0)' - length(mu) - logDetP - logDetCov0);
        end

    end
end
