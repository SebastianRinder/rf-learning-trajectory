classdef RepsBanditsGrad
    methods (Static)
        function sign = getSignature(optimizerInput)
            sign = ['RepsBanditsGrad' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.initVar)];
        end
        
        function [perf, kls] = optimizeStruct(optimizerInput)
            [perf, kls] = static_optimization_algs.RepsBanditsGrad.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbInitSamples, optimizerInput.nbIter, optimizerInput.maxIterReuse, optimizerInput.fun, optimizerInput.videoFile);
        end
        
        function [perf, kls] = optimize(initDistrib, epsiKL, nbSamplesPerIter, nbInitSamples, nbIter, maxIterReuse, fun, videoFile)
            perf = zeros(nbIter, 2);
            kls = zeros(nbIter, 1);
            nbSamplesForEval = 100;
            dim = initDistrib.getDim;
            lastPolicies{1} = initDistrib;
            allSamples = [];
            allVals = [];
           
            %% animation
            if(~isempty(videoFile))
                open(videoFile);
            end
            
            for iter = 1:nbIter
                %% evaluation
                numberOfSamplesEval = max(nbSamplesPerIter, nbSamplesForEval);
                if(iter == 1)
                    numberOfSamplesEval = max(numberOfSamplesEval, nbInitSamples);
                end
                newSamples = lastPolicies{end}.getSamples(numberOfSamplesEval);
                vals = fun.eval(newSamples);
                newPerf = mean(vals);
                
                %% samples
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

                %% current distrib plotting
                if(~isempty(videoFile))
                    frame = static_optimization_algs.RepsBanditsIProjection.getFrame(fun, lastPolicies{end}, allSamples);
                    writeVideo(videoFile, frame);
                end
                %% update
                if(~exist('eta', 'var'))
                    [eta, ~] = static_optimization_algs.RepsBandits.optimizeDual(allVals, epsiKL);
                end
                %%% optimizing the I-Projection
                initPrecision = lastPolicies{end}.getPrecision;
                initDetPrecision = lastPolicies{end}.getDetPrecision;
                initMu = lastPolicies{end}.getMu;
                params0 = [initMu'; initPrecision(tril(ones(dim)) == 1)];
                
                %%%Importance weights
                logImportanceProbas = zeros(length(allVals), length(lastPolicies)); 
                for k = 1:length(lastPolicies)
                    logImportanceProbas(:, k) = lastPolicies{k}.getLogProbas(allSamples);
                end
                maxLogProbas = max(logImportanceProbas, [], 2);
                logImportanceProbas = log(sum(exp(bsxfun(@minus, logImportanceProbas, maxLogProbas)), 2));
                logImportanceProbas = bsxfun(@plus, logImportanceProbas, maxLogProbas);
                
                optFun = @(eta)@(params) static_optimization_algs.IProjection.eval(params, initPrecision, initMu, allSamples, allVals, logImportanceProbas, eta);
                klFun = @(covOptim, detCovOptim, muOptim) static_optimization_algs.RepsBanditsGradReverse.KL(initPrecision, initDetPrecision, initMu, covOptim, detCovOptim, muOptim);
                [~, eta, muOptim, covOptim, kl] = static_optimization_algs.RepsBanditsGrad.loopOverEta(optFun, klFun, dim, params0, eta, epsiKL, false);
                
                lastPolicies{end+1} = static_optimization_algs.Normal(muOptim, covOptim);
                %%% delete old policies
                if(length(lastPolicies) > maxIterReuse)
                    lastPolicies = lastPolicies(2:end);
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
                if(lastPolicies{end}.getEntropy < -2)
                    perf(iter+1:end, :) = [((iter-1:nbIter-2) * nbSamplesPerIter + nbInitSamples)' (perf(iter, 2) * ones(nbIter-iter, 1))];
                    break;
                end
            end
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        function [params, eta, muOptim, covOptim, kl] = loopOverEta(fun, klFun, dim, params0, eta, epsiKL, klFromPrec)
            %%% IProjection optimization option
            optOptim =  optimset('GradObj','on', 'DerivativeCheck', 'off', 'Display', 'off');
%             optOptim.TolFun = .01;
            klLb = inf; klUp = -inf;
            etaLb = eta; etaUp = eta;
            for iterEta = 1:10
%                 params = fmincon(fun(eta), params0, [], [], [], [], [], [], [], optOptim);
                params = fminunc(fun(eta), params0, optOptim);
                muOptim = params(1:dim)';
                precisionOptim = zeros(dim);
                precisionOptim(tril(ones(dim)) == 1) = params(dim+1:end);
                precisionOptim = precisionOptim + precisionOptim' - diag(diag(precisionOptim));
                covOptim = inv(precisionOptim);
                detCovOptim = det(covOptim);
                if(klFromPrec)
                    detPrcisionOptim = 1/detCovOptim;
                    kl = klFun(precisionOptim, detPrcisionOptim, muOptim);
                else
                    kl = klFun(covOptim, detCovOptim, muOptim);
                end
                if(abs(kl - epsiKL) / epsiKL < .1)
                    break;
                elseif(kl < epsiKL)
                    klLb = kl;
                    etaUp = eta;
                else
                    klUp = kl;
                    etaLb = eta;
                end
                
                if(epsiKL < klLb)
                    etaUp = etaUp * 2;
                end
                if(epsiKL > klUp)
                    etaLb = etaLb / 2;
                end
                eta = (etaUp - etaLb) / 2 + etaLb;                
%                 params0 = [muOptim'; precisionOptim(tril(ones(dim)) == 1)];
            end
            iterEta
%             eta
%             result = [klLb kl klUp]
        end
        
    end
end
