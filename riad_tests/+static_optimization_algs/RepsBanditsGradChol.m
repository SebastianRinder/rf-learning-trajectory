classdef RepsBanditsGradChol
    methods (Static)
        function sign = getSignature(optimizerInput)
            sign = ['RepsBanditsGradChol' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.initVar)];
        end
        
        function [perf, kls] = optimizeStruct(optimizerInput)
            [perf, kls] = static_optimization_algs.RepsBanditsGradChol.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbIter, optimizerInput.fun, optimizerInput.videoFile);
        end
        
        function [perf, kls] = optimize(initDistrib, epsiKL, nbSamplesPerIter, nbIter, fun, videoFile)
            perf = zeros(nbIter, 2);
            kls = zeros(nbIter, 1);
            nbSamplesForEval = 100;
            currentDistrib = initDistrib;
            dim = currentDistrib.getDim;
           
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
%                 if(~exist('eta', 'var'))
                    [eta, ~] = static_optimization_algs.RepsBandits.optimizeDual(vals, epsiKL);
%                 end
                
                %%% optimizing the I-Projection
                initPrecision = currentDistrib.getPrecision;
                initDetPrecision = currentDistrib.getDetPrecision;
                initMu = currentDistrib.getMu;
                initChol = currentDistrib.getCholP;
                params0 = [initMu'; initChol(triu(ones(dim)) == 1)];
                logImportanceProbas = currentDistrib.getLogProbas(newSamples);
                
                optFun = @(eta)@(params) static_optimization_algs.IProjection.eval_chol(params, initPrecision, initMu, newSamples, vals, logImportanceProbas, eta);
                klFun = @(covOptim, detCovOptim, muOptim) static_optimization_algs.RepsBanditsGradReverse.KL(initPrecision, initDetPrecision, initMu, covOptim, detCovOptim, muOptim);
                [~, eta, muOptim, cholOptim] = static_optimization_algs.RepsBanditsGradChol.loopOverEta(optFun, klFun, dim, params0, eta, epsiKL, false);
                
                currentDistrib = static_optimization_algs.Normal(muOptim, cholOptim, true);
                kls(iter) = klFun(currentDistrib.getCov, currentDistrib.getDetCov, muOptim);
            end
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        function [params, eta, muOptim, cholOptim] = loopOverEta(fun, klFun, dim, params0, eta, epsiKL, klFromPrec)
            %%% IProjection optimization option
            optOptim =  optimset('GradObj','on', 'DerivativeCheck', 'on', 'Display', 'off');
%             optOptim.TolFun = .01;
            klLb = inf; klUp = -inf;
            etaLb = eta; etaUp = eta;
            lbParams = -inf * ones(size(params0));
            diagCholIndices = zeros(dim);
            diagCholIndices(triu(ones(dim)) == 1) = 1:(dim * (dim+1) / 2);
            lbParams(dim + diag(diagCholIndices)) = 0;
            for iterEta = 1:10
                params = fmincon(fun(eta), params0, [], [], [], [], lbParams, [], [], optOptim);
                muOptim = params(1:dim)';
                cholOptim = zeros(dim);
                cholOptim(triu(ones(dim)) == 1) = params(dim+1:end);
                precisionOptim = cholOptim' * cholOptim;
                invCholP = cholOptim \ eye(dim);
                covOptim = invCholP *  invCholP';
                detCovOptim = prod(diag(invCholP))^2;
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
%                 params0 = [muOptim'; cholOptim(triu(ones(dim)) == 1)];
            end
%             iterEta
%             eta
%             result = [klLb kl klUp]
        end
        
    end
end
