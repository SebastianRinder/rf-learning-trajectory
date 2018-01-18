classdef RepsBanditsGradReverse
    methods (Static)
        function sign = getSignature(optimizerInput)
            sign = ['RepsBanditsGradReverse' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.initVar)];
        end
        
        function [perf] = optimizeStruct(optimizerInput)
            perf = static_optimization_algs.RepsBanditsGradReverse.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbIter, optimizerInput.fun, optimizerInput.videoFile);
        end
        
        function [perf] = optimize(initDistrib, epsiKL, nbSamplesPerIter, nbIter, fun, videoFile)
            perf = zeros(nbIter, 2);
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
                if(~exist('eta', 'var'))
                    [eta, ~] = static_optimization_algs.RepsBandits.optimizeDual(vals, epsiKL);
                end
                
                %%% optimizing the Lagrangian
                initPrecision = currentDistrib.getPrecision;
                initCov = currentDistrib.getCov;
                initMu = currentDistrib.getMu;
                initDetCov = currentDistrib.getDetCov;
                params0 = [initMu'; initPrecision(tril(ones(dim)) == 1)];
                logImportanceProbas = currentDistrib.getLogProbas(newSamples);
                
                optFun = @(eta)@(params) static_optimization_algs.RepsBanditsGradReverse.Lagrangian(params, initCov, initDetCov, initMu, newSamples, vals, logImportanceProbas, epsiKL, eta);
                klFun = @(precisionOptim, detPreOptim, muOptim) static_optimization_algs.RepsBanditsGradReverse.KL(precisionOptim, detPreOptim, muOptim, initCov, initDetCov, initMu);
                [~, eta, muOptim, covOptim] = static_optimization_algs.RepsBanditsGrad.loopOverEta(optFun, klFun, dim, params0, eta, epsiKL, true);
                
                currentDistrib = static_optimization_algs.Normal(muOptim, covOptim);
            end
            
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        function [y, dx] = Lagrangian(params, cov0, detCov0, mu0, thetas, rewards, weights, epsiKL, eta)
            % Extracting the mu and the precision matrix
            dx = zeros(length(params), 1);
            dim = length(mu0);
            mu = params(1:dim)';
            precision = zeros(dim);
            precision(tril(ones(dim)) == 1) = params(dim+1:end);
            precision = precision + precision' - diag(diag(precision));
            [cholP, p] = chol(precision);
            if(p > 0)
                y = inf;
                dx = zeros(length(params), 1);
                return;
            end
            detP = prod(diag(cholP))^2;
            invCholP = cholP \ eye(dim);
            covc = invCholP *  invCholP';
            
            % Weighted sum
            [wgSum, parDeriv] = static_optimization_algs.IProjection.WeightedGaussPrecisionNormalized(...
                mu, precision, detP, thetas, rewards, weights);
            
            currKL = static_optimization_algs.RepsBanditsGradReverse.KL(precision, detP, mu, cov0, detCov0, mu0);
            y = eta * (currKL - epsiKL) - wgSum;
            
            % Gradients
            % Mu
            dx(1:dim) = eta * precision * (mu - mu0)';% - parDeriv(1:dim);
            % Precision
            dcov = eta * (cov0 - covc + (mu-mu0)' * (mu-mu0));
            dcov = dcov - .5 * dcov .* eye(dim);
            dx(dim+1:end) = dcov(tril(ones(dim)) == 1);%- parDeriv(dim+1:end);
            dx = dx - parDeriv;
        end
        
        function y = KL(precision, detP, mu, cov0, detCov0, mu0)
            y = .5 * (trace(precision * cov0) + (mu-mu0) * precision * (mu-mu0)' - length(mu) - log(detP) - log(detCov0));
        end
        
    end
end
