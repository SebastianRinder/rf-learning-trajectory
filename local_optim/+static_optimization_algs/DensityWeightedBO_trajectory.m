classdef DensityWeightedBO_trajectory
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['DensityWeightedBO_trajectory' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' optimizerInput.gpHyperOption '_' optimizerInput.samplingOption '_' optimizerInput.kernelType '_' optimizerInput.featureName...
                '_' optimizerInput.yCenteringType '_' num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.minProbaReuse) '_' num2str(optimizerInput.initVar)];
        end
        
        %%
        function [perf, kls] = optimizeStruct(optimizerInput, func)
            if(~isfield(optimizerInput, 'videoFile'))
                optimizerInput.videoFile = [];
            end
            
            [perf, kls] = static_optimization_algs.DensityWeightedBO_trajectory.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL, optimizerInput.entropyReduction, ...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbIter, optimizerInput.maxIterReuse, ...
                optimizerInput.nbPStarSamples, optimizerInput.minProbaReuse, func, optimizerInput.gpStuffPath, optimizerInput.kernelType, optimizerInput.featureFunction, optimizerInput.yCenteringType, optimizerInput.gpHyperOption, optimizerInput.samplingOption, optimizerInput.videoFile);
        end
        
        %% main function
        function [allVals, kls] = optimize(initDistrib, epsiKL, entropyReduc, nbSamplesPerIter, nbIter, maxIterReuse, nbPStarSamples, minProbaReuse, fun, gpStuffPath, kernelType, featureFunction, yCenteringType, gpHyperOption, samplingOption, videoFile)
            kls = zeros(nbIter, 1);
            nbPStarSamples = max(nbPStarSamples, nbSamplesPerIter);
            currentDistrib = initDistrib;
            storedSamples = [];
            storedVals = [];
            storedTrajectories = [];
            allVals = [];
            %%% initializing GP
            static_optimization_algs.GP.addGPStuffPath(gpStuffPath);
            lik = lik_gaussian('sigma2', 1e-6, 'sigma2_prior', prior_fixed());
            %             pn = prior_logunif();
            %             lik = lik_gaussian(lik,'sigma2_prior', pn);
            if(strcmp(kernelType, 'matern52'))
                gpcf = gpcf_matern52('lengthScale', 1, 'magnSigma2', 0.01);
                % hyper-param priors of GP
                pl = prior_unif();
                pm = prior_sqrtunif();
                gpcf = gpcf_matern52(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
            elseif(strcmp(kernelType, 'sexp'))
                gpcf = gpcf_sexp('lengthScale', 1, 'magnSigma2', 0.01);
                % hyper-param priors of GP
                pl = prior_unif();
                pm = prior_sqrtunif();
                gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
            elseif(strcmp(kernelType, 'linear'))
                gpcf = gpcf_linear('coeffSigma2', 1);
                % hyper-param priors of GP
                pl = prior_logunif();
                gpcf = gpcf_linear(gpcf, 'coeffSigma2_prior', pl);
            else
                error('unrecognized kernel in DensityWeightedBO_trajectory');
            end
            gp = gp_set('lik', lik, 'cf', gpcf);
            gp_rec = gp;
            
            %%% animation
            if(~isempty(videoFile))
                open(videoFile);
            end
            
            localSamples = [];
            localVals = [];
            localTrajectories = [];
            for iter = 1:nbIter
                neg_proba_thresh = -.8;
                [newSamples, newVals, newtrajectory, gp, gp_rec] = static_optimization_algs.DensityWeightedBO_core_trajectory.sample(currentDistrib, fun, nbSamplesPerIter, ...
                    localSamples, localVals, localTrajectories, gp, gp_rec, neg_proba_thresh, gpHyperOption, featureFunction, yCenteringType);
                
                % STEP 2: evaluate and manage dataset
                allVals = [allVals; newVals];
                storedSamples = [storedSamples; newSamples];
                storedVals = [storedVals; newVals];
                storedTrajectories = [storedTrajectories; newtrajectory];
                % delete old samples
                if(iter > maxIterReuse)
                    storedSamples = storedSamples(nbSamplesPerIter+1:end, :);
                    storedVals = storedVals(nbSamplesPerIter+1:end, :);
                    storedTrajectories = storedTrajectories(nbSamplesPerIter+1:end, :);
                end
                % select local samples according to search distrib
                if(~isinf(minProbaReuse))
                    mahDistMu = pdist2(allSamples, currentDistrib.getMu, 'mahalanobis', currentDistrib.getCov) .^ 2;
                    cutOffDist = chi2inv(minProbaReuse, size(allSamples, 2));
                    localMask = mahDistMu < cutOffDist;
                    localSamples = storedSamples(localMask, :);
                    localVals = storedVals(localMask);
                    localTrajectories = storedTrajectories(localMask);
                else
                    localSamples = storedSamples;
                    localVals = storedVals;
                    localTrajectories = storedTrajectories;
                end
                
                % STEP 3: update search distribution
                % GP hyper-param optim
                cholPrec = currentDistrib.getCholP()';
                x = bsxfun(@minus, localSamples, currentDistrib.mu) * cholPrec; %normalize data according to current search distribution
                yCentering = 0;
                if(strcmp(yCenteringType, 'min'))
                    yCentering = min(localVals);
                elseif(strcmp(yCenteringType, 'mean'))
                    yCentering = mean(localVals);
                elseif(strcmp(yCenteringType, 'max'))
                    yCentering = max(localVals);
                else
                    error('Unrecognized y centering type in LocalBayes');
                end
                y = localVals - yCentering;
                if(~isempty(featureFunction))
                    x = featureFunction(x);
                end
                
                %hyper param optim
%                 [gp, gp_rec] = static_optimization_algs.DensityWeightedBO_core.hyperParamOptim(x, y, gp, gp_rec, gpHyperOption, nbSamplesPerIter);
%                 gp = static_optimization_algs.DensityWeightedBO_core.copyHyperParam(gp, gp_rec, 1);
                
                % learn quad model from GP
%                 quadModelSamples = currentDistrib.getSamples(400);
%                 modeVals = static_optimization_algs.GP.gpPredTrans(gp, x, y, quadModelSamples, currentDistrib.mu, cholPrec, featureFunction);
%                 [~, quadModel] = static_optimization_algs.Quadratic.learnQuadModel(quadModelSamples, modeVals, 1e-6, []);

                % cma-es like: GP mean
%                 nbSamplesForTarget = 200;
%                 targetDistribSamples = currentDistrib.getSamples(2 * nbSamplesForTarget);
%                 modelVals = static_optimization_algs.GP.gpPredTrans(gp, x, y, targetDistribSamples, currentDistrib.mu, cholPrec, featureFunction);
%                 weights = log(nbSamplesForTarget + 1/2)-log(1:nbSamplesForTarget)'; % muXone array for weighted recombination
%                 [~, sorted_idx] = sort(modelVals, 1, 'descend');
%                 targetDistribSamples = targetDistribSamples(sorted_idx(1:nbSamplesForTarget), :);
%                 targetDistrib = currentDistrib.wmle(targetDistribSamples, weights);
%                 quadModel = static_optimization_algs.DensityWeightedBO.gaussToQuad(targetDistrib);
                
                % cma-es like: Sampled values
                nbSamplesForTarget = length(localVals) / 2;
                weights = log(nbSamplesForTarget + 1/2) - log(1:nbSamplesForTarget)'; % muXone array for weighted recombination
                [~, sorted_idx] = sort(localVals, 1, 'descend');
                targetDistribSamples = localSamples(sorted_idx(1:nbSamplesForTarget), :);
                targetDistrib = currentDistrib.wmle(targetDistribSamples, weights);
                quadModel = static_optimization_algs.DensityWeightedBO_trajectory.gaussToQuad(targetDistrib);

                % cma-es like: p_star samples
%                 nbSamplesForTarget = nbSamplesPerIter / 2;
%                 weights = log(nbSamplesForTarget + 1/2) - log(1:nbSamplesForTarget)'; % muXone array for weighted recombination
%                 [~, sorted_idx] = sort(newVals, 1, 'descend');
%                 targetDistribSamples = newSamples(sorted_idx(1:nbSamplesForTarget), :);
%                 targetDistrib = currentDistrib.wmle(targetDistribSamples, weights);
%                 quadModel = static_optimization_algs.DensityWeightedBO.gaussToQuad(targetDistrib);



                %                 % learn quad model from Data
                %                 logInitProbas = currentDistrib.getLogProbas(localSamples);
                %                 logImportanceProbas = static_optimization_algs.NonParamImportanceSampling.getLogProbas(localSamples, currentDistrib.getPrecision, logInitProbas);
%                                 sampleWeights = exp(logInitProbas-logImportanceProbas);
%                                 [~, quadModel] = static_optimization_algs.Quadratic.learnQuadModel(localSamples, localVals, 1e-6, sampleWeights);
                
                
%                 logInitProbas = currentDistrib.getLogProbas(localSamples);
%                 [~, quadModel] = static_optimization_algs.Quadratic.learnQuadModel(localSamples, localVals, 1e-6, exp(logInitProbas));
                
                % plotting current iteration if video recoreded
                if(~isempty(videoFile))
                    frame = static_optimization_algs.DensityWeightedBO_trajectory.getFrame(fun, currentDistrib, localSamples, x, y, gp, currentDistrib.mu, cholPrec, newSamples, iter, featureFunction);
                    writeVideo(videoFile, frame);
                end
                
                % solve the optim problem
                desiredEntrop = currentDistrib.getEntropy - entropyReduc;
                [etaStar, betaStar] = static_optimization_algs.More.optimizeDual(currentDistrib, quadModel, epsiKL, desiredEntrop);
                [currentDistrib, ~, kl] = static_optimization_algs.More.updatePolicy(currentDistrib, quadModel, [etaStar, betaStar], 1);
                kls(iter) = kl;
            end
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        function quadModel = gaussToQuad(normal)
            quadModel.b = 0;
            quadModel.A = -normal.getPrecision();
            quadModel.w = normal.getPrecision() * normal.mu()';
        end

        
        function frame = getFrame(fun, policy, samples, x, y, gp, mu, rotMat, newSamples, iter, featureFunction)
            persistent currPlot;
            if(isempty(currPlot))
                currPlot = figure;
            end
            
            %function and distrib
            subplot(1, 2, 1);
            fun.plot();
            hold on;
            plot(samples(:, 1), samples(:, 2), '*r');
            plot(newSamples(:, 1), newSamples(:, 2), '*b');
            policy.plot();
            hold off;
            
            %gp
            xLimits = xlim;
            center(1) = (xLimits(2) + xLimits(1)) / 2;
            length(1) = xLimits(2) - xLimits(1);
            yLimits = ylim;
            center(2) = (yLimits(2) + yLimits(1)) / 2;
            length(2) = yLimits(2) - yLimits(1);
            
            subplot(1, 2, 2);
            xlim(xLimits);
            ylim(yLimits);
            static_optimization_algs.GP.plotGP(gp, x, y, center, length, mu, rotMat, featureFunction);
            
            text(-5, -5, ['iteration = ' num2str(iter)]);
            frame = getframe(gcf,[0 0 560 420]);
        end
    end
end
