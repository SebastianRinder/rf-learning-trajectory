classdef LocalBayes
    methods (Static)        
        %%
        function sign = getSignature(optimizerInput)
            sign = ['LocalBayes' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' optimizerInput.gpHyperOption '_' optimizerInput.samplingOption '_' optimizerInput.kernelType '_' optimizerInput.featureName...
                '_' optimizerInput.yCenteringType '_' num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.minProbaReuse) '_' num2str(optimizerInput.initVar)];
        end
        
        %%
        function [perf, kls] = optimizeStruct(optimizerInput, func)
            if(~isfield(optimizerInput, 'videoFile'))
                optimizerInput.videoFile = [];
            end
            
            [perf, kls] = static_optimization_algs.LocalBayes.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL, optimizerInput.entropyReduction, ...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbIter, optimizerInput.maxIterReuse, ...
                optimizerInput.nbPStarSamples, optimizerInput.minProbaReuse, func, optimizerInput.gpStuffPath, optimizerInput.kernelType, optimizerInput.featureFunction, optimizerInput.yCenteringType, optimizerInput.gpHyperOption, optimizerInput.samplingOption, optimizerInput.videoFile);
        end
        
        %% main function
        function [perf, allEvalsReturn, kls] = optimize(initDistrib, epsiKL, entropyReduc, nbSamplesPerIter, nbIter, maxIterReuse, nbPStarSamples, minProbaReuse, fun, gpStuffPath, kernelType, featureFunction, yCenteringType, gpHyperOption, samplingOption, videoFile)
            perf = zeros(nbIter, 2); %perf: number of samples, mean reward
            kls = zeros(nbIter, 1);
            nbSamplesForEval = 200;
            nbSamplesForAcquisition = 1200;
            nbPStarSamples = max(nbPStarSamples, nbSamplesPerIter);
            currentDistrib = initDistrib;
            allSamples = [];
            allVals = [];
            allEvalsReturn = [];
            %%% initializing GP
            static_optimization_algs.GP.addGPStuffPath(gpStuffPath);
            lik = lik_gaussian('sigma2', 0.01);
            pn = prior_logunif();
            lik = lik_gaussian(lik,'sigma2_prior', pn);
            if(strcmp(kernelType, 'sexp'))
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
                error('unrecognized kernel in LocalBayes');
            end
            gp = gp_set('lik', lik, 'cf', gpcf);
            
            gp_rec = gp;
            
            samplesArgmax = [];
            %%% animation
            if(~isempty(videoFile))
                open(videoFile);
            end
            
            %!!debug
            momentum = .6;
            targetMu = [];
            
            for iter = 1:nbIter
                %%% evaluation
%                 currPerf = static_optimization_algs.LocalBayes.distEval(currentDistrib, @(x)fun.eval(x), nbSamplesForEval);
%                 perf(iter, :) = [((iter-1) * nbSamplesPerIter), currPerf];
                
                %%% samples (acquisition)                
                if(isempty(samplesArgmax))
                    newSamples = currentDistrib.getSamples(nbSamplesPerIter);
                else            
                    %newSamples = static_optimization_algs.LocalBayes.getLocalArgmaxSamples(nbSamplesPerIter, currentDistrib, gp_rec, gp, x, y, cholPrec, nbSamplesForAcquisition);
                    %newSamples = currentDistrib.getSamples(nbSamplesPerIter);
                    if(strcmp(samplingOption, 'Acquisition'))
                        lineShuffle = randperm(size(samplesArgmax, 1));
                        newSamples = samplesArgmax(lineShuffle(1:nbSamplesPerIter), :);
                    elseif(strcmp(samplingOption, 'Random'))
                        newSamples = currentDistrib.getSamples(nbSamplesPerIter);
                    else
                        error('wrong sampling option in LocalBayes');
                    end
                end                
                newVals = fun.eval(newSamples);
                allEvalsReturn = [allEvalsReturn; newVals];
                
                allSamples = [allSamples; newSamples];
                allVals = [allVals; newVals];
                
                if(~isinf(minProbaReuse))
                    mahDistMu = pdist2(allSamples, currentDistrib.getMu, 'mahalanobis', currentDistrib.getCov) .^ 2;
                    cutOffDist = chi2inv(minProbaReuse, size(allSamples, 2));
                    localMask = mahDistMu < cutOffDist;
                    localSamples = allSamples(localMask, :);
                    localVals = allVals(localMask);
                else
                    localSamples = allSamples;
                    localVals = allVals;
                end
                
                %%% update   
                % 1- Gaussian approx of p(x = x* | Dn)
                cholPrec = currentDistrib.getCholP()';
                x = bsxfun(@minus, localSamples, currentDistrib.mu) * cholPrec; %normalize data according to current search distribution
                if(strcmp(yCenteringType, 'min'))
                    y = localVals - min(localVals);                
                elseif(strcmp(yCenteringType, 'mean'))
                    y = localVals - mean(localVals);                
                elseif(strcmp(yCenteringType, 'max'))
                    y = localVals - max(localVals);                
                else
                    error('Unrecognized y centering type in LocalBayes');
                end
                if(~isempty(featureFunction))
                    x = featureFunction(x);
                end
                [targetDistrib, samplesArgmax, gp, gp_rec] = static_optimization_algs.LocalBayes.gaussianArgmax...
                    (currentDistrib, x, y, currentDistrib.mu, cholPrec, gp, gp_rec, gpHyperOption, nbPStarSamples, nbSamplesForAcquisition, featureFunction);
                
                %!! debug
%                 if(isempty(targetMu))
%                     targetMu = targetDistrib.mu;                
%                     targetCov = targetDistrib.getCov;
%                 else
%                     targetMu = (1 - momentum) * targetMu + momentum * targetDistrib.mu;                
%                     targetCov = (1 - momentum) * targetCov + momentum * targetDistrib.getCov;
%                 end
%                 targetDistrib = static_optimization_algs.Normal(targetMu, targetCov);
                
                % 2- Quad model
                quadModel = static_optimization_algs.LocalBayes.gaussToQuad(targetDistrib);
                
                %%% current distrib plotting
                if(~isempty(videoFile))
                    frame = static_optimization_algs.LocalBayes.getFrame(fun, currentDistrib, localSamples, x, y, gp_rec, currentDistrib.mu, cholPrec, quadModel, samplesArgmax, iter, featureFunction);
                    writeVideo(videoFile, frame);
                end
                
                % 3- Updating the policy
                desiredEntrop = currentDistrib.getEntropy - entropyReduc;
                [etaStar, betaStar] = static_optimization_algs.More.optimizeDual(currentDistrib, quadModel, epsiKL, desiredEntrop);
                [currentDistrib, ~, kl] = static_optimization_algs.More.updatePolicy(currentDistrib, quadModel, [etaStar, betaStar], 1);
                
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
        
        function frame = getFrame(fun, policy, samples, x, y, gp, mu, rotMat, quadModel, samplesArgmax, iter, featureFunction)
            persistent currPlot;
            if(isempty(currPlot))
                currPlot = figure;
            end
            
            %function and distrib
            subplot(1, 3, 1);
            fun.plot();
            hold on;
            plot(samples(:, 1), samples(:, 2), '*r');
            policy.plot();            
            hold off;
            
            xLimits = xlim;
            center(1) = (xLimits(2) + xLimits(1)) / 2;
            length(1) = xLimits(2) - xLimits(1);
            yLimits = ylim;
            center(2) = (yLimits(2) + yLimits(1)) / 2;
            length(2) = yLimits(2) - yLimits(1);
          
            %gp
            subplot(1, 3, 2);
            static_optimization_algs.GP.plotGP(gp, x, y, center, length, mu, rotMat, featureFunction);

            %quad model
            subplot(1, 3, 3);
            static_optimization_algs.Quadratic.plotQuad(quadModel, center, length);
            hold on;
            plot(samplesArgmax(:, 1), samplesArgmax(:, 2), '*r');
            plot(policy.mu(1), policy.mu(2), '*g', 'MarkerSize', 20);
            hold off;
            
            text(-5, -5, ['iteration = ' num2str(iter)]);
            frame = getframe(gcf,[0 0 560 420]);
        end
        
        function eval = distEval(dist, fun, nbSamples)
            samples = dist.getSamples(nbSamples);
            vals = fun(samples);
            eval = mean(vals);
        end
        
        function [targetDistrib, samplesArgmax, gp, gp_rec] = gaussianArgmax(currentDistrib, x, y, mu, cholPrec, gp, gp_previousRec, gpHyperOption, nbPStarSamples, nbSamplesForAcquisition, featureFunction)
            %%% a- Sampling from p(x = x* | Dn)
            %hyper-param
            if(strcmp(gpHyperOption, 'MCMC'))
                try
                    burnin = 50;
                    multThin = 2;
                    gp_rec = gp_mc(gp, x, y, 'nsamples', multThin * nbPStarSamples + burnin, 'display', 0);
                    gp_rec = thin(gp_rec, burnin, multThin); %delete first 50 samples (burn-in)
                catch
                    disp('error in mcmc sampling of gp parameters');
                    gp_rec = gp_previousRec;
                end
            elseif(strcmp(gpHyperOption, 'MAP'))
                opt = optimset('TolFun', 1e-3, 'TolX', 1e-3);
                gp_rec = gp_optim(gp, x, y, 'opt', opt);
            else
                error('wrong hyperparm option in LocalBayes');
            end
            
            % fmax sampling dist
            fmaxGauss = static_optimization_algs.Normal(currentDistrib.mu, currentDistrib.getCov);
            %                 fmaxGauss = currentDistrib;
            samplesArgmax = static_optimization_algs.LocalBayes.getLocalArgmaxSamples(nbPStarSamples, fmaxGauss, gp_rec, gp, x, y, mu, cholPrec, nbSamplesForAcquisition, featureFunction);
            %%% b- ML Gaussian fit of the samples
            targetDistrib = currentDistrib.wmle(samplesArgmax);
        end
        
        function samplesArgmax = getLocalArgmaxSamples(nbPStarSamples, distrib, gp_rec, gp, x, y, mu, rotMat, nbSamplesForAcquisition, featureFunction)
            nbHyper = length(gp_rec.lik.sigma2);
            samplesArgmax = zeros(nbPStarSamples, size(rotMat, 2));
            for k = 1:nbPStarSamples
                gp.lik.sigma2 = gp_rec.lik.sigma2(mod(k - 1, nbHyper) + 1);
                try %params of sexp covariance func
                    gp.cf{1}.lengthScale = gp_rec.cf{1}.lengthScale(mod(k - 1, nbHyper) + 1);
                    gp.cf{1}.magnSigma2 = gp_rec.cf{1}.magnSigma2(mod(k - 1, nbHyper) + 1);
                catch %params of linear covariance func
                    gp.cf{1}.coeffSigma2 = gp_rec.cf{1}.coeffSigma2(mod(k - 1, nbHyper) + 1);
                end
                samplesArgmax(k, :) = static_optimization_algs.GP.localArgmaxSample(gp, x, y, distrib, nbSamplesForAcquisition, mu, rotMat, featureFunction);
            end            
        end
    end
end
