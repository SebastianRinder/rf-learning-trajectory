classdef More
    methods (Static)        
        %%
        function sign = getSignature(optimizerInput)
            sign = ['More' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbInitSamples) '_' ...
                num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.initVar)];
        end
        
        %%
        function [perf, kls] = optimizeStruct(optimizerInput, func)
            if(~isfield(optimizerInput, 'videoFile'))
                optimizerInput.videoFile = [];
            end
            [perf, kls] = static_optimization_algs.More.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL, optimizerInput.entropyReduction, ...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbInitSamples, optimizerInput.nbIter, optimizerInput.maxIterReuse, ...
                optimizerInput.regularization, func, optimizerInput.useImportanceSampling, optimizerInput.videoFile);
        end
        
        %% main function
        function [perf, allEvalsReturn, kls] = optimize(initDistrib, epsiKL, entropyReduc, nbSamplesPerIter, nbInitSamples, nbIter, maxIterReuse, regularization, fun, useImportanceSampling, videoFile)
            perf = zeros(nbIter, 2);
            kls = zeros(nbIter, 1);
            nbSamplesForEval = 0;
            currentDistrib = initDistrib;
            allSamples = [];
            allVals = [];
            allEvalsReturn = [];
            lastDistribs = {currentDistrib};
            
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
                
                allEvalsReturn = [allEvalsReturn; vals];

                
                %%% current distrib plotting
                if(~isempty(videoFile))
                    frame = static_optimization_algs.More.getFrame(fun, currentDistrib, allSamples, iter);
                    writeVideo(videoFile, frame);
                end
                
                %%% update   
                %1- Importance Sampling
                if(useImportanceSampling)
                    %%% Importance weights
                    logInitProbas = currentDistrib.getLogProbas(allSamples);
%                     logImportanceProbas = static_optimization_algs.NonParamImportanceSampling.getLogProbas(allSamples, currentDistrib.getPrecision, logInitProbas);
                    logImportanceProbas = zeros(size(allSamples, 1), length(lastDistribs));
                    for k = 1:length(lastDistribs)                        
                        logImportanceProbas(:, k) = lastDistribs{k}.getLogProbas(allSamples); %biased if nbInitSamples != nbSamplesPerIter
                    end
                    maxLogip = max(logImportanceProbas, [], 2);
                    logImportanceProbas = maxLogip + log(sum(exp(bsxfun(@minus, logImportanceProbas, maxLogip)), 2));
                    
                    %%% Learn Concave Quadratic Reward Model
                    sampleWeights = exp(logInitProbas-logImportanceProbas);
                else
                    sampleWeights = ones(size(allVals));
                end
                
                %2- Learning quadratic model
%                 [~, quadModel] = static_optimization_algs.Quadratic.learnConcaveQuadModel(allSamples, allVals, regularization, sampleWeights, .999, 1e5);                 
                [~, quadModel] = static_optimization_algs.Quadratic.learnQuadModel(allSamples, allVals, regularization, sampleWeights); 
                
                %3- Updating the policy
                desiredEntrop = currentDistrib.getEntropy - entropyReduc;
                [etaStar, betaStar] = static_optimization_algs.More.optimizeDual(currentDistrib, quadModel, epsiKL, desiredEntrop);
                [currentDistrib, ~, kl] = static_optimization_algs.More.updatePolicy(currentDistrib, quadModel, [etaStar, betaStar], 1);
                lastDistribs{end + 1} = currentDistrib;
                
                %%% delete old samples
                if(iter > maxIterReuse)
                    lastDistribs = lastDistribs(2:end);
                    if(iter == maxIterReuse + 1)
                        allSamples = allSamples(nbInitSamples+1:end, :);
                        allVals = allVals(nbInitSamples+1:end, :);
                    else
                        allSamples = allSamples(nbSamplesPerIter+1:end, :);
                        allVals = allVals(nbSamplesPerIter+1:end, :);
                    end
                end
                kls(iter) = kl;
                
%                  %!debug                
%                  if(currentDistrib.getEntropy < -2)
%                      perf(iter+1:end, :) = [((iter-1:nbIter-2) * nbSamplesPerIter + nbInitSamples)' (perf(iter, 2) * ones(nbIter-iter, 1))];
%                      break;
%                  end
            end
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        %% Returns the dual param eta corresponding to max_p f(a) s.t. KL(p | initPolicy) < epsiAction
        % initPolicy is of Normal class.
        % quadFunction is a struct containing A, w, b s.t. the reward f(a) = 1/2 a' A a + a' w + b.     
        % epsiAction upper bound of the KL
        %        
        function [etaStar, betaStar] = optimizeDual(initPolicy, quadFunction, epsiAction, desiredEntrop, initEta, initBeta)
            % initial configuration
            initPolicyParams.mu = initPolicy.getMu';
            initPolicyParams.P = initPolicy.getPrecision;
            initPolicyParams.cholP = initPolicy.getCholP;
            initPolicyParams.logDetP = 2 * sum(log(diag(initPolicyParams.cholP)));
            if(~exist('desiredEntrop', 'var'))
                desiredEntrop = [];
            end
            fun = @(param) static_optimization_algs.More.dualFunction(param, initPolicyParams, quadFunction, epsiAction, desiredEntrop);
            if(~exist('initEta', 'var') || isempty(initEta))
                initEta = 1;
            end
            if(~exist('initBeta', 'var') || isempty(initBeta))
                initBeta = 0;
            end
            initParams =  initEta;
            paramsLb = 1e-16;
            if(~isempty(desiredEntrop))
                initParams(2) = initBeta;
                paramsLb(2) = -inf;
            end
            
            % optimizing the dual
            options = optimset('GradObj', 'on', 'Display', 'off', 'Algorithm', 'trust-region-reflective', 'TolX', 1e-16, 'TolFun', 1e-16, ...
                'DerivativeCheck', 'off');
            redoOptim = true;
            while(redoOptim == true)
                try
                    paramStar = fmincon(fun, initParams, [] ,[], [], [], paramsLb, [], [], options);
                    redoOptim = false;
                catch
                    initParams(1) = initParams(1) * 2;
                end
            end
            
            etaStar = paramStar(1);
            if(length(paramStar) > 1)
                betaStar = paramStar(2);
            end
        end
        
        %% Dual function and derivative
        % initPolicyParams is a struct containing: P
        % (precision), logDetP, cholP, and mu (mean) 
        % quadFunction is a struct containing A, w, b s.t. the reward f(a) = 1/2 a' A a + a' w + b.   
        % epsiAction upper bound of the KL
        function [g, dg] = dualFunction(params, initPolicyParams, quadFunction, epsiAction, desiredEntropy)
            eta = params(1);
            if(length(params) > 1)
                beta = params(2);
            else
                beta = 0;
                desiredEntropy = 0;
            end
            
            % Computing precision and mean of new policy
            F = eta * initPolicyParams.P - quadFunction.A;
            fc = eta * initPolicyParams.P  * initPolicyParams.mu + quadFunction.w;
            [cholF, notPsd] = chol(F);
            if(notPsd > 0)
                disp('--> F not positive definite');
                g = inf;
                dg = 0;
                return;
            end      
            
            % Computing the dual function g
            detTerm = 2 * sum(log(diag(cholF)));
            dim = length(quadFunction.w);
            logDetF = dim * log((eta + beta) * 2 * pi) - detTerm; % log |2 pi eta F^-1|
            logDetQ = dim * log(2 * pi) - initPolicyParams.logDetP; % log |2 pi C|
            
            Ffc = fc' / cholF;
            dualFTerm = Ffc * Ffc';
            Pw = initPolicyParams.mu' * initPolicyParams.cholP';
            dualPTerm = Pw * Pw';
            g = eta * epsiAction - beta * desiredEntropy + .5 * (- eta * dualPTerm + dualFTerm  ...
                - eta * logDetQ + (eta + beta) * logDetF);
            
            %dg
            if(nargout > 1)
                cFFfc = (cholF \ Ffc');
                temp = cFFfc' * initPolicyParams.cholP';
                cQF = initPolicyParams.cholP / cholF;
                dg = epsiAction - .5 * dualPTerm  - .5 * (temp * temp') + temp * Pw' - .5 * logDetQ ...
                    + .5 * logDetF - .5 * (eta + beta) * trace(cQF * cQF') + .5 * dim;
                if(length(params) > 1)
                    dg(2) = .5 * (dim + logDetF) - desiredEntropy;
                end
            end
        end
        
        
        %% compute the new policy 
        % initPolicy is of Normal class.
        % quadFunction is a struct containing A, w, b s.t. the reward f(a) = 1/2 a' A a + a' w + b.     
        function [newPol, rewardImprovement, kl, entropy]  = updatePolicy(initPolicy, quadFunction, params, verbose, targetEntrop)
            % policy update
            eta = params(1);
            P = initPolicy.getPrecision;
            F = P - quadFunction.A / eta;
            mu = initPolicy.getMu;
            fc = P * mu' + quadFunction.w / eta;
            cholF = chol(F);            
            Ffc = (fc' / cholF)'; 
            newMean = cholF \ Ffc;
            if(exist('targetEntrop', 'var') && ~isempty(targetEntrop))
                logDetF = 2 * sum(log(diag(cholF)));
                dim = length(mu);
                F = (2 * pi * exp(1 - 2 / dim * targetEntrop - logDetF / dim)) * F ;
            elseif(length(params) > 1)
                beta = params(2);
                F = (eta / (eta + beta)) * F;
            end
            newPol = static_optimization_algs.Normal(newMean, F, 2); %2 indicates that F is the precision matrix
            
            % computing additional variables
            if(nargout > 1 || (exist('verbose', 'var') && verbose == 1))
                initReward = mu * quadFunction.w + .5 * mu * quadFunction.A * mu' + .5 * trace(initPolicy.cholC * quadFunction.A * initPolicy.cholC');
                newReward = newPol.getMu * quadFunction.w + .5 * newPol.getMu * quadFunction.A * newPol.getMu' + .5 * trace(newPol.cholC * quadFunction.A * newPol.cholC');
                rewardImprovement = (newReward - initReward) / abs(initReward);
                entropy = newPol.getEntropy;
                kl = newPol.getKL(initPolicy);
                if(exist('verbose', 'var') && verbose == 1)
                    fprintf('relative reward change %g, kl %g, entropy %g, entropy change %g\n', newReward, rewardImprovement, kl, entropy, entropy - initPolicy.getEntropy);
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
