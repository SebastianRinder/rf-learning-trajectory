classdef DensityWeightedBO_core_trajectory_hyper
    methods (Static)
        %% main function
        function [newSamples, newVals, newTrajectories, gp, hyperTrace] = sample(dist, fun, lambda, xl, yl, trajectories , gp, hyperTrace, beta, gpHyperOption, ~, yCenteringType)
            % STEP 1: sample
            if(isempty(xl))
                newSamples = dist.getSamples(lambda);
                [newVals, newTrajectories] = fun.eval(newSamples);
            else
                % use current cov to standardize data
                cholPrec = dist.getCholP()';
                x = xl;
                % centering of the y values before learning the GP has an
                % impact on thompson sampling (exploration/exploitation)
                if(strcmp(yCenteringType, 'min'))
                    yCentering = min(yl);
                elseif(strcmp(yCenteringType, 'mean'))
                    yCentering = mean(yl);
                elseif(strcmp(yCenteringType, 'max'))
                    yCentering = max(yl);
                else
                    error('Unrecognized y centering type in LocalBayes');
                end
                std_yl = std(yl);
%                 std_yl = 1;
                y = (yl - yCentering) / std_yl;
                
                %hyper_param optim
                %[gp, gp_rec] = static_optimization_algs.DensityWeightedBO_core_trajectory.hyperParamOptim(x, y, gp, gp_rec, gpHyperOption, lambda);
                
                

                % use CMA-ES to maximize exp(thompson) * density acquisition
                newSamples = zeros(lambda, length(dist.mu));
                newVals = zeros(lambda, 1);
                newTrajectories = cell(lambda,1);
                
                D = fun.opts.distanceMat(x, x, trajectories, false, fun.opts);
                hyperTrace = [hyperTrace; static_optimization_algs.DensityWeightedBO_core_trajectory_hyper.optimizeHypers(hyperTrace,x,y,D, fun.opts)];                    
                %hyperTrace = [hyperTrace; 0,0];
                fun.opts.hyper = hyperTrace(end,:);
                
                for k = 1:lambda
                    if k > 1
                        D = fun.opts.distanceMat(x, x, trajectories, false, fun.opts);
                    end
                    K = fun.opts.scaleKernel(D, fun.opts.hyper);                    
                    [L, alpha] = getLowerCholesky(K, y, false);                    
                    
                    newSamples(k, :) = static_optimization_algs.DensityWeightedBO_core_trajectory_hyper.maxThompsonSamplingFmin(x,y, trajectories, L, alpha, dist, beta, fun.opts);
%                     newSamples(k, :) = static_optimization_algs.DensityWeightedBO_core_trajectory_hyper.maxExpectedImprovement(x,y, trajectories, L, alpha, dist, beta, fun.opts);
                    
                    [newVals(k), newTrajectories(k,1)] = fun.eval(newSamples(k, :));
                    x = [x; (newSamples(k, :) - dist.mu) * cholPrec];
                    y = [y; (newVals(k) - yCentering) / std_yl];
                    trajectories = [trajectories; newTrajectories(k,1)];
                end
            end
        end
        
        function newSample = maxExpectedImprovement(x,y, trajectories, L, alpha, dist, beta, opts)
            negEIFcn = @(testX) -expectedImprovement(testX, x, y, trajectories, L, alpha, opts);
            newSample = localMinSearch(negEIFcn, dist, beta);
        end
        
        function newSample = maxThompsonSampling(x,y, trajectories, L, alpha, dist, beta, opts)
            evalSamples = dist.getSamples(800);
            if (beta < 0)  % if beta < 0: disgard samples too far away from mean 
                probaThresh = -beta;
                mahDistMu = pdist2(evalSamples, dist.mu, 'mahalanobis', dist.getCov) .^ 2;
                cutOffDist = chi2inv(probaThresh, length(dist.mu));
                evalSamples = evalSamples(mahDistMu < cutOffDist, :);
                [meanVec, ~, covarianceMat] = gaussianProcess(evalSamples, x, y, trajectories, L, alpha, opts);
                cholCov = getLowerCholesky(covarianceMat, y, false);
                vals = meanVec + cholCov * randn(size(evalSamples,1), 1);
            else % if beta >= 0: weight TS value with beta * distance to mean
                keyboard;
%                         vals = static_optimization_algs.GP.gpRandTrans(gp, x, y, evalSamples, dist.mu, cholPrec, []);
%                         vals = static_optimization_algs.DensityWeightedBO_core_trajectory_hyper.trajectoryGP(x,y, trajectories, evalSamples, L, alpha, fun.opts);
%                         vals = vals + beta * dist.getLogProbas(evalSamples);
            end
            [~, argmax] = max(vals);
            newSample = evalSamples(argmax, :);
        end
        
        function newSample = maxThompsonSamplingFmin(x,y, trajectories, L, alpha, dist, beta, opts)
            negTSFcn = @(testX) -thompsonSample(testX, x, y, trajectories, L, alpha, opts);
            newSample = localMinSearch(negTSFcn, dist, beta);
        end
        
        function hypers = optimizeHypers(hyperTrace,x,y,D, opts)
            if isempty(hyperTrace)
                hyperLb(1:2) = -15;
                hyperUb(1:2) = 15;
                hyperTrace = [0,0];
            else
                hyperLb(1:2) = hyperTrace(end,1:2) - 10;
                hyperUb(1:2) = hyperTrace(end,1:2) + 10;
            end
%             hyperLb(1:2) = [-10, -40];
%             hyperUb(1:2) = [10,1];

            negHyperFcn = @(X) log(-findHypers(X, x, y, D, opts));
            hypers = globalMinSearch(negHyperFcn, hyperLb, hyperUb);
            
            hypers(1,hypers == 0) = hyperTrace(end,hypers == 0);
            
        end
    end
end
