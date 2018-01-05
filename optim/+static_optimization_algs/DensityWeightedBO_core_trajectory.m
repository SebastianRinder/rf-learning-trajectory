classdef DensityWeightedBO_core_trajectory
    methods (Static)
        %% main function
        function [newSamples, newVals, newTrajectories, gp, hyperTrace] = sample(dist, func, lambda, xl, yl, trajectories , gp, hyperTrace, beta, gpHyperOption, ~, yCenteringType)
            % STEP 1: sample
            if(isempty(xl))
                newSamples = dist.getSamples(lambda);
                [newVals, newTrajectories] = func.eval(newSamples);
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
                
                func.opts.sigmaNoiseSquared = 1/sqrt(2);
                
                D = func.opts.distanceMat(x, x, trajectories, false, func.opts);
%                 func.opts.hyper = [0,0];
                func.opts.hyper = optimizeHyper(x,y,D, func.opts);    
                hyperTrace = [hyperTrace; func.opts.hyper];
                
                
                for k = 1:lambda
                    if k > 1
                        D = func.opts.distanceMat(x, x, trajectories, false, func.opts);
                    end
                    [L, alpha] = getLowerCholesky(D, y, false, func.opts.sigmaNoiseSquared);
                    bestY = max(y);
                    
%                     newSamples(k, :) = static_optimization_algs.DensityWeightedBO_core_trajectory.maxThompsonSamplingFmin(x, trajectories, L, alpha, dist, beta, func, bestY);
                    newSamples(k, :) = static_optimization_algs.DensityWeightedBO_core_trajectory.maxExpectedImprovement(x, trajectories, L, alpha, dist, beta, func, bestY);
                    
                    [newVals(k), newTrajectories(k,1)] = func.eval(newSamples(k, :));
                    x = [x; (newSamples(k, :) - dist.mu) * cholPrec];
                    y = [y; (newVals(k) - yCentering) / std_yl];
                    trajectories = [trajectories; newTrajectories(k,1)];
                end
            end
        end

        function newSample = maxExpectedImprovement(x, trajectories, L, alpha, dist, beta, func, bestY)
            negEIFcn = @(testX) -expectedImprovement(testX, x, trajectories, L, alpha, func, bestY);
            newSample = localMinSearch(negEIFcn, dist, beta);
            if func.opts.acquisitionPlot
                ub = ones(1,4);
                lb = -ub;
                selectFigure('Expected Improvement values (sorted)');
                xplot = randBound(lb,ub,10000);
                yplot = negEIFcn(xplot);
                plot(sort(-yplot));
                pause(0.1);
            end
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
        
        
    end
end
