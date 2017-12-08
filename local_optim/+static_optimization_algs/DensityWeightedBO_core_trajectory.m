classdef DensityWeightedBO_core_trajectory
    methods (Static)
        %% main function
        function [newSamples, newVals, newTrajectories, gp, gp_rec] = sample(dist, fun, lambda, xl, yl, trajectories ,gp, gp_rec, beta, gpHyperOption, ~, yCenteringType)
            % STEP 1: sample
            if(isempty(xl))
                newSamples = dist.getSamples(lambda);
                [newVals, newTrajectories] = fun.eval(newSamples);
            else
                % use current cov to standardize data

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
                
                % use CMA-ES to maximize exp(thompson) * density acquisition
                newSamples = zeros(lambda, length(dist.mu));
                newVals = zeros(lambda, 1);
                newTrajectories = cell(lambda,1);
                
                for k = 1:lambda
                    D = trajectoryCovariance(x, x, trajectories, false, fun.opts);
                    K = scaleKernel(D, fun.opts.hyper);
                    
                    [L, alpha] = getLowerCholesky(K, y, false);
                    evalSamples = dist.getSamples(800);
                    if (beta < 0)  % if beta < 0: disgard samples too far away from mean 
                        probaThresh = -beta;
                        mahDistMu = pdist2(evalSamples, dist.mu, 'mahalanobis', dist.getCov) .^ 2;
                        cutOffDist = chi2inv(probaThresh, length(dist.mu));
                        evalSamples = evalSamples(mahDistMu < cutOffDist, :);
                        
                        vals = static_optimization_algs.DensityWeightedBO_core_trajectory.trajectoryGP(x,y, trajectories, evalSamples, L, alpha, fun.opts);
                    else % if beta >= 0: weight TS value with beta * distance to mean
                        vals = static_optimization_algs.DensityWeightedBO_core_trajectory.trajectoryGP(x,y, trajectories, evalSamples, L, alpha, fun.opts);
                        vals = vals + beta * dist.getLogProbas(evalSamples);
                    end
                    [~, argmax] = max(vals);
                    newSamples(k, :) = evalSamples(argmax, :);
                    [newVals(k), newTrajectories(k,1)] = fun.eval(newSamples(k, :));
                    x = [x; (newSamples(k, :) - dist.mu) * cholPrec];
                    y = [y; (newVals(k) - yCentering) / std_yl];
                    trajectories = [trajectories; newTrajectories(k,1)];
                end
            end            
        end
        
        function vals = trajectoryGP(x,y, trajectories, samples, L, alpha, opts)            
            [meanVec, ~, covarianceMat] = gaussianProcess(samples, x, y, trajectories, L, alpha, opts);
            cholCov = getLowerCholesky(covarianceMat, y, false);
            vals = meanVec + cholCov * randn(size(samples,1), 1);
        end
    end
end
