classdef DensityWeightedBO_core_trajectory
    methods (Static)
        %% main function
        function [newSamples, newVals, newTrajectories, gp, hyperTrace] = sample(dist, func, lambda, xl, yl, trajectories , gp, hyperTrace, beta, gpHyperOption, ~, yCenteringType)
            % STEP 1: sample
            if(isempty(xl))
                newSamples = dist.getSamples(lambda);
%                 newVals = func.eval(newSamples);
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
                std_yl = max(std(yl), 1e-6);
%                 std_yl = 1;
                y = (yl - yCentering) / std_yl;
                %y(isnan(y)) = 0;

                % use CMA-ES to maximize exp(thompson) * density acquisition
                newSamples = zeros(lambda, length(dist.mu));
                newVals = zeros(lambda, 1);
                newTrajectories = cell(lambda,1);
                
                
                D = func.opts.distanceMat(x, x, trajectories, false, func.opts);
%                 if func.opts.hyperOptimize
%                     func.opts.hyper = optimizeHyper(x,y,D, func.opts);
%                 end
                hyperTrace = [hyperTrace; func.opts.hyperl.uk, func.opts.hyperl.kk, func.opts.hyperl.uu];
                
                acq = 1;
                
                for k = 1:lambda
                    if k > 1
                        D = func.opts.distanceMat(x, x, trajectories, false, func.opts);
                    end
                    K = func.opts.scaleKernel(D, [0,func.opts.hyperl.kk]);
                    [L, alpha] = getLowerCholesky(K, y, false, func.opts.noiseVariance);
                    
                    if acq == 1
                        newSamples(k, :) = static_optimization_algs.DensityWeightedBO_core_trajectory.maxThompsonSampling(x, trajectories, L, alpha, dist, beta, func);
                    elseif acq == 2
                        %newSamples(k, :) = static_optimization_algs.DensityWeightedBO_core_trajectory.maxThompsonSamplingFmin(x, y, trajectories, dist, beta, func);
                    elseif acq == 3
                        bestY = max(y);
                        newSamples(k, :) = static_optimization_algs.DensityWeightedBO_core_trajectory.maxExpectedImprovement(x, trajectories, L, alpha, dist, beta, func, bestY);
                    end
                    
                    [newVals(k), newTrajectories(k,1)] = func.eval(newSamples(k, :));
                    x = [x; (newSamples(k, :) - dist.mu) * cholPrec];
                    y = [y; (newVals(k) - yCentering) / std_yl];
                    y(isnan(y)) = 0;
                    trajectories = [trajectories; newTrajectories(k,1)];
                end
            end
        end

        function newSample = maxExpectedImprovement(x, trajectories, L, alpha, dist, beta, func, bestY)
            negEIFcn = @(testX) -expectedImprovement(testX, x, trajectories, L, alpha, func, bestY);
            newSample = localMinSearch(negEIFcn, dist, beta);
            if func.opts.acquisitionPlot
                evalSamples = dist.getSamples(10000);
                probaThresh = -beta;
                mahDistMu = pdist2(evalSamples, dist.mu, 'mahalanobis', dist.getCov) .^ 2;
                cutOffDist = chi2inv(probaThresh, length(dist.mu));
                evalSamples = evalSamples(mahDistMu < cutOffDist, :);
                selectFigure('Expected Improvement values (sorted)');
                yplot = negEIFcn(evalSamples);
                plot(sort(-yplot));
                pause(0.1);
            end
        end

        function newSample = maxThompsonSamplingFmin(x, y, trajectories, dist, beta, func)
            clear thompsonSample;   %clear persistent variables 
            negTSFcn = @(testX) -thompsonSample(testX, x, y, trajectories, func);
            newSample = localMinSearch(negTSFcn, dist, beta);            
        end
        
        function newSample = maxThompsonSampling(x, trajectories, L, alpha, dist, beta, func)
            evalSamples = dist.getSamples(300);
            
            probaThresh = -beta;
            mahDistMu = pdist2(evalSamples, dist.mu, 'mahalanobis', dist.getCov) .^ 2;
            cutOffDist = chi2inv(probaThresh, length(dist.mu));
            evalSamples = evalSamples(mahDistMu < cutOffDist, :);
            [meanVec, ~, covarianceMat] = gaussianProcess(evalSamples, x, trajectories, L, alpha, func);
            cholCov = getLowerCholesky(covarianceMat, [], false, func.opts.noiseVariance);
            vals = meanVec + cholCov * randn(size(evalSamples,1), 1);

            [~, argmax] = max(vals);
            newSample = evalSamples(argmax, :);
            if func.opts.acquisitionPlot
                selectFigure('Thompson Sample values (sorted)');
                plot(sort(vals));
                pause(0.1);
            end
        end
    end
end
