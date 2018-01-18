classdef DensityWeightedBO_core_optim
    methods (Static)
        %% main function
        function [newSamples, newVals, gp, gp_rec] = sample(dist, fun, lambda, xl, yl, gp, gp_rec, beta, gpHyperOption, featureFunction, yCenteringType)
            % STEP 1: sample
            if(isempty(xl))
                newSamples = dist.getSamples(lambda);
                newVals = fun.eval(newSamples);
            else
                % optimizer (of Thompson sampling) options
                cma_opts = cmaes('defaults');
                cma_opts.DispFinal = 'off';
                cma_opts.DispModulo = 'Inf';
                cma_opts.Restarts = '0';
                cma_opts.CMA.active = '0';
                cma_opts.DiagonalOnly = '0';
                cma_opts.MaxFunEvals = 2000;
                cma_opts.TolFun = 1e-3;
                cma_xstart = zeros(size(dist.mu))';
                cma_insigma = exp(sum(log(diag(dist.getCholC))) / length(dist.mu));
                
                % use current cov to standardize data
                cholPrec = dist.getCholP()';
                x = bsxfun(@minus, xl, dist.mu) * cholPrec; %normalize data according to current search distribution
                
                % centering of the y values before learning the GP has an
                % impact on Thompson sampling (exploration/exploitation)
                if(strcmp(yCenteringType, 'min'))
                    yCentering = min(yl);
                elseif(strcmp(yCenteringType, 'mean'))
                    yCentering = mean(yl);
                elseif(strcmp(yCenteringType, 'max'))
                    yCentering = max(yl);
                else
                    error('Unrecognized y centering type in LocalBayes');
                end
                std_yl = std(yl) + 1e-6;
                %                 std_yl = 1;
                y = (yl - yCentering) / std_yl;
                if(~isempty(featureFunction))
                    x = featureFunction(x);
                end
                
                %hyper_param optim
                [gp, gp_rec] = static_optimization_algs.DensityWeightedBO_core.hyperParamOptim(x, y, gp, gp_rec, gpHyperOption, lambda);
                
                % use CMA-ES to maximize exp(thompson) * density acquisition
                newSamples = zeros(lambda, length(dist.mu));
                newVals = zeros(lambda, 1);
                for k = 1:lambda
                    gp = static_optimization_algs.DensityWeightedBO_core.copyHyperParam(gp, gp_rec, k);
                    add_x = dist.getSamples(100);
                    add_y = static_optimization_algs.GP.gpRandTrans(gp, x, y, add_x, dist.mu, cholPrec, featureFunction);
                    x_ts = [x; bsxfun(@minus, add_x, dist.mu) * cholPrec];
                    y_ts = [y; add_y];
                    ts_obs = static_optimization_algs.BoundedThompsonSampling(gp, x_ts, y_ts, -beta, cma_insigma);
                    [x_new, f_new] = cmaes(@(xt) -ts_obs.eval(xt'), cma_xstart, cma_insigma, cma_opts);
                    % check if value from cma is better than random samples (if hyper-param optim goes wrong, acquisition function can degenerate to zero function)
                    mahDistMu = pdist2(add_x, dist.mu, 'mahalanobis', dist.getCov) .^ 2;
                    cutOffDist = chi2inv(-beta, length(dist.mu));
                    add_x = add_x(mahDistMu < cutOffDist, :);
                    add_y = add_y(mahDistMu < cutOffDist);
                    [v, id] = max(add_y);
                    if(-f_new < v)
                        newSamples(k, :) = add_x(id, :);
                        x = [x; (newSamples(k, :) - dist.mu) * cholPrec];
                    else
                        newSamples(k, :) = (x_new' / cholPrec + dist.mu)';
                        x = [x; x_new'];
                    end
                    newVals(k) = fun.eval(newSamples(k, :));
                    y = [y; (newVals(k) - yCentering) / std_yl];
                end
            end
        end
    end
end
