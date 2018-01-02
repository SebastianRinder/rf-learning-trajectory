classdef DensityWeightedBO_core
    methods (Static)
        %% main function
        function [newSamples, newVals, gp, gp_rec] = sample(dist, fun, lambda, xl, yl, gp, gp_rec, beta, gpHyperOption, featureFunction, yCenteringType)
            % STEP 1: sample
            if(isempty(xl))
                newSamples = dist.getSamples(lambda);
                newVals = fun.eval(newSamples);
            else
                % use current cov to standardize data
                cholPrec = dist.getCholP()';
                x = bsxfun(@minus, xl, dist.mu) * cholPrec; %normalize data according to current search distribution
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
                    evalSamples = dist.getSamples(800);
                    if (beta < 0)  % if beta < 0: disgard samples too far away from mean 
                        probaThresh = -beta;
                        mahDistMu = pdist2(evalSamples, dist.mu, 'mahalanobis', dist.getCov) .^ 2;
                        cutOffDist = chi2inv(probaThresh, length(dist.mu));
                        evalSamples = evalSamples(mahDistMu < cutOffDist, :);
                        vals = static_optimization_algs.GP.gpRandTrans(gp, x, y, evalSamples, dist.mu, cholPrec, featureFunction);
                    else % if beta >= 0: weight TS value with beta * distance to mean
                        vals = static_optimization_algs.GP.gpRandTrans(gp, x, y, evalSamples, dist.mu, cholPrec, featureFunction);
                        vals = vals + beta * dist.getLogProbas(evalSamples);
                    end
                    [~, argmax] = max(vals);
                    newSamples(k, :) = evalSamples(argmax, :);
                    newVals(k) = fun.eval(newSamples(k, :));
                    x = [x; (newSamples(k, :) - dist.mu) * cholPrec];
                    y = [y; (newVals(k) - yCentering) / std_yl];
                end
            end            
        end
       
        function gp = copyHyperParam(gp, gp_rec, k)
            %             nbHyper = length(gp_rec.lik.sigma2);
            %             gp.lik.sigma2 = gp_rec.lik.sigma2(mod(k - 1, nbHyper) + 1);
            try %params of sexp covariance func
                nbHyper = length(gp_rec.cf{1}.lengthScale);
                gp.cf{1}.lengthScale = gp_rec.cf{1}.lengthScale(mod(k - 1, nbHyper) + 1);
                gp.cf{1}.magnSigma2 = gp_rec.cf{1}.magnSigma2(mod(k - 1, nbHyper) + 1);
            catch %params of linear covariance func
                nbHyper = length(gp_rec.cf{1}.coeffSigma2);
                gp.cf{1}.coeffSigma2 = gp_rec.cf{1}.coeffSigma2(mod(k - 1, nbHyper) + 1);
            end
        end
        
        function [gp, gp_rec] = hyperParamOptim(x, y, gp, gp_previousRec, gpHyperOption, nbSamplesPerIter)
            %%% a- Sampling from p(x = x* | Dn)
            % hyper-param
            if(strcmp(gpHyperOption, 'MCMC'))
                try
                    burnin = 20;
                    multThin = 8;
                    gp_rec = gp_mc(gp, x, y, 'nsamples', multThin * nbSamplesPerIter + burnin, 'display', 0);
                    gp_rec = thin(gp_rec, burnin, multThin); %delete first 50 samples (burn-in)
                catch
                    disp('error in mcmc sampling of gp parameters, taking previous params');
                    gp_rec = gp_previousRec;
                end
            elseif(strcmp(gpHyperOption, 'MAP'))
                try
                    opt = optimset('TolFun', 1e-3, 'TolX', 1e-3);
                    gp_rec = gp_optim(gp, x, y, 'opt', opt);
                catch
                    disp('error in map estimation of gp parameters, taking previous params');
                    gp_rec = gp_previousRec;
                end
            else
                error('wrong hyperparm option in DensityWeightedBO');
            end
        end
    end
end
