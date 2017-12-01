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
                [gp, gp_rec] = static_optimization_algs.DensityWeightedBO_core_trajectory.hyperParamOptim(x, y, gp, gp_rec, gpHyperOption, lambda);

                % use CMA-ES to maximize exp(thompson) * density acquisition
                newSamples = zeros(lambda, length(dist.mu));
                newVals = zeros(lambda, 1);
                newTrajectories = cell(lambda,1);
                
                D = trajectoryCovariance(x, x, trajectories, fun.opts);
                K = scaleKernel(D, fun.opts);
                [fun.opts.L, fun.opts.alpha] = getLowerCholesky(K, y, fun.opts, false);
                
                for k = 1:lambda
                    gp = static_optimization_algs.DensityWeightedBO_core_trajectory.copyHyperParam(gp, gp_rec, k); 
                    evalSamples = dist.getSamples(800);
                    if (beta < 0)  % if beta < 0: disgard samples too far away from mean 
                        probaThresh = -beta;
                        mahDistMu = pdist2(evalSamples, dist.mu, 'mahalanobis', dist.getCov) .^ 2;
                        cutOffDist = chi2inv(probaThresh, length(dist.mu));
                        evalSamples = evalSamples(mahDistMu < cutOffDist, :);
%                         vals = static_optimization_algs.GP.gpRandTrans(gp, x, y, evalSamples, dist.mu, cholPrec, []);                        
                        vals = static_optimization_algs.DensityWeightedBO_core_trajectory.trajectoryGP(x,y, trajectories, evalSamples, fun.opts);
                    else % if beta >= 0: weight TS value with beta * distance to mean
%                         vals = static_optimization_algs.GP.gpRandTrans(gp, x, y, evalSamples, dist.mu, cholPrec, []);
%                         vals = vals + beta * dist.getLogProbas(evalSamples);                        
                        vals = static_optimization_algs.DensityWeightedBO_core_trajectory.trajectoryGP(x,y, trajectories, evalSamples, fun.opts);
                    end
                    [~, argmax] = max(vals);
                    newSamples(k, :) = evalSamples(argmax, :);
                    [newVals(k), newTrajectories{k,1}] = fun.eval(newSamples(k, :));
                    x = [x; (newSamples(k, :) - dist.mu) * cholPrec];
                    y = [y; (newVals(k) - yCentering) / std_yl];
                    trajectories = [trajectories; newTrajectories{k,1}];
                end
            end            
        end
        
        function vals = trajectoryGP(x,y, trajectories, samples, opts)            
            [meanVec, ~, covarianceMat] = gaussianProcess(samples, x, y, trajectories, opts);
            cholCov = getLowerCholesky(covarianceMat, [], opts, false);
            vals = meanVec + cholCov * randn(size(samples,1), 1);
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
                error('wrong hyperparm option in DensityWeightedBO_trajectory');
            end
        end
    end
end
