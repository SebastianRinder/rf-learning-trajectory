classdef CMAESThompson < handle
    properties
        iter;
        allEvals;
        sampledPoints;
        opts;
        func;
        currPlot;
        gp;
        gp_rec;
        video;
        xl;
        yl;
        distrib;
    end
    
    methods
        %%
        function obj = CMAESThompson(optimizerInput, func)
            obj.allEvals = [];
            obj.sampledPoints = [];
            obj.video = [];
            obj.opts = optimizerInput;
            obj.func = func;
            static_optimization_algs.GP.addGPStuffPath('GPstuff-4.7/');
            lik = lik_gaussian('sigma2', 1e-6, 'sigma2_prior', prior_fixed());
            gpcf = gpcf_matern52('lengthScale', 1, 'magnSigma2', 0.01);
            % hyper-param priors of GP
            pl = prior_unif();
            pm = prior_sqrtunif();
            gpcf = gpcf_matern52(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
            obj.gp = gp_set('lik', lik, 'cf', gpcf);
            obj.gp_rec = obj.gp; % struct storing list of hyper_parames. Needed for possible compatibility with MCMC
            obj.distrib = obj.opts.initDistrib;
            obj.xl = [];
            obj.yl = [];
            obj.iter = 0;
        end
        
        function sign = getSignature(obj)
            sign = ['CMAESThompson' '_' num2str(obj.opts.initVar) '_' num2str(obj.opts.maxEvals)];
        end
        
        %% main function
        function optimize(obj)
            if(~isempty(obj.video))
                open(obj.video);
            end
            cma_opts = cmaes_modif('defaults');
            cma_opts.EvalInitialX = 'no';
            cma_opts.Restarts = '0';
            cma_opts.CMA.active = '0';
            cma_opts.DiagonalOnly = '0';
            cma_opts.ExternalSampling = '1';
%             cma_opts.TolFun = 1e-3;
            xinit = obj.opts.xinit;
            initSigma = obj.opts.initDistrib.getCholC;
            initSigma = initSigma(1, 1);
            cmaes_modif('static_optimization_algs.CMAESThompson.funCallback', xinit, initSigma, cma_opts, {obj.func, obj})
            if(~isempty(obj.video))
                close(obj.video);
            end
        end
    end
    
    methods (Static)
        function [arx, arz, cmaevals, isOver] = funCallback(x, muc, sigma, cholc, lambda, varargins)
            % extract function and object
            if(nargin < 6)
                warning('this call sould not happen during the core of the optimization');
                [f, ~] = muc{:};
                arx = -f.eval(x');
                return;
            else
                assert(isempty(x));
            end
            [f, obj] = varargins{:};
            
            % thompson sampling
            beta = -.8;
            covcurr = (sigma .^ 2) * (cholc * cholc');
            obj.distrib = static_optimization_algs.Normal(muc, covcurr);
            [samples, evals, obj.gp, obj.gp_rec] = static_optimization_algs.DensityWeightedBO_core.sample(obj.distrib, f, lambda, ...
                obj.xl, obj.yl, obj.gp, obj.gp_rec, beta, 'MAP', [], 'mean');
            
            obj.allEvals = [obj.allEvals; evals];            
            obj.xl = [obj.xl; samples];
            obj.yl = [obj.yl; evals];            
            if length(obj.yl) > obj.opts.maxIterReuse * lambda
                obj.xl = obj.xl(lambda+1:end, :);
                obj.yl = obj.yl(lambda+1:end);
            end

            arx = samples';
            arz = zeros(length(muc), lambda);
            for k = 1:lambda
                arz(:, k) = (cholc \ (arx(:, k) - muc)) ./ sigma;
            end
            cmaevals = -evals';

            %plot if needed
            obj.iter = obj.iter + 1;                        
            if(~isempty(obj.video))
                cholp = obj.distrib.getCholP()';
                x = bsxfun(@minus, obj.xl, muc') * cholp; %normalize data according to current search distribution
                y = obj.yl; % - min(obj.yl);
                gaussian = static_optimization_algs.Normal(muc, covcurr);
                frame = static_optimization_algs.DensityWeightedBO.getFrame(f, gaussian, obj.xl, x, y, obj.gp, muc', cholp, samples, obj.iter, []);
                writeVideo(obj.video, frame);
            end
            isOver = false;
            if obj.iter > obj.opts.maxIter
                isOver = true;
            end
        end
    end
end
