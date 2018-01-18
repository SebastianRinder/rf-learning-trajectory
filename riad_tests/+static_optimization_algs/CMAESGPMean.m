classdef CMAESGPMean < handle
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
        function obj = CMAESGPMean(optimizerInput, func)
            obj.allEvals = [];
            obj.sampledPoints = [];
            obj.video = [];
            obj.opts = optimizerInput;
            obj.func = func;
            static_optimization_algs.GP.addGPStuffPath('GPstuff-4.7');
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
            sign = ['CMAESGPMean' '_' num2str(obj.opts.initVar) '_' num2str(obj.opts.maxEvals)];
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
            cma_opts.PopSize = '200';
%             cma_opts.ParentNumber = '100';
            xinit = 'varargin{1}{2}.opts.initDistrib.getSamples(1)';
            initSigma = obj.opts.initDistrib.getCholC;
            initSigma = initSigma(1, 1);
            cmaes_modif('static_optimization_algs.CMAESGPMean.funCallback', xinit, initSigma, cma_opts, {obj.func, obj})
            if(~isempty(obj.video))
                close(obj.video);
            end
        end
        
        function frame = getFrame(obj, policy, samples, newSamples, iter)
            if(isempty(obj.currPlot))
                obj.currPlot = figure;
            end
            
            %function and distrib
            obj.func.plot();
            hold on;
            if(~isempty(samples))
                plot(samples(:, 1), samples(:, 2), '*r');
            end
            plot(newSamples(:, 1), newSamples(:, 2), '*b');
            policy.plot();
            hold off;
            text(-5, -5, ['iteration = ' num2str(iter)]);
            frame = getframe(gcf,[0 0 560 420]);
        end
    end
    
    methods (Static)
        function [arx, arz, cmaevals, isOver] = funCallback(x, muc, sigma, cholc, lambda, varargins)
            if(nargin < 6)
                warning('this call sould not happen during the core of the optimization');
                [f, ~] = muc{:};
                val = -f.eval(x');
            else
                assert(isempty(x));
            end
            % evaluation
            [f, obj] = varargins{:};
            
            % obtain the sampled points
            beta = 1e-4;
            covcurr = (sigma .^ 2) * (cholc * cholc');
            obj.distrib = static_optimization_algs.Normal(muc, covcurr);
            [samples, evals, obj.gp, obj.gp_rec] = static_optimization_algs.DensityWeightedBO_core.sample(obj.distrib, f, obj.opts.nbSamplesPerIter, ...
                obj.xl, obj.yl, obj.gp, obj.gp_rec, beta, 'MAP', [], 'mean');
            
            obj.allEvals = [obj.allEvals; evals];
            
            obj.xl = [obj.xl; samples];
            obj.yl = [obj.yl; evals];
            
            if length(obj.yl) > obj.opts.maxIterReuse * obj.opts.nbSamplesPerIter
                obj.xl = obj.xl(obj.opts.nbSamplesPerIter+1:end, :);
                obj.yl = obj.yl(obj.opts.nbSamplesPerIter+1:end);
            end
            
            cholp = obj.distrib.getCholP()';
            x = bsxfun(@minus, obj.xl, muc') * cholp; %normalize data according to current search distribution
            y = obj.yl; % - min(obj.yl);
            [obj.gp, obj.gp_rec] = static_optimization_algs.DensityWeightedBO_core.hyperParamOptim(x, y, obj.gp, obj.gp_rec, 'MAP', obj.opts.nbSamplesPerIter);
            obj.gp = static_optimization_algs.DensityWeightedBO_core.copyHyperParam(obj.gp, obj.gp_rec, 1);
            
            
            arz = randn(length(muc), lambda);
            arx = bsxfun(@plus, sigma * (cholc * arz), muc);
            cmaevals = -(static_optimization_algs.GP.gpPredTrans(obj.gp, x, y, arx', muc', cholp, []) + beta * obj.distrib.getLogProbas(arx'))';
            
            obj.iter = obj.iter + 1;
            
            % perform special operations (e.g. plotting) after each lambda
            % points (each generation)
            if(~isempty(obj.video))
                gaussian = static_optimization_algs.Normal(muc, sigma * (cholc * cholc'));
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
