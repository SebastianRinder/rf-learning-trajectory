classdef CMAESCallbackWrapper < handle
    properties
        allEvals;
        sampledPoints;
        optimizerInput;
        func;
        currPlot;
        video;
    end
    
    methods
        %%
        function obj = CMAESCallbackWrapper(optimizerInput, func)
            obj.allEvals = [];
            obj.sampledPoints = [];
            obj.video = [];
            obj.optimizerInput = optimizerInput;
            obj.func = func;
        end
        
        function sign = getSignature(obj)
            sign = ['CMAESCallbackWrapper' '_' num2str(obj.optimizerInput.initVar) '_' num2str(obj.optimizerInput.maxEvals)];
        end
        
        %% main function
        function optimize(obj)
            if(~isempty(obj.video))
                open(obj.video);
            end
            opts = cmaes_modif('defaults');
            opts.EvalInitialX = 'no';
            opts.Restarts = '0';
            opts.CMA.active = '0';
            opts.DiagonalOnly = '0';
            opts.MaxFunEvals = obj.optimizerInput.maxEvals;
%             xinit = ['rand(' num2str(length(obj.optimizerInput.initDistrib.mu)) ', 1) * 20 - 10'];
            opts.TolFun = 1e-3;
%             opts.ParentNumber = 'floor((4 + floor(3*log(N))) / 2)';
%             xinit = 'varargin{1}{2}.optimizerInput.initDistrib.getSamples(1)';
            xinit = obj.optimizerInput.initDistrib.mu;
            initSigma = obj.optimizerInput.initDistrib.getCholC;
            initSigma = initSigma(1,1);
            obj.optimizerInput.initDistrib = static_optimization_algs.Normal(obj.optimizerInput.initDistrib.mu, 4 * obj.optimizerInput.initDistrib.getCovariance);
            [xmin, fxmin, a, b, c, d] = cmaes_modif('static_optimization_algs.CMAESCallbackWrapper.funCallback', xinit, initSigma, opts, {obj.func, obj})
            if(~isempty(obj.video))
                close(obj.video);
            end
        end
        
        function endOfGenerationCallBack(obj, muc, cholc, lambda)
            if(~isempty(obj.video))
                gaussian = static_optimization_algs.Normal(muc, cholc * cholc');
                frame = obj.getFrame(gaussian, [], ...
                    obj.sampledPoints', length(obj.allEvals) / lambda);
                writeVideo(obj.video, frame);
                obj.sampledPoints = [];
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
        function val = funCallback(x, muc, cholc, nbEvals, lambda, varargins)
            if(nargin < 6)
                warning('this call sould not happen during the core of the optimization');
                [f, ~] = muc{:};
                val = -f.eval(x');
                return;
            end
            % evaluation
            [f, valStorage] = varargins{:};
            val = -f.eval(x');
            
            % perform special operations (e.g. plotting) after each lambda
            % points (each generation)
            valStorage.allEvals = [valStorage.allEvals; val];
            valStorage.sampledPoints = [valStorage.sampledPoints x];
            if(mod(nbEvals + 1, lambda) == 0)
                valStorage.endOfGenerationCallBack(muc, cholc, lambda);
            end
        end
    end
end
