classdef CMAESWrapper < handle
    properties
        opts;
        fun;
        xstart;
        insigma;
    end
    
    methods
        %%
        function obj = CMAESWrapper(fun, optimizerInput)
            obj.fun = fun;
            opts = cmaes('defaults');
            opts.EvalInitialX = 'no';
            opts.Restarts = '0';
            opts.CMA.active = '0';
            opts.DiagonalOnly = '0';
            opts.MaxFunEvals = optimizerInput.maxEvals;
            opts.TolFun = 1e-3;
            obj.opts = opts;
            obj.xstart = optimizerInput.xstart;
            obj.insigma = optimizerInput.insigma;

        end    
       
        function sign = getSignature(obj)
            sign = ['CMAESWrapper' '_' num2str(obj.insigma)];
        end
        
        %% main function
        function [xmin, fxmin] = optimize(obj)
            [xmin, fxmin] = cmaes(obj.fun, obj.xstart, obj.insigma, obj.opts);
        end
    end
end
