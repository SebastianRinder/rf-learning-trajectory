classdef gpOptions < handle
    properties
        opts;
    end
    
    methods
        function obj = gpOptions(kernelName, environment)            
            obj.opts = environmentSettings(environment, 'none');
            
            if strcmp(kernelName, 'sexp')
                obj.opts.distanceMat = @obj.squaredDist;
                obj.opts.scaleKernel = @obj.sexpScaleKernel;
            elseif strcmp(kernelName, 'matern52')
                obj.opts.distanceMat = @obj.squaredDist;
                obj.opts.scaleKernel = @obj.matern52ScaleKernel;
            elseif strcmp(kernelName, 'trajectory')
                obj.opts.distanceMat = trajectoryDistance;
                obj.opts.scaleKernel = @obj.trajectoryScaleKernel;
            else
                error(['no covariance function for ', kernelName]);
            end
        end
        
        
        function [sigmaf, sigmal] = expHyper(obj, hyper)
            sigmaf = exp(hyper(1,1));
            sigmal = exp(hyper(1,2));
        end
        
        function D = squaredDist(obj, Xm, Xn, ~, ~, ~)
            D = pdist2(Xm,Xn).^2;
        end
        
        function K = sexpScaleKernel(obj, D, hyper)
            [sigmaf, sigmal] = expHyper(obj, hyper);
            
            K = sigmaf.^2 .* exp(-D./(2.*sigmal.^2));
        end
        
        function K = matern52ScaleKernel(obj, D, hyper)
            [sigmaf, sigmal] = expHyper(obj, hyper);
            
            K1 = sqrt(5) .* D ./ sigmal;
            K2 = 5 .* D .^ 2 / (3.*sigmal^2);
            K3 = -sqrt(5) .* D ./ sigmal;
            K = sigmaf .^ 2 .* (1 + K1 + K2) .* exp(K3);
        end
        
        function K = trajectoryScaleKernel(obj, D, hyper)
            [sigmaf, sigmal] = expHyper(obj, hyper);
            
            K = sigmaf .* exp (-D .*sigmal);
        end
        
        function [vals, traj] = eval(obj, samples)
            vals = zeros(size(samples,1),1);
            traj = cell(size(samples,1),1);
            for i=1:size(samples,1)
                [vals(i,1), traj{i,1}] = objectiveFcn(samples(i,:), obj.opts);
            end
        end
    end
    
end