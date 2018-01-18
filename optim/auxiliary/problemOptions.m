%% author: Sebastian Rinder

classdef problemOptions < handle
    properties
        maxExp;
        opts;
    end
    
    methods
        function obj = problemOptions(kernelName, platform, environment)    
            obj.maxExp = 100;
            obj.opts = environmentSettings(environment, platform, 'none');
            obj.opts.platform = platform;
            
            if strcmp(kernelName, 'sexp')
                obj.opts.distanceMat = @obj.euclideanDistance;
                obj.opts.scaleKernel = @obj.sexpScaleKernel;
            elseif strcmp(kernelName, 'matern52')
                obj.opts.distanceMat = @obj.euclideanDistance;
                obj.opts.scaleKernel = @obj.matern52ScaleKernel;
            elseif strcmp(kernelName, 'trajectory')
                obj.opts.distanceMat = @obj.trajectoryDistance;
                obj.opts.scaleKernel = @obj.trajectoryScaleKernel;
            else
                error(['no covariance function for ', kernelName]);
            end
        end
        
        function D = euclideanDistance(obj, Xm, Xn, trajectories, isCovMat, opts)
            D = euclideanDistance(Xm, Xn);
        end
        
        function D = trajectoryDistance(obj, Xm, Xn, trajectories, isCovMat, opts)
            D = trajectoryDistance(Xm, Xn, trajectories, isCovMat, opts);
        end        
      
        function [K,norml] = sexpScaleKernel(obj, D, hyper)
            if nargout == 1
                sigmaf = exp(hyper(1,1));
                sigmal = exp(hyper(1,2));
                
                K = sigmaf.^2 .* exp(-0.5 .* (D./sigmal).^2);
            else
                K = [];
                norml = max(max(D)) / sqrt(2*obj.maxExp);
                norml = max(norml,1);
            end
        end
        
        function [K,norml] = matern52ScaleKernel(obj, D, hyper)
            if nargout == 1
                sigmaf = exp(hyper(1,1));
                sigmal = exp(hyper(1,2));

                K1 = sqrt(5) .* D ./ sigmal;
                K2 = 5 .* D .^ 2 / (3.*sigmal^2);
                K3 = -sqrt(5) .* D ./ sigmal;
                K = sigmaf .^ 2 .* (1 + K1 + K2) .* exp(K3);
            else
                K = [];                
                norml = sqrt(5) * max(max(D)) / obj.maxExp;
                norml = max(norml,1);
            end
        end
        
        function [K,norml] = trajectoryScaleKernel(obj, D, hyper)
            if nargout == 1
                sigmaf = exp(hyper(1,1));
                sigmal = exp(hyper(1,2));

                K = sigmaf .* exp (-D ./sigmal);
            else
                K = [];                
                norml = max(max(D)) / obj.maxExp;
                norml = max(norml,1);
            end
        end
        
        function [vals, traj] = eval(obj, samples)
            if isequal(obj.opts.platform, 'matlab')
                vals = zeros(size(samples,1),1);
                traj = cell(size(samples,1),1);
                for i=1:size(samples,1)
                    for j=1:obj.opts.trajectoriesPerSample
                        [vals(i,j), traj{i,j}] = objectiveFcn(samples(i,:), obj.opts);
                    end
                end
            elseif isequal(obj.opts.platform, 'pygym')
                %TODO
            end
        end
    end
    
end