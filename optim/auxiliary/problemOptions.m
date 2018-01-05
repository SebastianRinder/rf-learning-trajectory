%% author: Sebastian Rinder

classdef problemOptions < handle
    properties
        opts;
    end
    
    methods
        function obj = problemOptions(kernelName, environment)            
            obj.opts = environmentSettings(environment, 'none');
            
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
        
        function D = euclideanDistance(obj, Xm, Xn, ~, ~, ~)
            m = size(Xm,1);
            n = size(Xn,1);

            D = zeros(m,n);
            
            if m==n && all(all(Xm == Xn))
                for i=1:size(Xm,1)
                    D(i,i) = 0;
                    k = i+1;
                    for j=k:size(Xn,1)
                        r = Xm(i,:) - Xn(j,:);
                        D(i,j) = r*r';
                    end
                end
                D = sqrt(D + D');
            else
                for i=1:size(Xm,1)
                    for j=1:size(Xn,1)
                        r = Xm(i,:) - Xn(j,:);
                        D(i,j) = r*r';
                    end
                end
                D = sqrt(D);
            end

            %pdist2(Xm,Xn); %licence error on cluster
        end
        
        function D = trajectoryDistance(obj, Xm, Xn, trajectories, isCovMat, opts)
            D = trajectoryDistance(Xm, Xn, trajectories, isCovMat, opts);
        end        
      
        function K = sexpScaleKernel(obj, D, hyper)
            sigmaf = exp(hyper(1,1));
            sigmal = exp(hyper(1,2));
            
            K = sigmaf.^2 .* exp(-0.5 .* (D./sigmal).^2);
        end
        
        function K = matern52ScaleKernel(obj, D, hyper)
            sigmaf = exp(hyper(1,1));
            sigmal = exp(hyper(1,2));
            
            K1 = sqrt(5) .* D ./ sigmal;
            K2 = 5 .* D .^ 2 / (3.*sigmal^2);
            K3 = -sqrt(5) .* D ./ sigmal;
            K = sigmaf .^ 2 .* (1 + K1 + K2) .* exp(K3);
        end
        
        function K = trajectoryScaleKernel(obj, D, hyper)
            sigmaf = exp(hyper(1,1));
            sigmal = exp(hyper(1,2));
            
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