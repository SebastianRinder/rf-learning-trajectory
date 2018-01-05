classdef RandomGaussianMixtureFunction < handle
    properties
        normals;
        signature;
    end
    
    methods
        function obj = RandomGaussianMixtureFunction(nbGauss, rangeMu, minVar, dim, seed)
            if(~exist('dim', 'var') || isempty(dim))
                dim = 2;
            end
            if(exist('seed', 'var'))
                rng(seed);
            end
            
            seed = rng;
            obj.signature = ['func' num2str(seed.Seed) '_' num2str(nbGauss) '_' num2str(rangeMu) '_' num2str(minVar)];
            
            obj.normals = {};
            for k = 1:nbGauss
                covC = ones(dim, dim) - 2 * rand(dim, dim);
                covC = covC + diag(minVar * ones(dim,1));
                covC = covC' * covC;
                mu = rangeMu * (ones(1,dim) - 2 * rand(1,dim));
                obj.normals{end+1} = Normal(mu, covC);
            end
        end
        
        function str = getSignature(obj)
            str = obj.signature;
        end
        
        function plot(obj)
            %only valid if dim == 2
            [limitInf, limitSup] = obj.getRange();
            nbPoint = 99;
            line1 = limitInf(1):(limitSup(1)-limitInf(1))/nbPoint:limitSup(1);
            line2 = limitInf(2):(limitSup(2)-limitInf(2))/nbPoint:limitSup(2);
            [x1, x2] = meshgrid(line1, line2);
            y = obj.eval([x1(:) x2(:)]);
            %h = surf(line1, line2, reshape(y, length(line1), length(line2))');
            %set(h,'LineStyle','none');
            contour(line1, line2, reshape(y, length(line2), length(line1)));
        end
        
        function [limitInf, limitSup] = getRange(obj)
            infLimits = zeros(length(obj.normals), obj.normals{1}.getDim);
            supLimits = infLimits;
            for k = 1:length(obj.normals)
                [infLimits(k, :), supLimits(k, :)] = obj.normals{k}.getLimits;
            end
            
            limitInf = min(infLimits, [], 1);
            limitSup = max(supLimits, [], 1);
        end
        
        function plotNormals(obj)
            for k = 1:length(obj.normals)
                obj.normals{k}.plot();
            end
        end
        
        
        function vals = eval(obj, samples)
            vals = zeros(size(samples, 1), 1);
            for k = 1:length(obj.normals)
                vals = vals + exp(obj.normals{k}.getLogProbas(samples));
            end
        end
    end
end

