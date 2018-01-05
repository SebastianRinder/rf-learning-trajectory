classdef Normal
    %NORMAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        mu;
        cholC;
        cholP;
        invC; % precision matrix
    end
    
    methods
        function obj = Normal(mu, covC, matType)
            if(size(mu, 2) == length(mu))
                obj.mu = mu;
            else
                obj.mu = mu';
            end
            
            if(exist('matType', 'var'))
                if(matType == 1) %%% covC is cholesky of the precision
                    obj.cholP = covC;
                    obj.cholC = (covC\eye(length(mu)))';
                    obj.invC = covC' * covC;
                elseif(matType == 2) %%% covC is the precision
                    obj.invC = covC;
                    obj.cholP = chol(obj.invC);
                    obj.cholC = (obj.cholP \ eye(length(obj.mu)))';
                else
                    error('matType unrecognized');
                end
            else
                obj.cholC = chol(covC);
                obj.cholP = (obj.cholC\eye(length(mu)))';
                obj.invC = obj.cholP' * obj.cholP;
            end
        end
        
        function samples = getSamples(obj, numSamples)
            samples = randn(numSamples, obj.getDim) * obj.cholC + repmat(obj.mu, numSamples, 1);
        end
        
        function probas = getLogProbas(obj, samples)
%             diff = bsxfun(@minus, samples, obj.mu);
%             probas = dot(diff * obj.invC, diff, 2);
%             probas = -.5 * bsxfun(@plus, probas, obj.getDim * log(2*pi) + log(prod(diag(obj.cholC))^2));
            
            diff = bsxfun(@minus, samples, obj.mu) / obj.cholC;
            probas = dot(diff, diff , 2);
            probas = -.5 * (probas + obj.getDim * log(2*pi)) - sum(log(diag(obj.cholC)));

        end
        
        function probas = getLogProbasP(obj, samples)
            diff = bsxfun(@minus, samples, obj.mu);
            probas = dot(diff * obj.invC, diff, 2);
            probas = -.5 * (probas + (length(obj.mu) * log(2*pi) - obj.getLogDetPrecision) * ones(size(probas)));
        end

        
        % KL(obj||q)
        function kl = getKL(obj, q) 
            precisionOld = q.getPrecision;
            logDetPOld = q.getLogDetPrecision;
            muOld = q.mu;
            logDetP = obj.getLogDetPrecision;
            cov = obj.getCov;            
            kl = .5 * (trace(precisionOld * cov) + (obj.mu-muOld) * precisionOld * (obj.mu-muOld)' - length(muOld) - logDetPOld + logDetP);
        end
        
        
        % Entropy
        function h = getEntropy(obj) 
             h = .5 * (length(obj.mu) * log(2 * pi * exp(1)) - obj.getLogDetPrecision);
        end
        
        function obj = wmle(obj, samples, weights)
            if(~exist('weights', 'var'))
                weights = ones(size(samples, 1), 1) / size(samples, 1);
            end                   
            [z, weights] = static_optimization_algs.Normal.normalizeWeights(weights);
            %%%mean estimation
            obj.mu = sum(bsxfun(@times, samples, weights));
            
            %%%covariance estimation
            if(z > 1e-3)
                diff = bsxfun(@minus, samples, obj.mu);
                covC = bsxfun(@times, diff, weights)' * diff / z;
                regul = 1e-20;
                p = 1;
                while(p > 0)
                    [obj.cholC, p] = chol(covC + regul * eye(size(covC)));
                    regul = regul * 10;
                end
            else
                warning('only one sample')
                obj.cholC = eye(size(obj.cholC)) * 1e-16;
            end            
            obj.cholP = (obj.cholC \ eye(length(obj.mu)))';
            obj.invC = obj.cholP' * obj.cholP;
        end
        
        function covC = getCovariance(obj)
            covC = obj.cholC' * obj.cholC;
        end
        
        function plot(obj, nbPoint)
            % for dim = 2 only
            if(obj.getDim == 2)
                if(~exist('nbPoint', 'var'))
                    nbPoint = 50;
                end
                [limitInf, limitSup] = obj.getLimits;
                line1 = limitInf(1):(limitSup(1)-limitInf(1))/nbPoint:limitSup(1);
                line2 = limitInf(2):(limitSup(2)-limitInf(2))/nbPoint:limitSup(2);
                [x1, x2] = meshgrid(line1, line2);
                y = exp(obj.getLogProbas([x1(:) x2(:)]));
                %                 h = surf(line1, line2, reshape(y, length(line2), length(line1)));
                %                 set(h,'LineStyle','none');
                contour(line1, line2, reshape(y, length(line2), length(line1)));
            else
                warning('plot not implemented for univariate gaussian');
            end
        end
        
        function [limInf, limSup] = getLimits(obj)
            % for plotting
            sigma = diag(obj.getCovariance)'.^.5;
            sigma = max(sigma, 0.1 * ones(size(sigma)));
            limInf = obj.mu - sigma * 2;
            limSup = obj.mu + sigma * 2;
        end
        
        function dim = getDim(obj)
            dim = length(obj.mu);
        end
        
        function mu = getMu(obj)
            mu = obj.mu;
        end
        
        function precision = getPrecision(obj)
            precision = obj.invC;
        end
        
        function covRet = getCov(obj)
            covRet = obj.cholC' * obj.cholC;
        end
        
        function detCov = getDetCov(obj)
            detCov = prod(diag(obj.cholC))^2;
        end
        
        function cholC = getCholC(obj)
            cholC = obj.cholC;
        end
        
        function cholP = getCholP(obj)
            cholP = obj.cholP;
        end
        
        function logDetPrecision = getLogDetPrecision(obj)
            logDetPrecision = 2 * sum(log(diag(obj.cholP)));
        end

        
        function detPrecision = getDetPrecision(obj)
            detPrecision = 1 / obj.getDetCov();
        end
        
        function obj = setMu(obj, mu)
            obj.mu = mu;
        end
        
        function obj = setPrecision(obj, P)
            obj.invC = P;
            obj.cholP = chol(P);
            obj.cholC = (obj.cholP \ eye(length(obj.mu)))';
        end

    end
    
    methods (Static)
        function probas = getLogProbasPrecision(samples, mu, precision, detP)
            diff = bsxfun(@minus, samples, mu);
            probas = dot(diff * precision, diff, 2);
            probas = -.5 * (probas + (length(mu) * log(2*pi) - log(detP)) * ones(size(probas)));
        end

        function probas = getLogProbasPrecisionLogDet(samples, mu, precision, logDetP)
            diff = bsxfun(@minus, samples, mu);
            probas = dot(diff * precision, diff, 2);
            probas = -.5 * (probas + (length(mu) * log(2*pi) - logDetP) * ones(size(probas)));
        end

        
        %%% u ~ N(X * K, covC)
        function [K, covC] = wmleLinearMean(X, U, weights, regularization)
            if(~exist('weights', 'var') || isempty(weights))
                weights = ones(size(X, 1), 1) / size(X, 1);
            end  
            if(~exist('regularization', 'var') || isempty(regularization))
                regularization = 0;
            end 
            [z, weights] = static_optimization_algs.Normal.normalizeWeights(weights);
            
            %%% mean estimation
            DhalfX = bsxfun(@times, X, sqrt(weights));
            cholM = chol((DhalfX' * DhalfX) + regularization * eye(size(X, 2)));
            K = cholM \ (bsxfun(@times, X, weights) / cholM)' * U;

            %%% covariance
            if(nargout > 1)
                if(z > 1e-3)
                    diff = U - X * K;
                    covC = bsxfun(@times, diff, weights)' * diff / z;
                else
                    warning('only one sample')
                    covC = eye(size(U, 2)) * 1e-16;
                end
            end
        end
        
        
        function [z, weights] = normalizeWeights(weights)
            sumW = sum(weights);
            if(sumW == 0)
                weights = ones(size(samples, 1), 1) / size(samples, 1);
            else
                weights = weights / sumW;
            end
            z = 1 - sum(weights.^2);
        end
    end
end

