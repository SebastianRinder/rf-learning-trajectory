classdef IProjection
    methods (Static)
        
        function [y, dy] = eval_chol(params, initPrecision, initMu, thetas, rewards, logWeights, eta, doNormalized)
            if(~exist('doNormalized', 'var'))
                doNormalized = true;
            end

            % Extracting the mu and the precision matrix
            dim = length(initMu);
            mu = params(1:dim)';
            cholP = zeros(dim);
            cholP(triu(ones(dim)) == 1) = params(dim+1:end);
            precision = cholP' * cholP;
            logDetP = 2 * sum(log(diag(cholP)));
            invCholP = cholP \ eye(dim);
            covc = invCholP *  invCholP';
            
            % Weighted sum
            if(doNormalized)
                [wgSum, parDeriv, dcovRaw] = static_optimization_algs.IProjection.WeightedGaussPrecisionNormalized(...
                    mu, precision, logDetP, thetas, rewards, logWeights);
            else
                [wgSum, parDeriv, dcovRaw] = static_optimization_algs.IProjection.WeightedGaussPrecision(...
                    mu, precision, logDetP, covc, thetas, rewards, logWeights);
            end
            
            y = .5 * (logDetP + trace(initPrecision * covc) +  ...
                (mu-initMu) * initPrecision * (mu-initMu)') - wgSum / eta;
            
            if(nargout > 1)                
                % Gradient Precision
                dy = zeros(length(params), 1);
                dy(1:dim) = initPrecision * (mu - initMu)' - parDeriv(1:dim) / eta;
                dcov = covc * (precision - initPrecision) * covc - dcovRaw / eta;
                dchol = cholP * dcov;
                dy(dim+1:end) = dchol(triu(ones(dim)) == 1);
            end
        end
        
        function [y, dy] = eval(params, initPrecision, initMu, thetas, rewards, logWeights, eta, doNormalized)
            % 
            if(~exist('doNormalized', 'var')) 
                doNormalized = true;
            end
            % Extracting the mu and the precision matrix
            dim = length(initMu);
            mu = params(1:dim)';
            precision = zeros(dim);
            precision(tril(ones(dim)) == 1) = params(dim+1:end);
            precision = precision + precision' - diag(diag(precision));
            [cholP, p] = chol(precision);
            if(p > 0)
                y = inf;
                dy = zeros(length(params), 1);
                return;
            end
            detP = prod(diag(cholP))^2;
            invCholP = cholP \ eye(dim);
            covc = invCholP *  invCholP';
            
            % Weighted sum
            if(doNormalized)
                [wgSum, parDeriv] = static_optimization_algs.IProjection.WeightedGaussPrecisionNormalized(...
                    mu, precision, detP, thetas, rewards, logWeights);
            else
                [wgSum, parDeriv] = static_optimization_algs.IProjection.WeightedGaussPrecision(...
                    mu, precision, detP, covc, thetas, rewards, logWeights);
            end
            
            y = .5 * (log(detP) + trace(initPrecision * covc) +  ...
                (mu-initMu) * initPrecision * (mu-initMu)') - wgSum / eta;
            
            if(nargout > 1)
                % Gradient Precision
                dy = zeros(length(params), 1);
                dy(1:dim) = initPrecision * (mu - initMu)';
                dcov = covc * (precision - initPrecision) * covc;
                dcov = dcov - .5 * dcov .* eye(dim);
                dy(dim+1:end) = dcov(tril(ones(dim)) == 1);
                dy = dy - parDeriv / eta;
            end
            
            if(any(~isfinite(dy)) || ~isreal(dy))
                warning('you''re in troubles.');
            end
            
        end
        
        function [y, dy, dcovRaw] = WeightedGaussPrecisionNormalized(mu, precision, logDetP, thetas, rewards, logWeights)
            % Weighted sum
            dim = length(mu);
            logProbas = static_optimization_algs.Normal.getLogProbasPrecisionLogDet(thetas, mu, precision, logDetP);
            probaRatios = exp(logProbas - logWeights);            
            infProbas = isinf(probaRatios);
            if(any(infProbas))
                probaRatios = zeros(size(probaRatios));
                probaRatios(infProbas) = 1;
            else
                maxProbas = max(probaRatios);
                if(~isfinite(probaRatios(1) / maxProbas)) % division by zero
                    probaRatios = ones(size(probaRatios)) / length(probaRatios);                    
                else
                    probaRatios = probaRatios / maxProbas;
                end
            end
            
            weightedRewards = probaRatios .* rewards;
            sumToN = sum(probaRatios);
            
            y = sum(weightedRewards) / sumToN;
            
            if(nargout > 1)                
                % Gradient Mu
                dy = zeros(dim + (dim * (dim + 1))/2, 1);
                gradMu = bsxfun(@times, thetas, weightedRewards)...
                    - bsxfun(@times, thetas, probaRatios) * y;
                dy(1:dim) =  (sum(gradMu, 1) / sumToN) * precision;
                
                % Gradient Precision
                diff = bsxfun(@minus, thetas, mu);
                dcov = bsxfun(@times, diff, y * probaRatios - weightedRewards)' * diff / sumToN;
                dcovRaw = dcov;
                dcov = dcov - .5 * dcov .* eye(dim);
                dy(dim+1:end) = dcov(tril(ones(dim)) == 1);
            end
        end
        
        function [y, dy, dcovRaw] = WeightedGaussPrecision(mu, precision, detP, covc, thetas, rewards, logWeights)
            % Weighted sum
            dim = length(mu);
            logProbas = static_optimization_algs.Normal.getLogProbasPrecision(thetas, mu, precision, detP);
            allWeights = exp(logProbas - logWeights) .* rewards;
            y = sum(allWeights) / length(logWeights);
            
            if(nargout > 1)
                % Gradient Mu
                dy = zeros(dim + (dim * (dim + 1))/2, 1);
                diff = bsxfun(@minus, thetas, mu);
                gradMu = diff * precision;
                gradMu = bsxfun(@times, gradMu, allWeights);
                dy(1:dim) = sum(gradMu, 1) / length(logWeights);
                
                % Gradient Precision
                dcov = y * covc - bsxfun(@times, diff, allWeights)' * diff / length(logWeights);
                dcovRaw = dcov;
                dcov = dcov - .5 * dcov .* eye(dim);
                dy(dim+1:end) = dcov(tril(ones(dim)) == 1);
            end
        end
        
        % obsolete 
        function [y, dy] = WeightedGauss(params, dim, thetas, weights)
            % Extracting the mu and the precision matrix
            mu = params(1:dim)';
            precision = zeros(dim);
            precision(tril(ones(dim)) == 1) = params(dim+1:end);
            precision = precision + precision' - diag(diag(precision));
            detP = det(precision);
            if(detP <= 0)
                y = inf;
                dy = zeros(length(params), 1);
                return;
            end
            
            % Weighted sum
            probas = exp(static_optimization_algs.Normal.getLogProbasPrecision(thetas, mu, precision, detP));
            y = -probas' * weights;
            
            % Gradient Mu
            dy = zeros(length(params), 1);
            gradMu = bsxfun(@minus, thetas, mu) * precision;
            gradMu = bsxfun(@times, gradMu, probas .* weights);
            dy(1:dim) = sum(gradMu, 1);
            
            % Gradient Precision
            covc = inv(precision);
            dcov = zeros(dim, dim);
            for k = 1:length(weights) % can be done with bsxfun, see above
                dcov = dcov + .5 * probas(k) * weights(k) * (covc - (thetas(k, :) - mu)' * (thetas(k, :) - mu));
            end
            dcov = 2 * dcov - dcov .* eye(dim);
            dy(dim+1:end) = dcov(tril(ones(dim)) == 1);
            dy = -dy;
        end
        
        function [y, dy] = QuadFunc(a, r)
            A = reshape(a, sqrt(length(a)), sqrt(length(a)));
            y = r' * A * r;
            dy = reshape(r * r', length(a), 1);
        end
    end
end

