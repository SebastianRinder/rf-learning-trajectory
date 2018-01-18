classdef Quadratic 
    %%%
    %%% quadModel contains A, w and b
    %%%
    methods (Static)
        function params = packParams(quadModel, matIndices)            
            dim = size(quadModel.A, 1);            
            if(~exist('matIndices', 'var') || isempty(matIndices))
                %default behavior is to take the upper triangular part of the matrix
                matIndices = find(triu(ones(dim))); 
            end
            lengthParams = 1 + length(matIndices);
            if(isfield(quadModel, 'w') && ~isempty(quadModel.w))
                %has linear part 
                lengthParams = lengthParams + dim;
            end
            params = zeros(lengthParams, 1);
            params(1) = quadModel.b;
            if(isfield(quadModel, 'w') && ~isempty(quadModel.w))
                params(2:dim+1) = quadModel.w;
            end
            params((end-length(matIndices)+1):end) = quadModel.A(matIndices);
        end
        
        %%
        function quadModel = unpackParams(params, dim, matIndices, isSymmetric)
            if(~exist('matIndices', 'var') || isempty(matIndices))
                %default behavior is to take the upper triangular part of the matrix
                matIndices = find(triu(ones(dim))); 
            end
            quadModel.b = params(1);
            if(length(params) > 1 + length(matIndices))
                %has linear part
                quadModel.w = params(2:dim+1);            
            end
            A = zeros(dim);
            A(matIndices) = params((end-length(matIndices)+1):end);     
            if(exist('isSymmetric', 'var') && isSymmetric) %assuming f(x) = 1/2 x' A x
                A = A + A';
            end
            quadModel.A = A;
        end
        
        %%
        function phix = getQuadFeatures(x, addBias, addLinear)
            if(~exist('addBias', 'var') || isempty(addBias))
                addBias = 1;
            end
            if(~exist('addLinear', 'var') || isempty(addLinear))
                addLinear = 1;
            end
            quadPart = reshape(bsxfun(@times, permute(x, [2 3 1]), permute(x, [3 2 1])), [], size(x,1))';
            keep = nonzeros(triu(reshape(1:size(x,2)^2, size(x,2), size(x,2))));
            if(addBias && addLinear)
                phix = [ones(size(x,1), 1) x quadPart(:, keep)];
            elseif(addLinear)
                phix = [x quadPart(:, keep)];
            elseif(addBias)
                phix = [ones(size(x,1), 1) quadPart(:, keep)];
            else
                phix = quadPart(:, keep);
            end
        end
        
        %% trim dataset and only keep trimRatio * 100% of the weight mass
        function [x, y, weights] = trimDataSet(x, y, weights, trimRatio)
            weights = static_optimization_algs.Quadratic.normalizeWeights(weights);
            [sWeights, is] = sort(weights, 'descend');
            firstId = find(cumsum(sWeights / sum(sWeights)) > trimRatio, 1, 'first');
            if(~isempty(firstId) && firstId < length(is))
                is = is(1:firstId);
                x = x(is, :);
                y = y(is, :);
                weights = sWeights(1:firstId);
            end
            weights = static_optimization_algs.Quadratic.normalizeWeights(weights); %renormalize weights
        end
        
        %% learn mean squared error quad model f(x) = 1/2 x' A x + w'x + b
        function [thetas, quadModel, x, phix, y, weights, regularization] = learnQuadModel(x, y, regularization, weights, trimRatio)
            if(~exist('weights', 'var') || isempty(weights))
                weights = ones(size(x, 1), 1) / size(x, 1);
            end
            if(exist('trimRatio', 'var') && ~isempty(trimRatio))
                [x, y, weights] = static_optimization_algs.Quadratic.trimDataSet(x, y, weights, trimRatio);
            end
            phix = static_optimization_algs.Quadratic.getQuadFeatures(x, 1); %1: has bias
            if(isstruct(regularization)) %do cross validation            
                regularization = static_optimization_algs.Quadratic.hyperParamOptim...
                    (phix, y, weights, regularization.hyperParamList, regularization.validationSetRatio, regularization.nbCrossValidation);
            end
            thetas = static_optimization_algs.Quadratic.MSEEstimator(phix, y, weights, regularization);
            quadModel = static_optimization_algs.Quadratic.unpackParams(thetas, size(x, 2), [], 1);
        end
        
        %% learn mean squared error quad model f(x) = 1/2 x' A x + w'x + b, where A is n.d.
        function [thetas, quadModel, phix, y, weights] = learnConcaveQuadModel(x, y, regularization, weights, trimRatio, additionalGradientSteps)
            if(~exist('weights', 'var'))
                weights = [];
            end
            if(~exist('trimRatio', 'var'))
                trimRatio = [];
            end
            [thetas, quadModel, x, phix, y, weights, regularization] = static_optimization_algs.Quadratic.learnQuadModel(x, y, regularization, weights, trimRatio);
            %%% Ensure that A is negative definite (theoretically n.s.d. is enough, but chol is much faster than eig on matlab)
            [isNd, quadModel.A] = static_optimization_algs.Quadratic.isNegativeDefinite(quadModel.A);
            if(~isNd)                         
                %%% relearn linear and bias
                linTarget = y - .5 * dot(x * quadModel.A, x, 2); %get rid of the quad part
                [quadModel.w, quadModel.b] = static_optimization_algs.Quadratic.MSEEstimator(x, linTarget, weights, regularization);
                %%% Improve over thetas by gradient descent
                if(exist('additionalGradientSteps', 'var') && additionalGradientSteps > 0)
                    [~, quadModel] = static_optimization_algs.Quadratic.optimizeThetas...
                        (quadModel, size(x, 2), phix, y, regularization, weights, additionalGradientSteps);
                end
                thetas = static_optimization_algs.Quadratic.packParams(quadModel, []);
            end
        end
        
        %% hyper-param optimization: monte carlo cross validation
        function [hyperParam] = hyperParamOptim(x, y, weights, hyperParamList, validRatio, nbCv)
            validScores = zeros(nbCv, length(hyperParamList));
            for icv = 1:nbCv
                shuffledI = randperm(length(y));
                validI = shuffledI(1:round(validRatio * length(y)));
                learnI = shuffledI(round(validRatio * length(y))+1:end);
                for ihp = 1:length(hyperParamList)
                    [~, ~, validScores(icv, ihp)] = static_optimization_algs.Quadratic.MSEEstimator...
                        (x(learnI, :), y(learnI, :), weights(learnI, :), hyperParamList(ihp),...
                        x(validI, :), y(validI, :), weights(validI, :));
                end
            end
            [~, ihpStar] = min(mean(validScores, 1));
            hyperParam = hyperParamList(ihpStar);
        end
        
        %% check if A is n.d. if not, project to the space of n.d. matrices
        function [isNd, An] = isNegativeDefinite(A)
            [~, notPsd] = chol(-A);
            isNd = ~notPsd;
            if(nargout > 1)
                if(notPsd)
                    [V, D]  = eig(A);
                    %%% Project A to the space of n.s.d matrices
                    negativeEigs = find(diag(D) < 0);
                    A = V(:, negativeEigs) * D(negativeEigs, negativeEigs) * V(:, negativeEigs)';                
                    
                    %%% substract small epsi until n.d.
                    epsiSub = 1e-8;
                    An = A;
                    [~, notPsd] = chol(-An);
                    while(notPsd)
                        An = A - epsiSub * eye(size(A));
                        [~, notPsd] = chol(-An);
                        epsiSub = epsiSub * 10;
                    end
                else
                    An = A;
                end
            end
        end
        
        %% minimize concaveWeightedMeanSquare function
        function [thetas, quadModel] = optimizeThetas(quadModel0, dim, phix, y, regularization, weights, maxFunEval)
            matIndices = find(triu(ones(dim)));
            quadModel0.A = chol(-quadModel0.A);
            thetas0 = static_optimization_algs.Quadratic.packParams(quadModel0, matIndices);
            fun = @(params) static_optimization_algs.Quadratic.concaveWeightedMeanSquareChol(params, dim, phix, y, regularization, weights, matIndices);
            optOptim =  optimset('GradObj','on', 'DerivativeCheck', 'off', 'Display', 'off');
            optOptim.MaxFunEvals = maxFunEval;
            thetas = fminunc(fun, thetas0, optOptim);
            if(nargout > 1)
                quadModel = static_optimization_algs.Quadratic.unpackParams(thetas, dim, matIndices, 1);                
            end
        end
        
        %% Computes the weighted mean square error and its gradient. If the function is not concave the function returns +inf.
        function [g, dg] = concaveWeightedMeanSquareChol(thetas, dim, phix, y, regularization, weights, matIndices)
            quadModel = static_optimization_algs.Quadratic.unpackParams(thetas, dim, matIndices, 0);
            cholA = quadModel.A;
            quadModel.A = -cholA' * cholA;
            thetas = static_optimization_algs.Quadratic.packParams(quadModel, matIndices);
            [g, dg] = static_optimization_algs.Quadratic.concaveWeightedMeanSquare(thetas, dim, phix, y, regularization, weights, matIndices, 0);
            gradQuad = static_optimization_algs.Quadratic.unpackParams(dg, dim, matIndices, 1);
            gradQuad.A = -cholA * (gradQuad.A);
            dg = static_optimization_algs.Quadratic.packParams(gradQuad, matIndices);
        end

        
        %% Computes the weighted mean square error and its gradient. If the function is not concave the function returns +inf.
        function [g, dg] = concaveWeightedMeanSquare(thetas, dim, phix, y, regularization, weights, matIndices, checkPsd)
            if(checkPsd)
                quadModel = static_optimization_algs.Quadratic.unpackParams(thetas, dim, matIndices, 1);
                if(~static_optimization_algs.Quadratic.isNegativeDefinite(quadModel.A))
                    g = inf;
                    dg = 0;
                    return;
                end
            end            
            % compute regularized mean-squared error
            g = sum(weights .* ((phix * thetas - y).^2)) / size(phix, 1) + regularization * (thetas' * thetas) ;
            if(nargout > 1)
                dg = 2 * (sum(bsxfun(@times, phix, weights .* (phix * thetas - y)))' / size(phix, 1) + regularization * thetas);
            end
        end
        
        %%
        function plotQuadParams(params, center, width)
            dim = length(center);
            if(dim == 2)
                nbPoint = 100;
                limitInf = center - width;
                limitSup = center + width;
                line1 = limitInf(1):(limitSup(1)-limitInf(1))/nbPoint:limitSup(1);
                line2 = limitInf(2):(limitSup(2)-limitInf(2))/nbPoint:limitSup(2);
                [x1, x2] = meshgrid(line1, line2);
                phix = static_optimization_algs.Quadratic.getQuadFeatures([x1(:) x2(:)], 1);
                y = phix * params;
                %                 h = surf(line1, line2, reshape(y, length(line2), length(line1)));
                %                 set(h,'LineStyle','none');
                contour(line1, line2, reshape(y, length(line2), length(line1)));
            else
                warning('plot not implemented for this dimension');
            end
        end
        %% 
        function y = evaluate(quadModel, X)
            y = .5 * dot(X, X * quadModel.A, 2) + X * quadModel.w + quadModel.b;
        end
        
        %%
        function plotQuad(quadModel, center, width)
            dim = length(center);
            if(dim == 2)
                nbPoint = 100;
                limitInf = center - width;
                limitSup = center + width;
                line1 = limitInf(1):(limitSup(1)-limitInf(1))/nbPoint:limitSup(1);
                line2 = limitInf(2):(limitSup(2)-limitInf(2))/nbPoint:limitSup(2);
                [x1, x2] = meshgrid(line1, line2);
                x = [x1(:) x2(:)];
                y = static_optimization_algs.Quadratic.evaluate(quadModel, x);
                %                 h = surf(line1, line2, reshape(y, length(line2), length(line1)));
                %                 set(h,'LineStyle','none');
                contour(line1, line2, reshape(y, length(line2), length(line1)));
            else
                warning('plot not implemented for this dimension');
            end
            
        end
        
        function weights = normalizeWeights(weights)
            infWeights = isinf(weights);
            if(any(infWeights))
                weights = zeros(size(weights));
                weights(infWeights) = 1;
            else
                maxWeights = max(weights);
                if(~isfinite(weights(1) / maxWeights)) % division by zero
                    weights = ones(size(weights)) / length(weights);
                else
                    weights = weights / maxWeights;
                end
            end            
        end
        
        % Ridge Regression. Last 3 params only necessary for validation
        % score computation
        function [param, bias, validationScore] = MSEEstimator(X, y, weights, regul, vX, vy, vweights)
            % normalize weights
            weights = static_optimization_algs.Quadratic.normalizeWeights(weights);
            
            % standardize data
            rangeInput = std(X,[], 1);                        
            meanRangeInput = mean(X, 1);    
            notBias = find(rangeInput ~= 0);            
            X = bsxfun(@rdivide, bsxfun(@minus, X(:, notBias), ...
                meanRangeInput(notBias)), rangeInput(notBias));      
            hasNoBias = (length(notBias) == length(meanRangeInput)) | ...
                isempty(find(rangeInput == 0 & meanRangeInput ~= 0, 1));
            if(hasNoBias)
                aBiasIndex = [];
            else
                aBiasIndex = find(rangeInput == 0 & meanRangeInput ~= 0, 1);
            end
            X = [ones(size(X, 1), 1) X];
            
            % learn the model
            DhalfX = bsxfun(@times, X, sqrt(weights));
            try
                cholM = chol((DhalfX' * DhalfX) + regul * eye(size(X, 2)));
            catch %cholM not invertible
                warning('matrix not invertible; regularization term too small?');
                param = zeros(length(rangeInput), 1);
                validationScore = inf;
                return;
            end
            tparam = cholM \ (bsxfun(@times, X, weights) / cholM)' * y;
            
            
            % reverse the normalization
            if(~hasNoBias) % has a bias in the data
                bias = 0;
                param = zeros(length(rangeInput), 1);
                param(notBias) = tparam(2:end) ./ rangeInput(notBias)';
                param(aBiasIndex) = (tparam(1) - meanRangeInput(notBias) * param(notBias)) ...
                    / meanRangeInput(aBiasIndex);
            else
                param = zeros(length(rangeInput), 1);
                param(notBias) = tparam(2:end) ./ rangeInput(notBias)';
                bias = tparam(1) - meanRangeInput(notBias) * param(notBias);
            end
            
            % compute score on validation set if asked for
            if(nargout > 2)
                validationScore = sum(vweights .* ((vX * param - vy) .^ 2));
            end
        end
        
        %% learn mean squared error quad model f(x) = 1/2 (x-mu)' A (x-mu) + b
        function [thetas, quadModel, phix] = learnCentredQuadModel(x, y, mu, regularization, weights)
            if(~exist('weights', 'var') || isempty(weights))
                weights = ones(size(X, 1), 1) / size(X, 1);
            end
            if(length(mu) == size(mu, 1))
                mu = mu';
            end
            dim = size(x, 2);            
            phix = static_optimization_algs.Quadratic.getQuadFeatures(bsxfun(@minus, x, mu), 1, 0); % do not add the linear part
            thetas = static_optimization_algs.Quadratic.MSEEstimator(phix, y, weights, regularization);
            quadModel = static_optimization_algs.Quadratic.unpackParams(thetas, dim, [], 1);
            %%% expanding the quadratic expression
            quadModel.w = -quadModel.A*mu';
            quadModel.b = -.5 * mu * quadModel.w + quadModel.b;
        end

        %% learn mean squared error quad model f(x) = 1/2 (x-mu)' A (x-mu) + b, where A is n.d.
        function [thetas, quadModel, phix] = learnCentredConcaveQuadModel(x, y, mu, regularization, weights)
            if(~exist('weights', 'var') || isempty(weights))
                weights = ones(size(X, 1), 1) / size(X, 1);
            end
            if(length(mu) == size(mu, 1))
                mu = mu';
            end

            dim = size(x, 2);
            [thetas, quadModel, phix] = static_optimization_algs.Quadratic.learnCentredQuadModel(x, y, mu, regularization, weights);
            %%% Ensure that A is negative definite (theoretically n.s.d. is enough, but chol is much faster than eig on matlab)
            [isNd, quadModel.A] = static_optimization_algs.Quadratic.isNegativeDefinite(quadModel.A);            
            if(~isNd)
                %%% put back the model in centred form
                quadModel.b = quadModel.b + .5 * mu * quadModel.w;
                quadModel.w = [];
                %%% Improve over thetas by gradient descent
                [thetas, quadModel] = static_optimization_algs.Quadratic.optimizeThetas(quadModel, dim, phix, y, regularization, weights);
                %%% Expand the quad expression again
                quadModel.w = -quadModel.A*mu';
                quadModel.b = -.5 * mu * quadModel.w + quadModel.b;
            end
        end

    end
end
