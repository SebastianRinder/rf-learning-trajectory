classdef NonParamImportanceSampling
    methods (Static)
        function [y, dy] = nbEffSamples(distMat, initLogProbas, rho)
            importanceProbas = sum(exp(-rho * distMat), 2);
            wis = exp(initLogProbas) ./ importanceProbas;
            wis2 = wis .^ 2;
            what = sum(wis);
            w2hat = sum(wis2);
            y = what / w2hat * what;
            dImportanceProbas = -sum(exp(-rho * distMat) .* distMat, 2);
            dy = 2 * y * sum( ((wis2 / w2hat) - (wis / what)) .* dImportanceProbas ./ importanceProbas);
        end
        
        function [y, dy] = minusFun(f, param)
            [y, dy] = f(param);
            y = -y;
            dy = -dy;
        end

        
        function [rho] = optimizeRho(distMat, initLogProbas, rho0)
            fun = @(rho) static_optimization_algs.NonParamImportanceSampling.nbEffSamples(distMat, initLogProbas, rho);
            minusFun = @(param) static_optimization_algs.NonParamImportanceSampling.minusFun(fun, param);
            infRho = 1e-16;
            options =  optimset('GradObj','on', 'DerivativeCheck', 'off', 'Display', 'off');
            rho = fmincon(minusFun, rho0, [], [], [], [], infRho, inf,[],options);
        end
        
        function [d] = mahDist(x, y, D) %x is 1 x n vector and y is m x n matrix
            diff = bsxfun(@minus, y, x);
            d = dot(diff * D, diff, 2);
        end

        
        function [logProbas] = getLogProbas(samples, initPrecision, initLogProbas)
            distMat = squareform(pdist(samples, @(x, y) static_optimization_algs.NonParamImportanceSampling.mahDist(x, y, initPrecision)));
            rho = static_optimization_algs.NonParamImportanceSampling.optimizeRho(distMat, initLogProbas, 1);
            probas = sum(exp(-rho * distMat), 2);
            logProbas = log(probas / sum(probas));
        end
    end
    
end

