classdef GradientCheck
    
    methods (Static)
        function [diff, numGrad] = do(func, x, d)
            numGrad = zeros(size(x));
            [~, dx] = func(x);
            xk = x;
            for k = 1:length(x)
                xk(k) = x(k) + d;
                ykp = func(xk);
                xk(k) = x(k) - d;
                ykn = func(xk);
                xk(k) = x(k);
                numGrad(k) = (ykp - ykn) / (2*d);
            end
            diff = dx - numGrad;
        end
    end
    
end

