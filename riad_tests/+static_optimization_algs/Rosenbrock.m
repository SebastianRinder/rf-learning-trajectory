classdef Rosenbrock
    methods
        function y = f(obj, x)
            y = (ones(size(x, 1), 1) - x(:, 1)).^2 + 100 * (x(:, 2) - x(:, 1).^2).^2;
        end
    end
end