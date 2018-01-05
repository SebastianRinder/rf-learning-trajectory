function [xMin,yMin] = globalMin(fcn,lb,ub,toHyper,debugPlotting)
    if ~toHyper
        xRand = randPolicy(lb,ub,10000);
        yRand = fcn(xRand);
        if debugPlotting
            selectFigure('Expected Improvement values (sorted)');
            plot(sort(-yRand));
            pause(0.1);
        end
    else
        xRand = randPolicy(lb,ub,10000);
        yRand = zeros(10000,1);
        %yRand = fcn(xRand);
        for i=1:10000
            yRand(i,1) = fcn(xRand(i,:));
        end
        if debugPlotting
            selectFigure('Hyper Optimization (sorted)');
            plot(sort(-yRand));
            pause(0.1);
        end
    end
    
    [~, rows] = sort(yRand, 'ascend');
    x0 = xRand(rows(1:10), :);
    %y0 = yRand(rows(1:10), :);
    
    opts = optimset('MaxIter', 10, 'TolX', Inf, 'TolFun', 0.001, 'Display', 'off');

    for row = 10:-1:1
        [x(row,:), y(row)] = fminsearch(@boundedFcn, x0(row,:) ,opts);        
    end
    
    [yMin, i] = min(y);
    xMin = x(i,:);
    
    function y = boundedFcn(x)
        if any(x < lb | x > ub)
            y = Inf;
        else
            y = fcn(x);
        end
    end
end

