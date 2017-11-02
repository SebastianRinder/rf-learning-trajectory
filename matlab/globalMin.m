function [xMin,yMin] = globalMin(fcn,lb,ub,toHyper)
    if ~toHyper
        xRand = randPolicy(lb,ub,10000);
        yRand = fcn(xRand);
        plot(sort(yRand));
        pause(0.1);
    else
        xRand = 10.^(randPolicy(log10(lb),log10(ub),10000));
        for i=10000:-1:1
            yRand(i,1) = fcn(xRand(i,:));
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

