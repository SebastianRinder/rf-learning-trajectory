function [EI, mean, variance] = expectedImprovement(testX, knownX, knownY, opts)
    [mean, variance] = gaussianProcessModel(testX, knownX, knownY, opts);
    if any(variance <= 0)
        if size(variance,1) > 1
            keyboard;
        end
        EI = 0;
    else        
        stdY = sqrt(variance);    
        v = (mean - opts.bestY) ./ stdY;
        EI = stdY .* (v .* normcdf(v) + normpdf(v));
    end
end