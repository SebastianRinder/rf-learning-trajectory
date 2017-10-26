function EI = expectedImprovement(testX, knownX, knownY, maxObjective, L)
    [mean, variance] = gaussianProcessModel(testX, knownX, knownY, L);
%     if any(variance <= 0)
%         keyboard;
%     end
    stdY = sqrt(variance);
    
    v = (mean - maxObjective) ./ stdY;
    EI = stdY .* (v .* normcdf(v) + normpdf(v));
end