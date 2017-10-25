function boundedVal = applyBound(x, bounds)
    minBound = min(bounds);
    maxBound = max(bounds);
    boundedVal = min(max(x, minBound), maxBound);    
end