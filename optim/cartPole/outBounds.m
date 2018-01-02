function isOutBounds = outBounds(x, bounds)
    minBound = min(bounds);
    maxBound = max(bounds);
    isOutBounds = x < minBound || x > maxBound;
end