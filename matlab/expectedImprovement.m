function [EI, mean, variance] = expectedImprovement(testX, knownX, knownY, opts)
    if ~isempty(opts.L)
        [mean, variance] = gaussianProcessModel(testX, knownX, knownY, opts);

        stdY = sqrt(max(0,variance));    
        v = (mean - opts.bestY) ./ stdY; %Brochu 2010
        EI = stdY .* (v .* normcdf(v) + normpdf(v));

        EI(isnan(EI)) = 0;
    else
        EI(1:size(testX,1),1) = 0;
    end
end

% [FMean, YSD, ~] = predict(ObjectiveFcnGP, X);
% FSD = sqrt(max(0, YSD.^2 - ObjectiveFcnGP.Sigma.^2));
% GammaX = (FBest - FMean)./FSD;
% PI = normcdf(GammaX, 0, 1);