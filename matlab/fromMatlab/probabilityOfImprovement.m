function [PI, FSD, GammaX, FMean] = probabilityOfImprovement(X, ObjectiveFcnGP, FBest)
    [FMean, YSD, ~] = predict(ObjectiveFcnGP, X);
    FSD = sqrt(max(0, YSD.^2 - ObjectiveFcnGP.Sigma.^2));
    GammaX = (FBest - FMean)./FSD;
    PI = normcdf(GammaX, 0, 1);
end
