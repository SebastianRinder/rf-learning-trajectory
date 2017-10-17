function [EI, FMean, FSD] = expectedImprovement(X, ObjectiveFcnGP, FBest)
    [PI, FSD, GammaX, FMean] = probabilityOfImprovement(X, ObjectiveFcnGP, FBest);
    EI = FSD.*(GammaX.*PI + normpdf(GammaX, 0, 1));
end
