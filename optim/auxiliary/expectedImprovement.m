function [EI, mean, variance] = expectedImprovement(testX, knownX, trajectories, L, alpha, func, bestY)
    if ~isempty(L)
        [mean, variance] = gaussianProcess(testX, knownX, trajectories, L, alpha, func);
        
        stdY = sqrt(max(0,variance)); %avoid complex numbers
        v = (mean - bestY - 0.001) ./ stdY; %Brochu 2010
        v(stdY <= 0) = 0;
        %EI(stdy > 0) = stdY .* (v .* normcdf(v) + normpdf(v));
        EI = (mean - bestY - 0.001).*normcdf(v) + stdY.*normpdf(v);
        EI(stdY <= 0) = 0;

        %EI(isnan(EI)) = -Inf;
    else
        disp('ei fail');
        EI(1:size(testX,1),1) = 0;
    end
    if ~isa(EI,'double')
        keyboard;
    end
end

% [FMean, YSD, ~] = predict(ObjectiveFcnGP, X);
% FSD = sqrt(max(0, YSD.^2 - ObjectiveFcnGP.Sigma.^2));
% GammaX = (FBest - FMean)./FSD;
% PI = normcdf(GammaX, 0, 1);