function [EI, mean, variance] = expectedImprovement(testX, knownX, trajectories, L, alpha, func, bestY)
    [mean, variance] = gaussianProcess(testX, knownX, trajectories, L, alpha, func);

    stdY = sqrt(max(0,variance)); %avoid complex numbers
    v = (mean - bestY - 0.01) ./ stdY; %Brochu 2010
    v(stdY <= 0) = 0;
    EI = (mean - bestY - 0.01).*normcdf(v) + stdY.*normpdf(v);
    EI(stdY <= 0) = 0;   
end