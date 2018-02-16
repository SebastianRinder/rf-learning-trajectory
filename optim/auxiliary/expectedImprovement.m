function [EI, mean, variance] = expectedImprovement(testX, knownX, trajectories, L, alpha, func, bestY)
    [mean, variance] = gaussianProcess(testX, knownX, trajectories, L, alpha, func);
    
    stdY = sqrt(max(0,variance)); %avoid complex numbers
    tau = 0.00;
    v = (mean - bestY - tau) ./ stdY; %Brochu 2010
    v(stdY <= 0) = 0;
    EI = (mean - bestY - tau).*normcdf(v) + stdY.*normpdf(v);
    EI(stdY <= 0) = 0;
    EI(EI<0)=0;
end