function [EI, meanY, stdvY] = acquisition(postX, priorY, priorX, covarianceFun, K, hyper, traj)
    [meanY, varY] = gaussianProcess(postX, priorY, priorX, covarianceFun, K, hyper, traj);
    stdvY = sqrt(max(0, varY));
   
    if stdvY > 0    
        maxY = max(priorY) + 1e-6;
        v = (meanY - maxY) / stdvY;
        PI = normcdf(v);                    %PI: probability of improvement
        EI = stdvY *(v * PI + normpdf(v));   %EI: expected improvement
    else
        EI = 0;
    end
    
%     minY = min(priorY);
%     gammaX = (minY - meanY) / stdvY;
%     PI = normcdf(gammaX);
%     EI = stdvY * (gammaX * PI + normpdf(gammaX));
    

%     minY = min(priorY);
%     PI = normcdf(minY,meanY,stdvY);
%     EI = (minY - meanY) * PI + stdvY * normpdf(minY, meanY, stdvY);
end