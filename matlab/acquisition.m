function [EI, meanMdl, stdMdl] = acquisition(X, trainX, trainY, minMeanMdl, covarianceFcn, K, hyper, traj)
    [meanMdl, stdMdl] = gaussianProcess(X, trainX, trainY, covarianceFcn, K, hyper, traj);
    
    stdv = sqrt(max(0, stdMdl.^2 - hyper.noise.^2));
    gamma = (minMeanMdl - meanMdl)./stdv;
    PI = normcdf(gamma, 0, 1);
    EI = stdv.*(gamma.*PI + normpdf(gamma, 0, 1));


%     [FMean, YSD, ~] = predict(ObjectiveFcnGP, X);
%     FSD = sqrt(max(0, YSD.^2 - ObjectiveFcnGP.Sigma.^2));         	% Want SD of F, not Y. Need max to avoid complex sqrt.
%     GammaX = (FBest - FMean)./FSD;
%     PI = normcdf(GammaX, 0, 1);
%     EI = FSD.*(GammaX.*PI + normpdf(GammaX, 0, 1));

%     [meanY, varY] = gaussianProcess(postX, priorY, priorX, covarianceFun, K, hyper, traj);
%     stdvY = sqrt(varY);
    
%     [meanF, stdvY] = predict(gprMdl, postX);

    %stdvF = sqrt(max(0,stdvY.^2 - gprMdl.Sigma.^2));
%     stdvF = max(0,stdvY);
%     
%     bestF = min(priorY); % + 1e-6;
    
%     gamma = (bestF - meanF) ./ stdvF;
    
%     bestF = max(priorY); % + 1e-6;
%     gamma = (meanF - bestF) ./ stdvF;

%     PI = normcdf(gamma);                    %PI: probability of improvement
%     EI = stdvF .* (gamma .* PI + normpdf(gamma));   %EI: expected improvement

    %EI2 = expectedImprovement(postX, gprMdl, bestF);
    
%     minY = min(priorY);
%     gammaX = (minY - meanY) / stdvY;
%     PI = normcdf(gammaX);
%     EI = stdvY * (gammaX * PI + normpdf(gammaX));
    

%     minY = min(priorY);
%     PI = normcdf(minY,meanY,stdvY);
%     EI = (minY - meanY) * PI + stdvY * normpdf(minY, meanY, stdvY);
end