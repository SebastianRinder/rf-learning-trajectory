function [xMin, yMin] = localMinSearch(fcn, dist, beta)
    c = 1;
%     while c
%         x0 = dist.getSamples(1);
%         c = nonlconFcn(x0);
%     end
    
    x0 = dist.mu;
    
    %options = optimoptions('fmincon','Display','iter');
    
%     try        
%         warning off MATLAB:singularMatrix
%         warning off MATLAB:nearlySingularMatrix
        [xMin, yMin] = fmincon(fcn, x0, [],[],[],[],[],[], @nonlconFcn); %, options);
%         warning on MATLAB:singularMatrix
%         warning on MATLAB:nearlySingularMatrix
%     catch me
%         xMin = [];
%         yMin = [];
%     end
    
    function [c, ceq] = nonlconFcn(x)
        probaThresh = -beta;
        mahDistMu = pdist2(x, dist.mu, 'mahalanobis', dist.getCov) .^ 2;
        cutOffDist = chi2inv(probaThresh, length(dist.mu));
        
%         c = double(mahDistMu >= cutOffDist);
        c = mahDistMu - cutOffDist;
        ceq = [];
    end
end

