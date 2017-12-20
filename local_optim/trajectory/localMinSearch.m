function [xMin, yMin] = localMinSearch(fcn, dist, beta)
    
    options = optimoptions('fmincon','Display','off');
    
%     try        
%         warning off MATLAB:singularMatrix
%         warning off MATLAB:nearlySingularMatrix
        [xMin, yMin] = fmincon(fcn, dist.mu, [],[],[],[],[],[], @nonlconFcn, options);
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

