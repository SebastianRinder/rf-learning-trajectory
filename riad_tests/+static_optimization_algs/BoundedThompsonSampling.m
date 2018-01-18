classdef BoundedThompsonSampling
    properties
        gp;
        x;
        Cy;
        thresh_dist;
        sigma;
    end
    
    methods
        function obj = BoundedThompsonSampling(gp, x, y, proba_thresh, sigma)
            obj.gp = gp;
            obj.x = x;
            [~, C] = gp_trcov(gp, x);
            obj.Cy = C \ y;
            obj.thresh_dist = chi2inv(proba_thresh, size(x, 2));
            obj.sigma = sigma + 1e-10;
        end
        
        function val = eval(obj, xtest)
            if(xtest * xtest' / obj.sigma > obj.thresh_dist)
                val = nan;
                return;
            end
            kx = gp_cov(obj.gp, xtest, obj.x)';
            val = kx' * obj.Cy;
        end
    end
    
end

