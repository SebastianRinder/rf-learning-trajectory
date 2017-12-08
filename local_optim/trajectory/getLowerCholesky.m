function [L, alpha] = getLowerCholesky(K, knownY, toHyper)    
    alpha = [];    
    hyperSigma = mean(std(knownY,0,2).^2);
    if hyperSigma == 0
        hyperSigma = 1;
    end
    
    p = 1;
    doublings = 0;
    while p > 0 && doublings < 100
        Knoise = K + hyperSigma.*eye(size(K));
        [L,p] = chol(Knoise, 'lower');
        if p > 0
            L = [];
            if ~toHyper
                doublings = doublings + 1;
                hyperSigma = hyperSigma * 2;
            else
                break;
            end
        end        
    end
    
    if nargout > 1 && ~isempty(L)
        alpha = L'\(L \ knownY);        %Rasmussen, Williams 2005
    end
end