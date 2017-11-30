function [L, alpha] = getLowerCholesky(K, knownY, opts, toHyper)    
    alpha = [];
    hyperSigma = opts.hyperSigma;
    p = 1;
    doublings = 0;
    while p > 0 && doublings < 10
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