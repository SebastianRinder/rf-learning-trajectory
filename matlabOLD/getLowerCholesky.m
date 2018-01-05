function [L, alpha] = getLowerCholesky(D, knownY, opts, toHyper)
    K = scaleKernel(D, opts);
    alpha = [];
    hyperN = opts.hyperN;
    p = 1;
    doublings = 0;
    while p > 0 && doublings < 10
        Knoise = K + hyperN.*eye(size(K));
        [L,p] = chol(Knoise, 'lower');
        if p > 0
            L = [];
            if ~toHyper
                doublings = doublings + 1;
                hyperN = hyperN * 2;
            else
                break;
            end
        end        
    end
    
    
    
    if nargout > 1 && ~isempty(L)
        alpha = L'\(L \ knownY);        %Rasmussen, Williams 2005
    end
end