function [L, alpha, hyperN] = preComputeK(knownX, knownY, opts)
%     noise = std(knownY)/sqrt(2);
    K = opts.covarianceFcn(knownX, knownX, opts);
    
    hyperN = opts.hyper(3);
    p = 1;
    while p > 0 && hyperN <= 2e1
        Knoise = K + hyperN.*eye(size(K));
        [L,p] = chol(Knoise, 'lower');
        if p > 0 && nargout == 1
            L = [];
            break;
        end
        hyperN = hyperN * 2;
    end
    
    if nargout > 1
        try
            alpha = L'\(L \ knownY);        %Rasmussen, Williams 2005
        catch me
            disp('fail');
        end
    end
end