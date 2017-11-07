function [L, alpha] = getLowerCholesky(D, knownY, opts)
    K = scaleKernel(D, opts);
    alpha = [];
    hyperN = opts.hyperN;
    %p = 1;
    %while p > 0 && hyperN <= 20
        Knoise = K + hyperN.*eye(size(K));
        [L,p] = chol(Knoise, 'lower');
        if p > 0
            L = [];
            %break;
        end
        %hyperN = hyperN * 2;
    %end
    
    if nargout > 1 && ~isempty(L)
        alpha = L'\(L \ knownY);        %Rasmussen, Williams 2005
    end
end