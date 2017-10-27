function [L,alpha] = preComputeK(knownX, knownY, opts)

%     noise = std(knownY)/sqrt(2);
    noise = 1e-6;

    K = opts.covarianceFcn(knownX, knownX, opts);
    Knoise = K + noise.*eye(size(K));
    L = chol(Knoise, 'lower');
    alpha = L'\(L \ knownY);        %Rasmussen, Williams 2005
end