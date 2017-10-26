function [L,alpha] = preComputeK(covarianceFcn, knownX, knownY)

    K = covarianceFcn(knownX, knownX);
    Knoise = K + 1e-6*eye(size(K));
    L = chol(Knoise, 'lower');
    alpha = L'\(L \ knownY);        %Rasmussen, Williams 2005
end