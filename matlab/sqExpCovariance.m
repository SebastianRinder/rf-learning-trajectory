function covMatrix = sqExpCovariance(Xi, Xj, theta, ~)
    hyper.l = theta(1);
    hyper.f = theta(2);
    if isempty(Xi)
        Xi = Xj;
    end
    for i=size(Xi,1):-1:1
        for j=size(Xj,1):-1:1
            r = Xi(i,:) - Xj(j,:);
            d(i,j) = r*r';            
        end
    end
    covMatrix = hyper.f ^ 2 .* exp( d ./ -(2 * hyper.l ^ 2));
end