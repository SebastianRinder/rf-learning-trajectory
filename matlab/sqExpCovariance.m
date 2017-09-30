function covMatrix = sqExpCovariance(Xi, Xj, theta, ~, ~)
    hyper.l = theta(1);
    hyper.f = theta(2);
    if isempty(Xi)
        Xi = Xj;
    end
    covMatrix = zeros(size(Xi,1), size(Xj,1));
    for i=1:size(Xi,1)
        for j=1:size(Xj,1)
            r = Xi(i,:) - Xj(j,:);
            covMatrix(i,j) = hyper.f ^ 2 * exp( -(r*r') / (2 * hyper.l ^ 2));
        end
    end
end