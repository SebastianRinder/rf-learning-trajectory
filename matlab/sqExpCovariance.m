function Knm = sqExpCovariance(Xm, Xn, theta, ~)
    hyper.l = max(theta(1), 1e-6);
    hyper.f = max(theta(2), 1e-6);
    if isempty(Xm)
        Xm = Xn;
    end
    for i=size(Xm,1):-1:1
        for j=size(Xn,1):-1:1
            r = Xm(i,:) - Xn(j,:);
            d(i,j) = r*r';            
        end
    end
    Knm = d ./ (hyper.l^2);
    Knm = (hyper.f^2) * exp(-0.5*Knm);
end