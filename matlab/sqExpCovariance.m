function covMatrix = sqExpCovariance(thetai, thetaj, ~, hyper)
    if isempty(thetai)
        thetai = thetaj;
    end
    covMatrix = zeros(size(thetai,1), size(thetaj,1));
    for i=1:size(thetai,1)
        for j=1:size(thetaj,1)
            r = thetai(i,:) - thetaj(j,:);
            covMatrix(i,j) = hyper(1) * exp( -hyper(2) * (r*r') / 2);
        end
    end
end