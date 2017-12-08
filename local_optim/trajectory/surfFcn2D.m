function surfFcn2D(fcn,lb,ub,nrPoints)
    x1 = linspace(lb(1),ub(1), nrPoints);
    x2 = linspace(lb(2),ub(2), nrPoints);
    [X1,X2] = meshgrid(x1,x2);

    fVals = zeros(nrPoints,nrPoints);
    for j = 1:nrPoints
        for k = 1:nrPoints
            fVals(j,k) = fcn([X1(j,k), X2(j,k)]);
        end
    end

    surf(X1,X2,fVals);
    pause(0.1);
end