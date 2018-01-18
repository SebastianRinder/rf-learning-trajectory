function plotFcn(fcn,lb,ub,nrPoints,dim)
% global nv;
    if dim == 1
        X = linspace(lb,ub);
        
        fVals = zeros(nrPoints,1);
        for j = 1:nrPoints
            fVals(j) = fcn(X(j));
        end
        plot(X, fVals);
    elseif dim == 2        
%         nvAll = zeros(nrPoints,nrPoints);
                
        x1 = linspace(lb(1),ub(1), nrPoints);
        x2 = linspace(lb(2),ub(2), nrPoints);
        [X1,X2] = meshgrid(x1,x2);

        fVals = zeros(nrPoints,nrPoints);
        for j = 1:nrPoints
            for k = 1:nrPoints
                fVals(j,k) = fcn([X1(j,k), X2(j,k)]);
%                 nvAll(j,k) = nv;
            end
        end

        surf(X1,X2,fVals);
        alpha 0.5
        
%         figure;
%         nvAll(nvAll>0.1) = 0.1;
%         surf(X1,X2,nvAll);
    end
end