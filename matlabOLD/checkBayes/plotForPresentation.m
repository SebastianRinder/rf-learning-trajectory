% nbGauss = 8;
% muRange = 10;
% minSigma = 0.01;
% func = RandomGaussianMixtureFunction(nbGauss, muRange, minSigma);
load('func');
[lb, ub] = func.getRange();

points = 100;   
x1 = linspace(lb(1),ub(1), points);
x2 = linspace(lb(2),ub(2), points);
[X1,X2] = meshgrid(x1,x2);
% 
obj = zeros(points,points);
for j = 1:points
    for k = 1:points
        obj(j,k) = func.eval([X1(j,k), X2(j,k)]);
    end
    if mod(j,100) == 0
        disp([num2str(j*100 / points), ' %']);
    end
end

% load('obj.mat');

surf(X1,X2,obj,'EdgeColor','none');