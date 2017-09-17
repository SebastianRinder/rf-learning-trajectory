function checkBayes2dRandGauss()
    addpath('..');
    
    isGaussMix = 1;
    [func, lb, ub, hyper] = objective(isGaussMix);    
    if isGaussMix
        objFun = @(x) func.eval(x);
    else
        objFun = @(x) func(x);
    end
            
	theta(1,:) = randTheta(lb, ub);
    eta(1,1) = objFun(theta(1,:));
    disp(['step : ',num2str(1) ,'  |  cumulative reward: ', num2str(eta(1,1))]);    

    points = 100;
    x1 = linspace(lb(1),ub(1), points);
    x2 = linspace(lb(2),ub(2), points);
    [X1,X2] = meshgrid(x1,x2);
       
    for j = 1:points
        for k = 1:points
            obj(j,k) = objFun([X1(j,k), X2(j,k)]);
        end
    end

    close all;
    figure;
    for i=2:100
        clf;
        K = sqExpCovariance([], theta(1:i-1,:), [], hyper);        
        toMinimize = @(x) -acquisition(x, eta(1:i-1,:), theta(1:i-1,:), @sqExpCovariance, K, hyper, []);
        theta(i,:) = globalMin(toMinimize, [], lb, ub);
        eta(i,1) = objFun(theta(i,:));

        disp(['step : ',num2str(i) ,'  |  cumulative reward: ', num2str(eta(i,1))]);        

        if mod(i,1) == 0
            for j = 1:100
                for k = 1:100
                    [EI(j,k), meanY(j,k), var(j,k)] = acquisition([X1(j,k), X2(j,k)], eta(1:i-1,1),...
                        theta(1:i-1,:), @sqExpCovariance, K, hyper, []);
                end
            end

            subplot(2,2,1);
            contour(X1,X2,obj);
            colorbar;

            subplot(2,2,2);
            contour(X1,X2,meanY);
            colorbar;

            subplot(2,2,3);
            contour(X1,X2,EI);
            hold on;
            plot(theta(:,1),theta(:,2),'r+');
            plot(theta(end,1),theta(end,2),'k*');
            colorbar;

            subplot(2,2,4);
            surf(X1,X2,EI);
            colorbar;
        end
    end

end

function [obj, lb, ub, hyper] = objective(isGaussMix)
    if isGaussMix
        nbGauss = 25;
        muRange = 10;
        minSigma = 2;
        obj = RandomGaussianMixtureFunction(nbGauss, muRange, minSigma);
        [lb, ub] = obj.getRange();
        hyper = [0.1, 1];
    else
        obj = @(x) sin(x(1))+cos(x(2));
        lb = [-2*pi,0];
        ub = [2*pi, 4*pi];
        hyper = [1 1];
    end
end

