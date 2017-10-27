function checkBayes2dRandGauss()
    addpath('..');
    
    useMaxMean = 0;
    %[func, lb, ub] = objFcn();
    load('gaussMix.mat');
        
    opts.covarianceFcn = @sqExpCovariance;
    objectiveFcn = @(x) func.eval(x)+randn*1e-6;
       
%     points = 100;   
%     x1 = linspace(lb(1),ub(1), points);
%     x2 = linspace(lb(2),ub(2), points);
%     [X1,X2] = meshgrid(x1,x2);
%     
%     obj = zeros(points,points);
%     for j = 1:points
%         for k = 1:points
%             obj(j,k) = objectiveFcn([X1(j,k), X2(j,k)]);
%         end
%     end

    dim = 2;
    bayOptSteps = 100;
    initialPolicies = 5;
    hypers = logspace(-2,2,5);
    [h1,h2] = meshgrid(hypers);
    hypers = [h1(:), h2(:)];
    
     
    for i=1:initialPolicies
        opts.trajectory.policy(i,:) = randPolicy(lb, ub, 1);
        objective(i,1) = objectiveFcn(opts.trajectory.policy(i,:));
    end

    for i = initialPolicies+1:initialPolicies+bayOptSteps
        knownX = opts.trajectory.policy;
        knownY = objective;

%         opts.hyper.f = mean(std(knownX));
%         opts.hyper.l = std(knownY)/sqrt(2);
        
        opts.hyper.f = mean(std(knownX))/2000;
        opts.hyper.l = std(knownY)*10;
        
        [opts.L, opts.alpha] = preComputeK(knownX, knownY, opts);

        if useMaxMean
            negGPModel = @(testX) -gaussianProcessModel(testX, knownX, knownY, opts);
            [~,negMaxMean] = globalMin(negGPModel, lb, ub);
            opts.bestY = -negMaxMean;                
        else
            opts.bestY = max(knownY);
        end
        
        negAcqFcn = @(testX) -expectedImprovement(testX, knownX, knownY, opts);
        opts.trajectory.policy(i,:) = globalMin(negAcqFcn, lb, ub);
        
        objective(i,1) = objectiveFcn(opts.trajectory.policy(i,:));

        disp(['trial : ',num2str(0) ,...
            ' | step : ',num2str(i-initialPolicies) ,...
            ' | cr: ',num2str(objective(i,1)),...
            ' | f: ',num2str(opts.hyper.f),...
            ' | l: ',num2str(opts.hyper.l)...
            ]);
        
        knownX = opts.trajectory.policy;
        knownY = objective;

        if mod(i,10) == 0
            plotting(knownX, knownY, opts);
        end
        
    end

end

function [func, lb, ub] = objFcn()
    nbGauss = 25;
    muRange = 10;
    minSigma = 2;
    func = RandomGaussianMixtureFunction(nbGauss, muRange, minSigma);
    [lb, ub] = func.getRange();
end

