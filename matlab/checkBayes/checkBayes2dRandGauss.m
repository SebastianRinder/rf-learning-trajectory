function checkBayes2dRandGauss()
    addpath('..');
    addpath('../matlab');
    load('vararginForFit.mat','varargin');
    
    %[func, lb, ub] = objectiveFcn(isGaussMix);
    load('gaussMix.mat');
    useMatlabBayes = 0;
    isGaussMix = 1;
    if isGaussMix
        if useMatlabBayes
            objFun = @(x) -1.*func.eval(table2array(x));
        else
            objFun = @(x) -1.*func.eval(x);
        end
    else
        objFun = @(x) func(x);
    end
%         
%     points = 100;   
%     x1 = linspace(lb(1),ub(1), points);
%     x2 = linspace(lb(2),ub(2), points);
%     [X1,X2] = meshgrid(x1,x2);
%     
%     obj = zeros(points,points);
%     for j = 1:points
%         for k = 1:points
%             obj(j,k) = objFun([X1(j,k), X2(j,k)]);
%         end
%     end

    dim = 2;
    acqFcn = @expectedImprovement;
    isDeterministic = false;
    bayOptSteps = 100;
    initialPolicies = 4;
    objective = zeros(bayOptSteps,1);
    policy = randPolicy(lb, ub, initialPolicies);
    for i = 1:initialPolicies
        if useMatlabBayes
            objective(i,1) = objFun(array2table(policy(i,:)));
        else
            objective(i,1) = objFun(policy(i,:));
        end
    end     
    
    if ~useMatlabBayes
        for i=1+initialPolicies:bayOptSteps+initialPolicies
            gprMdl = fitGP(policy(1:i-1,:), objective(1:i-1), false);
            
            fMean = @(x) predict(gprMdl, x);
            [~,minFMean] = globalMin(fMean, lb, ub);
            
            negAcqFcn = @(x) -acqFcn(x, gprMdl, minFMean);
            policy(i,:) = globalMin(negAcqFcn, lb, ub);
            objective(i,1) = objFun(policy(i,:));

            disp(['step : ',num2str(i-initialPolicies) ,'  |  cumulative reward: ', num2str(objective(i,1))]);        

            if mod(i-initialPolicies,bayOptSteps) == 0
                if ~exist('myFig', 'var')
                    myFig = figure();
                end
                plotting(policy(1:i,:), acqFcn, gprMdl, minFMean, myFig);
            end
        end
    else
%         global fitVars;
%         fitVars = [];
        for i=1:dim
            X(1,i) = optimizableVariable(['X',int2str(i)], [lb(1,i) ub(1,i)]);
        end
        %fun = @(x) -objFun(x);
        ret = bayesopt(objFun,X,'IsObjectiveDeterministic',isDeterministic,...
        'AcquisitionFunctionName','expected-improvement',...
        'MaxObjectiveEvaluations',100);
    end

end

function [obj, lb, ub, hyper] = objectiveFcn(isGaussMix)
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

