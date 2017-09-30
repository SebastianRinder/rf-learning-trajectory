function ret = mainCartPole()
    addpath('..');
    addpath('../matlab');
    
    global state0;
    global executionTimeSteps;
    global worldBounds;
    
    bayOptSteps = 100;
    executionTimeSteps = 1000;
    trajectoriesPerPolicy = 1;
    dim = 4;
    initialPolicies = 40;
    trials = 1;

	worldBounds.position = [-10, 10];
    worldBounds.angle = [-pi/2, pi/2];
    worldBounds.rewardPosition = [-0.5, 0.5];
    worldBounds.rewardAngle = [-pi/9, pi/9];
    worldBounds.action = [-1,1];
    
    state0.position = 0;
    state0.velocity = 0;
    state0.acceleration = 0;
    state0.angle = 0;
    state0.angleVelocity = 0;
    
    ub = 1.*ones(1,dim);
    lb = -1.*ub;

    isDeterministic = 0;
    
    useMatlabBayes = 1;
    useMatlabGP = 1;
    
    optimizeNoiseHyper = 0;
    acqFcn = @expectedImprovement;
    
    %covarianceFcn = @sqExpCovariance;
    covarianceFcn = @trajectoryCovariance;
    actionSelectionFcn = @actionSelectionCartPole;
    
    if ~useMatlabBayes
        ret = cell(trials,1);
        for trial = 1:trials            
            policy = zeros(bayOptSteps + initialPolicies, dim);
            negObjective = zeros(bayOptSteps + initialPolicies, 1);
            tempObjective = zeros(1, trajectoriesPerPolicy);
            trajectory = cell(bayOptSteps + initialPolicies, trajectoriesPerPolicy);
            %hypers = zeros(bayOptSteps + initialPolicies, 6);

            for i=1:initialPolicies
                policy(i,:) = randPolicy(lb, ub, 1);
                for j = 1:trajectoriesPerPolicy
                    [tempObjective(1,j), trajectory{i,j}] = execPolicyCartPole(policy(i,:), state0, executionTimeSteps, worldBounds);
                end
                negObjective(i,1) = -mean(tempObjective);
            end
           
            for i = initialPolicies+1:initialPolicies+bayOptSteps
                xTrain = policy(1:i-1,:);
                yTrain = negObjective(1:i-1,:);
                hyper = setHypers(yTrain, lb, ub, isDeterministic);
                traj = trajectory(1:i-1,:);
                
%                 K = covarianceFcn(xTrain, xTrain, hyper, traj, actionSelectionFcn);
%                 histogram(K);
                
                if useMatlabGP
                    kernelFcn = @(Xm, Xn, theta) covarianceFcn(Xm, Xn, theta, traj, actionSelectionFcn);
                    gprMdl = fitGP(xTrain, yTrain, hyper, kernelFcn, optimizeNoiseHyper);
                    meanMdlFcn = @(X) predict(gprMdl, X);
                    [~,minMeanMdl] = globalMinSearch(meanMdlFcn, lb, ub, useMatlabGP);

                    negAcqFcn = @(X) -acqFcn(X, gprMdl, minMeanMdl);
                    policy(i,:) = globalMinSearch(negAcqFcn, lb, ub, useMatlabGP);
                else
                    K = covarianceFcn([], xTrain, hyper, traj);
                    meanMdlFcn = @(X) gaussianProcess(X, xTrain, yTrain, covarianceFcn, K, hyper, traj);
                    [~,minMeanMdl] = globalMin(meanMdlFcn, lb, ub, useMatlabGP);
                    
                    negAcqFcn = @(X) -acquisition(X, xTrain, yTrain, minMeanMdl, covarianceFcn, K, hyper, traj);
                    policy(i,:) = globalMin(negAcqFcn, lb, ub, useMatlabGP);
                end
                
                for j = 1:trajectoriesPerPolicy
                    [tempObjective(1,j), trajectory{i,j}] = execPolicyCartPole(policy(i,:), state0, executionTimeSteps, worldBounds);
                end
                negObjective(i,1) = -mean(tempObjective);

                disp(['trial : ',num2str(trial) ,...
                    ' | step : ',num2str(i-initialPolicies) ,...
                    ' | negCR: ',num2str(negObjective(i,1))...
                    ]);
                %visCartPole(traj{i,m},worldBounds);
            end
            
            ret{trial,1}.policy = policy;
            ret{trial,1}.objective = negObjective;

            save('ret.mat','ret');
        end
    else
        for i=1:dim
            policy(1,i) = optimizableVariable(['policy',int2str(i)], [lb(1,i) ub(1,i)]);
        end
        %fun = @(x) -objectiveFun(x, state0, executionTimeSteps, worldBounds);
        ret = bayesopt(@negObjectiveFcn,policy,actionSelectionFcn,covarianceFcn,...
            'IsObjectiveDeterministic',isDeterministic,...
            'AcquisitionFunctionName','expected-improvement',...
            'MaxObjectiveEvaluations',bayOptSteps);
    end
end

function [objective, constraints, trajectory] = negObjectiveFcn(x) %, state0, executionTimeSteps, worldBounds)
    global state0;
    global executionTimeSteps;
    global worldBounds;
    dim = size(x,2);
    
    nextTheta = zeros(1,dim);
    for i=1:dim
        nextTheta(1,i) = x.(['policy',int2str(i)]);
    end
    [objective, trajectory] = execPolicyCartPole(nextTheta, state0, executionTimeSteps, worldBounds);
    objective = -objective;
    constraints = [];
end