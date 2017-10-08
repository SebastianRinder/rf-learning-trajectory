function ret = mainCartPole()
    addpath('..');
    addpath('../fromMatlab');
    
    dim = 4;
    trials = 1;
    bayOptSteps = 200;
    ub = 1.*ones(1,dim);
    lb = -1.*ub;
    isDeterministic = 0;    
    trajectoriesPerPolicy = 1;
    initialPolicies = 4;
    
    useMatlabBayes = 1;
    useMatlabGP = 1;    
    optimizeNoiseHyper = 0;
    acqFcn = @expectedImprovement;
    
	simOpts.bounds.position = [-2.4, 2.4];
    simOpts.bounds.angle = [-12 * pi / 180, 12 * pi / 180];
    %simOpts.bounds.rewardPosition = [-0.5, 0.5];
    %simOpts.bounds.rewardAngle = [-pi/9, pi/9];    
    simOpts.bounds.action = [-1,1]; %force applied to the cart
    %simOpts.bounds.acceleration = [-1,1];
    %simOpts.bounds.velocity = [-1,1];
    
    simOpts.state0.position = 0;
    simOpts.state0.velocity = 0;
    simOpts.state0.acceleration = 0;
    simOpts.state0.angle = 0.02;
    simOpts.state0.angleVelocity = 0;
    
    simOpts.timeSteps = 1000;
    simOpts.execPolicyFcn = @execPolicyCartPole;
    
    global customKernel;
    customKernel.covarianceFcn = @sqExpCovariance;
    %customKernel.covarianceFcn = @trajectoryCovariance;
    customKernel.actionSelectionFcn = @actionSelectionCartPole;
    customKernel.trajectoriesPerPolicy = trajectoriesPerPolicy;
    customKernel.simOpts = simOpts;

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
                    [tempObjective(1,j), trajectory{i,j}] = execPolicyCartPole(policy(i,:), state0, executionTimeSteps, bounds);
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
                    [~,minMeanMdl] = globalMin(meanMdlFcn, lb, ub, useMatlabGP);

                    negAcqFcn = @(X) -acqFcn(X, gprMdl, minMeanMdl);
                    policy(i,:) = globalMin(negAcqFcn, lb, ub, useMatlabGP);
                else
                    K = covarianceFcn([], xTrain, hyper, traj);
                    meanMdlFcn = @(X) gaussianProcess(X, xTrain, yTrain, covarianceFcn, K, hyper, traj);
                    [~,minMeanMdl] = globalMin(meanMdlFcn, lb, ub, useMatlabGP);
                    
                    negAcqFcn = @(X) -acquisition(X, xTrain, yTrain, minMeanMdl, covarianceFcn, K, hyper, traj);
                    policy(i,:) = globalMin(negAcqFcn, lb, ub, useMatlabGP);
                end
                
                for j = 1:trajectoriesPerPolicy
                    [tempObjective(1,j), trajectory{i,j}] = execPolicyCartPole(policy(i,:), state0, executionTimeSteps, bounds);
                end
                negObjective(i,1) = -mean(tempObjective);

                disp(['trial : ',num2str(trial) ,...
                    ' | step : ',num2str(i-initialPolicies) ,...
                    ' | negCR: ',num2str(negObjective(i,1))...
                    ]);
                %visCartPole(traj{i,m},worldBounds);
                
                ret{trial,1}.policy = policy;
                ret{trial,1}.objective = negObjective;

                save('ret.mat','ret');
            end
        end
    else
        for i=1:dim
            policy(1,i) = optimizableVariable(['policy',int2str(i)], [lb(1,i) ub(1,i)]);
        end
        close all;
        tic;

        ret = customBayesopt(@negObjectiveFcn, policy, customKernel,...
            'IsObjectiveDeterministic',isDeterministic,...
            'AcquisitionFunctionName','expected-improvement',...
            'MaxObjectiveEvaluations',bayOptSteps,...
            'NumSeedPoints',initialPolicies);
        toc;
    end
end

function [negObjective, constraints, trajectory] = negObjectiveFcn(x)
    global customKernel;
    dim = size(x,2);
    
    policy = zeros(1,dim);
    for i=1:dim
        policy(1,i) = x.(['policy',int2str(i)]);
    end
    for i = customKernel.trajectoriesPerPolicy:-1:1
        [tempObjective(1,i), trajectory.data{1,i}] = execPolicyCartPole(policy, customKernel.simOpts);
    end
    
    negObjective = -mean(tempObjective);
    constraints = [];
    trajectory.policy = policy;
end