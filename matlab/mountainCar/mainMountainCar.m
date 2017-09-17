function ret = mainMountainCar()
    addpath('..');
    
    global execPolicyFun;
    global actionSelectionFun;    
    execPolicyFun = @execPolicyMountainCar;
    actionSelectionFun = @actionSelectionMountainCar;
    
    global state0;
    global executionTimeSteps;
    global worldBounds;
    global trajectoriesPerPolicy;

    bayOptSteps = 200;
    dim = 18;
    initialPolicies = 30;
    trials = 40;
    executionTimeSteps = 400;
    trajectoriesPerPolicy = 1;
    
    worldBounds.position = [-1.2, 0.5];
    worldBounds.velocity = [-0.07, 0.07];
    
    state0.position = -0.5;
    state0.velocity = 0;
    
    policyUb = ones(1,dim);
    policyLb = -1.*policyUb;
    
    hyper = 0.05;
    covarianceFun = @trajectoryCovariance;
    %covarianceFun = @covarianceSqExp;
    
    ret = cell(trials,1);
    %close all;
    %figure;
    
    for trial = 1:trials
        theta = zeros(bayOptSteps, dim);
        eta = zeros(bayOptSteps, 1);
        tempEta = zeros(1, trajectoriesPerPolicy);
        traj = cell(bayOptSteps, trajectoriesPerPolicy);
        
        for i=1:initialPolicies
            theta(i,:) = randTheta(policyLb, policyUb);
            for j = 1:trajectoriesPerPolicy
                [tempEta(1,j), traj{i,j}] = execPolicyMountainCar(theta(i,:), state0, executionTimeSteps, worldBounds);
            end
            eta(i,1) = mean(tempEta);
        end

        for i=initialPolicies+1:initialPolicies+bayOptSteps
            K = covarianceFun([], theta(1:i-1,:), traj(1:i-1,:), hyper);
            histogram(K);
            toMinimize = @(x) -acquisition(x, eta(1:i-1,1), theta(1:i-1,:), covarianceFun, K, hyper, traj(1:i-1,:));
            theta(i,:) = globalMin(toMinimize, policyLb, policyUb);
            for j = 1:trajectoriesPerPolicy
                [tempEta(1,j), traj{i,j}] = execPolicyMountainCar(theta(i,:), state0, executionTimeSteps, worldBounds);
            end
            eta(i,1) = mean(tempEta);

            disp(['trial : ',num2str(trial) ,...
                '  |  step : ',num2str(i-initialPolicies) ,...
                '  |  cumulative reward: ', num2str(eta(i,1))...
                ]);        
            %visMountainCar(traj,worldBounds);
        end

        ret{trial,1}.theta = theta;
        ret{trial,1}.eta = eta;
        
        save('ret.mat','ret');
    end
    
end