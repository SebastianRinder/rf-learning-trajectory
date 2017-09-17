function ret = mainCartPole()
    addpath('..');
    
    global execPolicyFun;
    global actionSelectionFun;    
    execPolicyFun = @execPolicyCartPole;
    actionSelectionFun = @actionSelectionCartPole;
    
    global state0;
    global executionTimeSteps;
    global worldBounds;
    global trajectoriesPerPolicy;

    bayOptSteps = 140;
    executionTimeSteps = 1000;
    trajectoriesPerPolicy = 1;
    trials = 40;
    dim = 4;
    initialPolicies = 50;

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
    
    ub = ones(1,dim);
    lb = -1.*ub;
    
%     policyUb = [7.3, 0.11, 9.8, 9.5];            %policy boundaries from random sampling
%     policyLb = [-0.2, -0.087, -2.6, -0.0082];    %with reward > 2500
    
    hyper = 1e-0;
    covarianceFun = @trajectoryCovariance;
    %c%ovarianceFun = @sqExpCovariance;
    
    ret = cell(trials,1);
    %close all;
    %figure;
    for trial = 1:trials        
        theta = zeros(bayOptSteps, dim);
        eta = zeros(bayOptSteps, 1);
        tempEta = zeros(1, trajectoriesPerPolicy);
        traj = cell(bayOptSteps, trajectoriesPerPolicy);

        for i=1:initialPolicies
            theta(i,:) = randTheta(lb, ub);
            for j = 1:trajectoriesPerPolicy
                [tempEta(1,j), traj{i,j}] = execPolicyCartPole(theta(i,:), state0, executionTimeSteps, worldBounds);
            end
            %histogram(tempEta);
            eta(i,1) = max(tempEta);
        end
            
        for i = initialPolicies+1:initialPolicies+bayOptSteps
            K = covarianceFun([], theta(1:i-1,:), traj(1:i-1,:), hyper);
            figure('name',[num2str(hyper),'  ip: ',num2str(initialPolicies)]);
            histogram(K);
            toMinimize = @(x) -acquisition(x, eta(1:i-1,1), theta(1:i-1,:), covarianceFun, K, hyper, traj(1:i-1,:));
            theta(i,:) = globalMin(toMinimize, theta(i-1,:), lb, ub);
            for j = 1:trajectoriesPerPolicy
                [tempEta(1,j), traj{i,j}] = execPolicyCartPole(theta(i,:), state0, executionTimeSteps, worldBounds);
            end
            eta(i,1) = max(tempEta);
            
            m = find(tempEta == max(tempEta));            
            m = m(1);
            disp(['trial : ',num2str(trial) ,...
                ' | step : ',num2str(i-initialPolicies) ,...
                ' | cr: ', num2str(tempEta(1,m)...
                )]);
            %visCartPole(traj{i,m},worldBounds);
        end

        ret{trial,1}.theta = theta;
        ret{trial,1}.eta = eta;
        
        save('ret.mat','ret');
    end
    
end