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

    bayOptSteps = 100;
    executionTimeSteps = 1000;
    trajectoriesPerPolicy = 1;
    trials = 40;
    dim = 4;
    initialPolicies = 5;

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
    
%     policyUb = [7.3, 0.11, 9.8, 9.5];            %policy boundaries from random sampling
%     policyLb = [-0.2, -0.087, -2.6, -0.0082];    %with reward > 2500
    
    hyperList = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000];
    hyperList = hyperList ./ 2;
    [hyper1, hyper2] = meshgrid(hyperList);
    hyperList = [hyper1(:), hyper2(:)];
    %covarianceFun = @trajectoryCovariance;
    covarianceFun = @sqExpCovariance;
    
    trials = size(hyperList,1);
    
    if 1    
        ret = cell(trials,1);
        %close all;
        %figure;
        for trial = 1:trials        
            hyper = hyperList(trial, :);
            %hyper = [1, 1];
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
    %             figure('name',[num2str(hyper),'  ip: ',num2str(initialPolicies)]);
    %             histogram(K);
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
                    ' | cr: ', num2str(tempEta(1,m)),...
                    ' | hyper: ', num2str(hyper)...
                    ]);
                %visCartPole(traj{i,m},worldBounds);
            end

            ret{trial,1}.theta = theta;
            ret{trial,1}.eta = eta;
            ret{trial,1}.hyper = hyper;

            save('ret.mat','ret');
        end
    else
        for i=1:dim
            policy(1,i) = optimizableVariable(['policy',int2str(i)], [lb(1,i) ub(1,i)]);
        end
        fun = @(x) -objectiveFun(x, state0, executionTimeSteps, worldBounds);
        ret = bayesopt(fun,policy,'IsObjectiveDeterministic',false,...
        'AcquisitionFunctionName','expected-improvement',...
        'MaxObjectiveEvaluations',bayOptSteps);
    end
    
    
end

function [objective, trajectory] = objectiveFun(x, state0, executionTimeSteps, worldBounds)
    dim = size(x,2);
    
    nextTheta = zeros(1,dim);
    for i=1:dim
        nextTheta(1,i) = x.(['policy',int2str(i)]);
    end
    [objective, trajectory] = execPolicyCartPole(nextTheta, state0, executionTimeSteps, worldBounds);
%     global maxObjective;
    %if maxObjective < objective
%         maxObjective = objective;
%         visCartPole(trajectory,bounds);
    %end
end