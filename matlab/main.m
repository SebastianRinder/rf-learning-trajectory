function ret = main()
    addpath('fromMatlab');
    
    global opts;
    
    trials = 5;
    bayOptSteps = 250;
    
    initialPolicies = 40;
    isDeterministic = 0;
    useMatlabBayes = 1;
    useMatlabGP = 1;    
    optimizeNoiseHyper = 0;
    acqFcn = @expectedImprovement;
    opts.trajectoriesPerPolicy = 10;
        
    %opts.covarianceFcn = 'squaredexponential';
    %opts.covarianceFcn = @sqExpCovariance;
    opts.covarianceFcn = @trajectoryCovariance;
        
    %opts.problem = 'cartPole';
    opts.problem = 'mountainCar';
    addpath(opts.problem);
    
    if isequal(opts.problem, 'cartPole')
        dim = 4;
        opts.bounds.position = [-5, 5];
        opts.bounds.angle = [-90 * pi / 180, 90 * pi / 180];
        opts.bounds.rewardPosition = [-1, 1];
        opts.bounds.rewardAngle = [-12 * pi / 180, 12 * pi / 180];
        opts.state0 = zeros(1,5); %position, velocity, acceleration, angle, angularVelocity
        opts.timeSteps = 1000;
    
        %opts.execPolicyFcn = @execPolicyCartPole;
        opts.actionSelectionFcn = @actionSelectionCartPoleConti;
        opts.simFcn = @simCartPole;
    elseif isequal(opts.problem, 'mountainCar')
        dim = 18;        
        opts.bounds.position = [-1.2, 0.5];
        opts.bounds.velocity = [-0.07, 0.07];
        opts.state0 = [-0.5, 0]; %position, velocity
        opts.timeSteps = 400;
        
        %opts.execPolicyFcn = @execPolicyMountainCar;
        opts.actionSelectionFcn = @actionSelectionMountainCar;
        opts.simFcn = @simMountainCar;
    end
    
    ub = 1.*ones(1,dim);
    lb = -1.*ub;
    
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
        
        ret = zeros(trials,4);
        for trial=1:trials
            errorRatio = [];
            close all;
            
            for i=1:dim
                policy(1,i) = optimizableVariable(['policy',int2str(i)], [lb(1,i) ub(1,i)]);
            end 
            
            tic;
            tempRet = customBayesopt(@negObjectiveFcn, policy, opts,...
                'IsObjectiveDeterministic',isDeterministic,...
                'AcquisitionFunctionName','expected-improvement',...
                'MaxObjectiveEvaluations',bayOptSteps,...
                'NumSeedPoints',initialPolicies);
            toc;
            disp(['minutes: ', num2str(toc/60)]);

            bestIter = find(tempRet.MinObjective == tempRet.ObjectiveTrace);
            
            ret(trial,1) = mean(errorRatio);
            ret(trial,2) = toc/60;
            ret(trial,3) = tempRet.MinObjective;
            ret(trial,4) = bestIter(1);
            
%             ret{trial,1}.minutes = toc/60;
%             ret{trial,1}.MinObjective = tempRet.MinObjective;
%             ret{trial,1}.bestIter = bestIter(1);
            
            save('ret.mat','ret');
        end
    end
end

function [negObjective, constraints, trajectory] = negObjectiveFcn(x)
    global opts;
    dim = size(x,2);
    
    policy = zeros(1,dim);
    for i=1:dim
        policy(1,i) = x.(['policy',int2str(i)]);
    end
    for i = opts.trajectoriesPerPolicy:-1:1
        [tempObjective(1,i), trajectory.data{1,i}] = execPolicy(policy, opts);
    end
    
    negObjective = -mean(tempObjective);
    constraints = [];
    trajectory.policy = policy;
end