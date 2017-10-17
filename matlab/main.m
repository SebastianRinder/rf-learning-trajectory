function ret = main()
    addpath('fromMatlab');
    
    global opts;
    
    trials = 5;
    bayOptSteps = 250;
    
    initialPolicies = 10;
    isDeterministic = 0;
    useMatlabBayes = 0;
    useMatlabGP = 0;    
    optimizeNoiseHyper = 0;
    acqFcn = @expectedImprovement;
    opts.trajectoriesPerPolicy = 5;
        
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
            opts.trajectory.data = cell(0);
            opts.trajectory.policy = [];

            for i=1:initialPolicies
                opts.trajectory.policy(i,:) = randPolicy(lb, ub, 1);
                for j = 1:opts.trajectoriesPerPolicy
                    [tempObjective(1,j), opts.trajectory.data{i,j}] = execPolicy(opts.trajectory.policy(i,:), opts);
                end
                negObjective(i,1) = -mean(tempObjective);
            end
           
            for i = initialPolicies+1:initialPolicies+bayOptSteps
                xTrain = opts.trajectory.policy;
                yTrain = negObjective;
                hyper = setHypers(yTrain, lb, ub, isDeterministic);
                                
                if useMatlabGP
                    kernelFcn = @(Xm, Xn, theta) opts.covarianceFcn(Xm, Xn, theta, opts);
                    gprMdl = fitGP(xTrain, yTrain, hyper, kernelFcn, optimizeNoiseHyper);
                    meanMdlFcn = @(X) predict(gprMdl, X);
                    [~,minMeanMdl] = globalMin(meanMdlFcn, lb, ub, useMatlabGP);

                    negAcqFcn = @(X) -acqFcn(X, gprMdl, minMeanMdl);
                    opts.trajectory.policy(i,:) = globalMin(negAcqFcn, lb, ub, useMatlabGP);
                else
                    K = opts.covarianceFcn(xTrain, xTrain, hyper, opts);
                    meanMdlFcn = @(X) gaussianProcess(X, xTrain, yTrain, opts.covarianceFcn, K, hyper, opts);
                    [~,minMeanMdl] = globalMin(meanMdlFcn, lb, ub, useMatlabGP);
                    
                    negAcqFcn = @(X) -acquisition(X, xTrain, yTrain, minMeanMdl, opts.covarianceFcn, K, hyper, opts);
                    opts.trajectory.policy(i,:) = globalMin(negAcqFcn, lb, ub, useMatlabGP);
                end
                
                for j = 1:opts.trajectoriesPerPolicy
                    [tempObjective(1,j), opts.trajectory.data{i,j}] = execPolicy(opts.trajectory.policy(i,:), opts);
                end
                negObjective(i,1) = -mean(tempObjective);

                disp(['trial : ',num2str(trial) ,...
                    ' | step : ',num2str(i-initialPolicies) ,...
                    ' | negCR: ',num2str(negObjective(i,1))...
                    ]);

                ret{trial,1}.policy = opts.trajectory.policy;
                ret{trial,1}.negObjective = negObjective;
                ret{trial,1}.trajectory = opts.trajectory.data;

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