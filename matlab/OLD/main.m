function ret = main()
    addpath('fromMatlab');
    
    global opts;
    global failz;
       
    trials = 1;
    bayOptSteps = 200;
    
    initialPolicies = 10;
    isDeterministic = 0;
    useMatlabBayes = 1;
    useMatlabGP = 0;    
    optimizeNoiseHyper = 0;
    acqFcn = @expectedImprovement;
    opts.trajectoriesPerPolicy = 5;
        
    opts.covarianceFcn = 'squaredexponential';
    %opts.covarianceFcn = @sqExpCovariance;
%     opts.covarianceFcn = @trajectoryCovariance;
    
    opts.environment = 'acroBot';
%     opts.environment = 'cartPole';
%     opts.environment = 'mountainCar';
    opts.visualize = false;
    
    addpath(opts.environment);    
    
    if isequal(opts.environment, 'cartPole')
        dim = 4;
        opts.actionList = [];	%continuous action selection
        opts.bounds.position = [-5, 5];
        opts.bounds.angle = [-90 * pi / 180, 90 * pi / 180];
        opts.bounds.rewardPosition = [-1, 1];
        opts.bounds.rewardAngle = [-12 * pi / 180, 12 * pi / 180];
        opts.state0 = zeros(1,5); %position, velocity, acceleration, angle, angularVelocity
        opts.timeSteps = 1000;
    
        opts.actionSelectionFcn = @actionSelectionCartPole;
        opts.simFcn = @simCartPole;
        opts.rewardFcn = @rewardCartPole;
        opts.visFcn = @visCartPole;
    elseif isequal(opts.environment, 'mountainCar')
        dim = 9*2;
        opts.actionList = [-1,1];   %apply acceleration to the rear or forward
        opts.bounds.position = [-1.2, 0.5];
        opts.bounds.velocity = [-0.07, 0.07];
        opts.state0 = [-0.5, 0]; %position, velocity
        opts.timeSteps = 400;
                
        opts.actionSelectionFcn = @actionSelectionMountainCar;
        opts.simFcn = @simMountainCar;
        opts.rewardFcn = @rewardMountainCar;
        opts.visFcn = @visMountainCar;
    elseif isequal(opts.environment, 'acroBot')
        dim = 4*3;
        opts.actionList = [-1,0,1];       %apply torque to hip
        opts.bounds.angle1 = [-pi, pi];
        opts.bounds.angle2 = [-pi, pi];
        opts.bounds.velocity1 = [-4*pi, 4*pi];
        opts.bounds.velocity2 = [-9*pi, 9*pi];
        opts.state0 = [0, 0, 0, 0]; %angle1, angle2, angularVelocity1, angularVelocity2
        opts.timeSteps = 400;
        %opts.actionStep = 4;
        
        opts.actionSelectionFcn = @actionSelectionAcroBot;
        opts.simFcn = @simAcroBot;
        opts.rewardFcn = @rewardAcroBot;
        opts.visFcn = @visAcroBot;        
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
        for trial=1:trials
            failz = 0;
            for i=1:dim
                policy(1,i) = optimizableVariable(['policy',int2str(i)], [lb(1,i) ub(1,i)]);
            end 
            
            tic;
            close all;
            tempRet = customBayesopt(@negObjectiveFcn, policy, opts,...
                'IsObjectiveDeterministic',isDeterministic,...
                'AcquisitionFunctionName','expected-improvement',...
                'MaxObjectiveEvaluations',bayOptSteps,...
                'NumSeedPoints',initialPolicies);
            toc;
            disp(['minutes: ', num2str(toc/60)]);

            bestIter = find(tempRet.MinObjective == tempRet.ObjectiveTrace);
            
            ret(trial,1) = toc/60;
            ret(trial,2) = tempRet.MinObjective;
            ret(trial,3) = bestIter(1);
            ret(trial,4) = mean(tempRet.ObjectiveTrace);
            ret(trial,5) = failz;
            
            trajSave(trial,1) = tempRet.UserDataTrace(bestIter(1),1);
            
            disp(num2str(trial));
            if tempRet.MinObjective == bestIter(1)
                opts.visFcn(tempRet.UserDataTrace{bestIter(1),1},opts.bounds);
            end
            save('ret.mat','ret', 'trajSave');
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