function ret = main()
    addpath('fromMatlab');
    %bishop book maximize 6.69
    %compute chol instead of det
    %0.5*sum(log(diag(chol(Cn)
    %if p of chol is > 0 return 1e-20
    
    %write wilson mail
    
    %std of objectives for noise
    global opts;
           
    trials = 1;
    bayOptSteps = 200;
    
    initialPolicies = 10;
%     isDeterministic = 0;
    useMaxMean = 0; % if the environment is very noisy use maxMean instead of max(knownY) in the expected improvement (Brochu 2009)

    opts.trajectoriesPerPolicy = 5;
        
%     opts.covarianceFcn = 'squaredexponential';
%     opts.covarianceFcn = @sqExpCovariance;
    opts.covarianceFcn = @trajectoryCovariance;
    
%     opts.environment = 'acroBot';
    opts.environment = 'cartPole';
%     opts.environment = 'mountainCar';
    opts.visualize = 'none'; %'none', 'best', 'all
    
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
                
        opts.actionSelectionFcn = @actionSelectionAcroBot;
        opts.simFcn = @simAcroBot;
        opts.rewardFcn = @rewardAcroBot;
        opts.visFcn = @visAcroBot;        
    end
    
    
    ub = 1.*ones(1,dim);
    lb = -1.*ub;
    
    
    ret = cell(trials,1);
    for trial = 1:trials            
        opts.trajectory.data = cell(0);
        opts.trajectory.policy = [];
        
        for i=1:initialPolicies
            opts.trajectory.policy(i,:) = randPolicy(lb, ub, 1);
            for j = 1:opts.trajectoriesPerPolicy
                [tempObjective(1,j), opts.trajectory.data{i,j}] = objectiveFcn(opts.trajectory.policy(i,:), opts);
            end
            objective(i,1) = mean(tempObjective);
        end

        for i = initialPolicies+1:initialPolicies+bayOptSteps
            knownX = opts.trajectory.policy;
            knownY = objective;
            
%             opts.hyper.f = mean(std(knownX));
%             opts.hyper.l = std(knownY)/sqrt(2);
            if isequal(opts.covarianceFcn, @sqExpCovariance)
                opts.hyper.f = 500;
                opts.hyper.l = 10;
            else
                opts.hyper.f = 10000;
                opts.hyper.l = 0.1;
            end

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

            for j = 1:opts.trajectoriesPerPolicy
                [tempObjective(1,j), opts.trajectory.data{i,j}] = objectiveFcn(opts.trajectory.policy(i,:), opts);
            end
            objective(i,1) = mean(tempObjective);

            disp(['trial : ',num2str(trial) ,...
                ' | step : ',num2str(i-initialPolicies) ,...
                ' | cr: ',num2str(objective(i,1))...
                ]);
%             plot(objective(initialPolicies+1:end));
%             pause(0.1);
            
            ret{trial,1}.policy = opts.trajectory.policy;
            ret{trial,1}.negObjective = objective;
            ret{trial,1}.trajectory = opts.trajectory.data;

            save('ret.mat','ret');
        end
    end
end