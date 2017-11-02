function ret = main()
    addpath('fromMatlab');
    
    %write wilson mail
           
    trials = 1;
    bayOptSteps = 400;
    initialPolicies = 10;   

    %env = 'acroBot', 'cartPole', 'mountainCar'
    %visualize = 'none', 'best', 'all'    
    opts = environmentSettings('acroBot', 'none');
    opts.trajectoriesPerPolicy = 5;
    opts.covarianceFcn = @sqExpCovariance;
%     opts.covarianceFcn = @trajectoryCovariance;

    ub = 1.*ones(1,opts.dim);
    lb = -1.*ub;

    opts.hyper = [1,1,1];
    hyperLog = [];
    
    useMaxMean = false; % if the environment is very noisy use maxMean instead of max(knownY) in the expected improvement (Brochu 2009)
    
    ret = cell(trials,1);
    for trial = 1:trials            
        opts.trajectory.data = cell(0);
        opts.trajectory.policy = [];
        
        for i=1:initialPolicies
            opts.trajectory.policy(i,:) = randPolicy(lb, ub, 1);
            for j = 1:opts.trajectoriesPerPolicy
                [tempF(1,j), opts.trajectory.data{i,j}] = objectiveFcn(opts.trajectory.policy(i,:), opts);
            end
            knownY(i,1) = mean(tempF) + randn * opts.hyper(3);
        end

        for i = initialPolicies+1:initialPolicies+bayOptSteps
            
%             opts.hyper.f = mean(std(knownX));
%             opts.hyper.l = std(knownY)/sqrt(2);
%             if isequal(opts.covarianceFcn, @sqExpCovariance)
%                 opts.hyper.f = 500;
%                 opts.hyper.l = 0.1;
%             else
%                 opts.hyper.f = 10000;
%                 opts.hyper.l = 0.1;
%             end
            
            [opts.L, opts.alpha] = preComputeK(opts.trajectory.policy, knownY, opts);            
            
            if useMaxMean
%                 negGPModel = @(testX) -gaussianProcessModel(testX, knownX, knownY, opts);
%                 [~,negMaxMean] = globalMin(negGPModel, lb, ub, false);
%                 opts.bestY = -negMaxMean;                
            else
                opts.bestY = max(knownY);
            end
            
            negAcqFcn = @(testX) -expectedImprovement(testX, opts.trajectory.policy, knownY, opts);
            opts.trajectory.policy(i,:) = globalMin(negAcqFcn, lb, ub, false);

            for j = 1:opts.trajectoriesPerPolicy
                [tempF(1,j), opts.trajectory.data{i,j}] = objectiveFcn(opts.trajectory.policy(i,:), opts);
            end
            knownY(i,1) = mean(tempF) + randn * opts.hyper(3); %Bishop 2006

            disp(['trial : ',num2str(trial) ,...
                ' | step : ',num2str(i-initialPolicies) ,...
                ' | cr: ',num2str(knownY(i,1)),...
                ' | hyper: ',num2str(opts.hyper)...
                ]);
%             plot(objective(initialPolicies+1:end));
%             pause(0.1);

            ret{trial,1}.policy = opts.trajectory.policy;
            ret{trial,1}.objective = knownY;
            ret{trial,1}.trajectory = opts.trajectory.data;
            ret{trial,1}.hyper = hyperLog;
            save('ret.mat','ret');
            
            if mod(i-initialPolicies-1,1) ~= 0 %find Hyperparameters
                if i == initialPolicies+1
                    hyperLb(1:2) = 1e-10;
                    hyperUb(1:2) = 1e10;
                else
                    hyperLb(1:2) = 10.^(log10(hyperLog(end,1:2)) - 2);
                    hyperUb(1:2) = 10.^(log10(hyperLog(end,1:2)) + 2);
                end
                hyperLb(3) = 1e-1;
                hyperUb(3) = 2e1;
                
                negHyperFcn = @(X) -findHypers(X, opts.trajectory.policy, knownY, opts);
                opts.hyper = globalMin(negHyperFcn, hyperLb, hyperUb, true);

%                 hyper = [143902, 0.00069278, 6];
                hyperLog = [hyperLog; opts.hyper];
            end
        end
    end
end