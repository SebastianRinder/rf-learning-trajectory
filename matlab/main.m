function ret = main()
    addpath('fromMatlab');
    %one slide motivation
    %two slides bo
    %two/three slides traj kernel
    %two slides experiments
    %one slide future work
               
    trials = 40;
    bayOptSteps = 40;
    initialPolicies = 10;

    %env = 'acroBot', 'cartPole', 'mountainCar'
    %visualize = 'none', 'best', 'all'
    opts = environmentSettings('cartPole', 'none');
    opts.trajectoriesPerPolicy = 5;
%     opts.covarianceFcn = @sqExpCovariance;
    opts.covarianceFcn = @trajectoryCovariance;
    opts.debugPlotting = 0;

    ub = 1.*ones(1,opts.dim);
    lb = -1.*ub;

    opts.hyper = [1,1];
    opts.hyperN = 1;
        
    useMaxMean = 0; % if the environment is very noisy use maxMean instead of max(knownY) in the expected improvement (Brochu 2010)
    
    ret = cell(trials,1);
    for trial = 1:trials
        tic
        opts.trajectory.data = cell(0);
        opts.trajectory.policy = [];
        allY = [];
        knownY = [];
        hyperTrace = [];
        
        for i=1:initialPolicies
            opts.trajectory.policy(i,:) = randPolicy(lb, ub, 1);
            for j = 1:opts.trajectoriesPerPolicy
                [allY(i,j), opts.trajectory.data{i,j}] = objectiveFcn(opts.trajectory.policy(i,:), opts);
            end
%             errorNoise = randn(1, opts.trajectoriesPerPolicy) .* opts.hyperN;
%             knownY(i,1) = mean(allF(i,:) + errorNoise); %Bishop 2006
            knownY(i,1) = mean(allY(i,:));
        end

        for i = initialPolicies+1:initialPolicies+bayOptSteps            
            opts.D = opts.covarianceFcn(opts.trajectory.policy, opts.trajectory.policy, opts);
            
            if mod(i-initialPolicies-1,1) == 0 %find Hyperparameters
                if i == initialPolicies+1
                    hyperLb(1:2) = -10;
                    hyperUb(1:2) = 10;
                else
                    hyperLb(1:2) = log10(hyperTrace(end,1:2)) - 2;
                    hyperUb(1:2) = log10(hyperTrace(end,1:2)) + 2;
                end
                
                negHyperFcn = @(X) -findLogHypers(X, opts.trajectory.policy, knownY, opts);
                log10Hyper = globalMin(negHyperFcn, hyperLb, hyperUb, true, opts.debugPlotting);
                opts.hyper = 10.^(log10Hyper);
            end
            
            opts.hyperN = mean(std(allY,0,2).^2);
            if opts.hyperN == 0
                opts.hyperN = 1;
            end
            hyperTrace = [hyperTrace; opts.hyper, opts.hyperN];
                        
            [opts.L, opts.alpha] = getLowerCholesky(opts.D, knownY, opts);            
            
            if useMaxMean
                negGPModel = @(testX) -gaussianProcessModel(testX, opts.trajectory.policy, knownY, opts);
                [~,negMaxMean] = globalMin(negGPModel, lb, ub, false, opts.debugPlotting);
                opts.bestY = -negMaxMean;
            else
                opts.bestY = max(knownY);
            end           
                        
            negAcqFcn = @(testX) -expectedImprovement(testX, opts.trajectory.policy, knownY, opts);
            opts.trajectory.policy(i,:) = globalMin(negAcqFcn, lb, ub, false, opts.debugPlotting);

            for j = 1:opts.trajectoriesPerPolicy
                [allY(i,j), opts.trajectory.data{i,j}] = objectiveFcn(opts.trajectory.policy(i,:), opts);
            end
%             errorNoise = randn(1, opts.trajectoriesPerPolicy) .* opts.hyperN;
%             knownY(i,1) = mean(allF(i,:) + errorNoise); %Bishop 2006
            knownY(i,1) = mean(allY(i,:));
            if knownY(i,1) > 1300
                disp('crfail');
            end
            disp(['trial : ',num2str(trial) ,...
                ' | step : ',num2str(i-initialPolicies) ,...
                ' | cr: ',num2str(knownY(i,1)),...
                ' | hyper: ',num2str(hyperTrace(end,:))...
                ]);
%             plot(objective(initialPolicies+1:end));
%             pause(0.1);

            ret{trial,1}.policy = opts.trajectory.policy;
            ret{trial,1}.objective = knownY;
            ret{trial,1}.trajectory = opts.trajectory.data;
            ret{trial,1}.hyper = hyperTrace;
                        
            if ~isequal(opts.visualize, 'none')
                [~,bestTrajectoryIdx] = max(allY(i,:));
                [~,bestObjectiveIdx] = max(knownY);
                if isequal(opts.visualize, 'best') && i == bestObjectiveIdx ||...
                        isequal(opts.visualize, 'all')
                    traj = opts.trajectory.data{i,bestTrajectoryIdx};
                    opts.visFcn(traj,opts.bounds);
                end
            end
            save('retCartpoleTraj.mat','ret');
            if i >= initialPolicies+20 && max(knownY) < 200
                break;
            end
        end
        
        disp(num2str(toc/60));
    end
end