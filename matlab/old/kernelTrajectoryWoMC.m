function K = kernelTrajectoryWoMC(thetai, thetaj, traji, trajj, hyper)
    global execPolicyFun;
    global actionSelectionFun;
    
    global state0;
    global executionTimeSteps;
    global worldBounds;
    global trajectoriesPerPolicy;
    
    if isempty(traji)
        traji = cell(1,trajectoriesPerPolicy);

        for j = 1:trajectoriesPerPolicy
            [~, traji{1,j}] = execPolicyFun(thetai, state0, executionTimeSteps, worldBounds);
        end
    end
    
    D1 = 0; 
    
    for k = 1:size(traji,2)
        traj = traji{1,k};
        tempPi = 1;
        tempPj = 1;
        
        for t=1:size(traj,1)
            state = traj{t,1}.state;
            action = traj{t,1}.action;

            Pi = actionSelectionFun(thetai,state, action);
            Pj = actionSelectionFun(thetaj,state, action);

            tempPi = tempPi * Pi;
            tempPj = tempPj * Pj;
        end
        
        D1 = D1 + tempPi * log(tempPi/tempPj);
    end
    
    D2 = 0;
    
    for k = 1:size(trajj,2)
        traj = trajj{1,k};
        tempPi = 1;
        tempPj = 1;
        
        for t=1:size(traj,1)
            state = traj{t,1}.state;
            action = traj{t,1}.action;

            Pi = actionSelectionFun(thetai,state, action);
            Pj = actionSelectionFun(thetaj,state, action);

            tempPi = tempPi * Pi;
            tempPj = tempPj * Pj;
        end
        
        D2 = D2 + tempPj * log(tempPj/tempPi);
    end
    
    if D1<0 
        D1 = 0;
    end
    if D2<0 
        D2 = 0;
    end
    
    D = sqrt(D1)+sqrt(D2);
    
    K = exp(-hyper * D);
end

