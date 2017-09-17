function K = kernelTrajectoryWoIS(thetai, thetaj, traji, trajj, hyper)

    global execPolicyFun;
    global actionSelectionFun;
    
    global state0;
    global executionTimeSteps;
    global worldBounds;
    global trajectoriesPerPolicy;
    
    D = 0;
    toBreak = false;
        
    if isempty(traji)   %unknown trajectory -> sample trajectories
        traji = cell(1,trajectoriesPerPolicy);

        for j = 1:trajectoriesPerPolicy
            [~, traji{1,j}] = execPolicyFun(thetai, state0, executionTimeSteps, worldBounds);
        end
    end

    for k = 1:size(traji,2)
        traj = traji{1,k};

        for t=1:size(traj,1)
            state = traj{t,1}.state;
            action = traj{t,1}.action;

            [~, quadi] = actionSelectionFun(thetai,state, action);            
            [~, quadj] = actionSelectionFun(thetaj,state, action);

%             if Pj == 0
%                 toBreak = true;
%                 break;
%             end
            
            D = D + quadi - quadj;
        end
%         if toBreak
%             break;
%         end
    end
    
    if ~toBreak
        for k = 1:size(trajj,2)
            traj = trajj{1,k};

            for t=1:size(traj,1)
                state = traj{t,1}.state;
                action = traj{t,1}.action;

                [~, quadi] = actionSelectionFun(thetai,state, action);
                [~, quadj] = actionSelectionFun(thetaj,state, action);

%                 if Pi == 0
%                     toBreak = true;
%                     break;
%                 end

                D = D + quadj - quadi;
            end
%             if toBreak
%                 break;
%             end
        end
    end
    
    if D<0
        %D = 0;
    end
    
    if toBreak
        K = 0;
    else
        K = exp(-hyper * D);
    end
end

