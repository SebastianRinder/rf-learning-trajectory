function K = kernelTrajectory(xi, xj, traji, trajj, hyper, actionSelectionFcn)    
    trajectoriesPerPolicy = size(trajj,2);
    
    if isempty(traji)   %unknownTrajectory -> importance Sampling
        tempD = [1;0;0];
        D = 0;

        for k = 1:trajectoriesPerPolicy
            traj = trajj{1,k};

            for t=1:size(traj,1)
                state = traj{t,1}.state;
                action = traj{t,1}.action;
                [~,quadNew] = actionSelectionFcn(xi, state, action);
                [~,quadj] = actionSelectionFcn(xj, state, action);
                tempD(1,1) = tempD(1,1) * exp(quadNew - quadj);
                tempD(2,1) = tempD(2,1) + quadNew - quadj;
                tempD(3,1) = tempD(3,1) + quadj - quadNew;
            end

            D = D + tempD(1,1) * tempD(2,1) + tempD(3,1);
        end
        
        D = D / trajectoriesPerPolicy; %normalize
        
    else %known Trajectory        
        D1 = 0;
        D2 = 0;
        
        for k = 1:trajectoriesPerPolicy
            traj = traji{1,k};

            for t=1:size(traj,1)
                state = traj{t,1}.state;
                action = traj{t,1}.action;
                [~,quadi] = actionSelectionFcn(xi, state, action);
                [~,quadj] = actionSelectionFcn(xj, state, action);
                D1 = D1 + quadi - quadj;
            end
        end
        
        for k = 1:trajectoriesPerPolicy
            traj = trajj{1,k};

            for t=1:size(traj,1)
                state = traj{t,1}.state;
                action = traj{t,1}.action;
                [~,quadi] = actionSelectionFcn(xi, state, action);
                [~,quadj] = actionSelectionFcn(xj, state, action);
                D2 = D2 + quadj - quadi;
            end
        end
        
        D = (D1 + D2) / trajectoriesPerPolicy; %normalize
    end
    
    K = hyper.f .* exp(-hyper.l * D);
    %K = D;
end

