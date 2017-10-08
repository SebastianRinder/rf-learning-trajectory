function [eta, traj] = execPolicyCartPole(policy, simOpts)
    cumReward = 0;
    state = simOpts.state0;    
    traj = cell(simOpts.timeSteps,1);
    
    for i=1:simOpts.timeSteps
        action = actionSelectionCartPole(policy,state,[]);
        [nextState, reward, finished] = simCartPole(state, action, simOpts);
        
        if i==simOpts.timeSteps
            reward = min(0,reward) * 500; %500 reward if in one reward bound, 1000 if in both
        end        
        cumReward = cumReward + reward;
        traj{i,1}.state = state;
        traj{i,1}.action = action;
        traj{i,1}.cumReward = cumReward;        
        
        if finished
            traj(i+1:end,:) = [];
            break;
        end
        state = nextState;
    end
    eta = cumReward / simOpts.timeSteps;
end