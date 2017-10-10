function [eta, traj] = execPolicyCartPole(policy, opts)
    cumReward = 0;
    opts.state0.angle = min(max(randn * 0.01, -0.01), 0.01); %at least ~0.6 degree deviation from vertical
    state = opts.state0;    
    traj = cell(opts.timeSteps,1);
    
    %policy = reshape(policy, [size(policy,2)/2, 2])';
    
    for i=1:opts.timeSteps
        [action, probs] = opts.actionSelectionFcn(policy,state,[]);
        [nextState, reward, finished] = simCartPole(state, action, opts.bounds);
        
        if i==opts.timeSteps
            if reward == 1
                reward = 1001;
            end
        end
        cumReward = cumReward + reward;
        traj{i,1}.state = state;
        traj{i,1}.action = action;
        traj{i,1}.probs = probs;
        traj{i,1}.cumReward = cumReward;        
        
        if finished
            traj(i+1:end,:) = [];
            break;
        end
        state = nextState;
    end
    eta = cumReward / opts.timeSteps;
end