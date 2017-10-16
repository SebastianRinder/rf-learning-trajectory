function [finalReward, traj] = execPolicyMountainCar(policy, opts)
    traj.state = zeros(opts.timeSteps,5);
    traj.action = zeros(opts.timeSteps,1);
    traj.prob = zeros(opts.timeSteps,1);
    traj.cumReward = zeros(opts.timeSteps,1);
    
    cumReward = 0;
    state = opts.state0;
    
    for i=1:opts.timeSteps 
        [action, prob] = opts.actionSelectionFcn(policy, state, []);
        [nextState, reward, finished] = simMountainCar(state, action, worldBounds);
   
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
    finalReward = cumReward;% / episodes;
end