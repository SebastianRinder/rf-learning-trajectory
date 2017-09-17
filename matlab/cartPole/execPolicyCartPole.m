function [eta, traj] = execPolicyCartPole(theta, state0, episodes, worldBounds)
    cumReward = 0;
    state = state0;    
    traj = cell(episodes,1);
    
    for i=1:episodes
        action = actionSelectionCartPole(theta,state,[]);
        [nextState, reward, toBreak] = simCartPole(state, action, worldBounds);
        
        if i==episodes
            reward = reward + 1000;
        end        
        cumReward = cumReward + reward;
        traj{i,1}.state = state;
        traj{i,1}.action = action;
        traj{i,1}.cumReward = cumReward;        
        
        if toBreak
            traj(i+1:end,:) = [];
            break;
        end
        state = nextState;
    end
    eta = cumReward; % / episodes;
end