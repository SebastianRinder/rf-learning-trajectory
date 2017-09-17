function [eta, traj] = execPolicyMountainCar(theta, state0, episodes, worldBounds)
    cumReward = 0;
    state = state0;
    
    traj = cell(episodes,1);

    for i=1:episodes 
        action = actionSelectionMountainCar(theta, state, []);
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
    eta = cumReward;% / episodes;
end