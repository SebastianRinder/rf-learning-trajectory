function [reward, finished] = rewardMountainCar(state, bounds, ~)
    finished = false;
    
    if state(1) < bounds.position(2)
        reward = -1;          
    else %reached goal: out of right bound
        reward = 100;
        finished = true;
    end
end

