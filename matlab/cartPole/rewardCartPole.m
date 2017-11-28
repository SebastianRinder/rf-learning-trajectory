function [reward, finished] = rewardCartPole(state, bounds, goalReward)
    x = state(1);
    omega = state(4);

    reward = 1;
    if outBounds(x, bounds.rewardPosition)
        reward = reward - 1;
    end
    if outBounds(omega, bounds.rewardAngle)
        reward = reward - 1;
    end
    
    if reward == 1 && goalReward
        reward = 200;
    end

    finished = false;
    if outBounds(x, bounds.position) || outBounds(omega, bounds.angle)
        finished = true;   
    end
end