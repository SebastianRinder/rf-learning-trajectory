function [sNext, reward, finished] = simMountainCar(state, action, worldBounds)
    p = state.position;
    v = state.velocity;
    
    vNext = v + 0.001 * action - 0.0025 * cos(3 * p);
    
    % Apply velocity boundary
    vNext = min(max(vNext, worldBounds.velocity(1,1)), worldBounds.velocity(1,2));
    
    sNext.position = p + vNext;
    sNext.velocity = vNext;
    finished = false;
    
    if sNext(1,1).position < worldBounds.position(1,1) %out of left bound, inelastic wall to the left
        sNext.position = worldBounds.position(1,1);
        sNext.velocity = 0;
    end
    
    if sNext(1,1).position < worldBounds.position(1,2)
        reward = -1;          
    else %reached goal: out of right bound
        finished = true;
        reward = 1000;
    end
end