function [sNext, reward, finished] = simMountainCar(state, action, bounds)
    p = state(1);
    v = state(2);
    
    vNext = v + 0.001 * action - 0.0025 * cos(3 * p);
    
    % Apply velocity boundary
    vNext = min(max(vNext, bounds.velocity(1)), bounds.velocity(2));
    
    sNext(1,1) = p + vNext;
    sNext(1,2) = vNext;
    finished = false;
    
    if sNext(1,1) < bounds.position(1) %out of left bound, inelastic wall to the left
        sNext(1,1) = bounds.position(1);
        sNext(1,2) = 0;
    end
    
    if sNext(1,1) < bounds.position(2)
        reward = -1;          
    else %reached goal: out of right bound
        reward = 1;
        finished = true;
    end
end