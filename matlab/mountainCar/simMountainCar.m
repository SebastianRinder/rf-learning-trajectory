function sNext = simMountainCar(state, action, bounds)
    p = state(1);
    v = state(2);
    
    vNext = v + 0.001 * action - 0.0025 * cos(3 * p);
    
    % Apply velocity boundary
    vNext = applyBound(vNext, bounds.velocity);
    
    sNext(1,1) = p + vNext;
    sNext(1,2) = vNext;
    
    if state(1) < bounds.position(1) %out of left bound, inelastic wall to the left
        sNext(1) = bounds.position(1);
        sNext(2) = 0;
    end
end