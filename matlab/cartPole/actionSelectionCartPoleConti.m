function [aNext, quadTerm] = actionSelectionCartPoleConti(policy, state, action)
    aNext = [];
    quadTerm = [];
    error = 1e-3;
    
    feature(1,1) = state.velocity;
    feature(1,2) = state.acceleration;
    feature(1,3) = state.angle;
    feature(1,4) = state.angularVelocity;

    mu = policy * feature';
    
    global errorRatio;
    errorRatio = [errorRatio; abs(error/mu)];
    
    if isempty(action)
        noise = randn * error;        
        aNext = mu + noise;
    else
        quadTerm = (-(action - mu)^2) / (2*error^2);
    end
end

