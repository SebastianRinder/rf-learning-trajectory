function [aNext, quadTerm] = actionSelectionCartPoleConti(policy, state, action)
    aNext = [];
    quadTerm = [];
    global trackMu;
    sigma = 1e-2;
    
    feature(1,1) = state.velocity;
    feature(1,2) = state.acceleration;
    feature(1,3) = state.angle;
    feature(1,4) = state.angularVelocity;

%     f(1,1) = state.position;
%     f(1,2) = state.velocity;
%     f(1,3) = state.angle;
%     f(1,4) = state.angleVelocity;

    mu = policy * feature';
    trackMu = [trackMu; mu];
    
    if isempty(action)
        noise = randn * sigma;        
        aNext = mu + noise;
    else
        quadTerm = (-(action - mu)^2) / (2*sigma^2);
    end
end

