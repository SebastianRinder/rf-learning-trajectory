function [aNext, quadTerm] = actionSelectionCartPole(theta, state, action)
    aNext = [];
    quadTerm = [];
    global trackMu;
    sigma = 1e-4;
        
%     f(1,1) = state.velocity;
%     f(1,2) = state.acceleration;
%     f(1,3) = state.angle;
%     f(1,4) = state.angleVelocity;

    f(1,1) = state.position;
    f(1,2) = state.velocity;
    f(1,3) = state.angle;
    f(1,4) = state.angleVelocity;

    mu = theta * f';
    trackMu = [trackMu; mu];
    
    if isempty(action) %select action
        noise = randn * sigma;        
        aNext = mu + noise;
    else %compute quadratic Term of the probability for the next action
        quadTerm = (-(action - mu)^2) / (2*sigma^2);
    end
end

