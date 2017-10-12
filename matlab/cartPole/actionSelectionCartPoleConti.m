function [aNext, quadTerm] = actionSelectionCartPoleConti(policy, state, action)
    aNext = [];
    error = 1e-3;
    
    %position
    %velocity
    %acceleration
    %angle
    %angularVelocity
    feature = state(2:5);

    mu = policy * feature';
    
%     global errorRatio;
%     errorRatio = [errorRatio; abs(error/mu)];
    
    if isempty(action)
        noise = randn * error;        
        aNext = mu + noise;
        quadTerm = (-(aNext - mu)^2) / (2*error^2);
    else
        quadTerm = (-(action - mu)^2) / (2*error^2);
    end
end

