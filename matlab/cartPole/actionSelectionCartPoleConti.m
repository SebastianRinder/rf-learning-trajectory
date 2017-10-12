function [aNext, quad] = actionSelectionCartPoleConti(policy, state, action)
    aNext = [];
    error = 1e-3;
    
    %position
    %velocity
    %acceleration
    %angle
    %angularVelocity
    
    
%     global errorRatio;
%     errorRatio = [errorRatio; abs(error/mu)];

    feature = state(:,2:5);
    mu = feature * policy';
    
    if isempty(action)
        
        noise = randn * error;        
        aNext = mu + noise;
        quad = (-(aNext - mu)^2) / (2*error^2);
    else
        
        
        quad = (-(action - mu).^2) ./ (2*error^2);
    end
end

