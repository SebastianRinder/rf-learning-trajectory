function [actionNext, quad] = actionSelectionCartPole(policy, state, action, ~)
    error = 1e-4;
       
    %position
    %velocity
    %acceleration
    %angle
    %angularVelocity
    feature = state(:,2:5);
    mu = feature * policy';
    mu(mu > 1) = 1;
    mu(mu < -1) = -1;
    
    if isempty(action)        
        noise = randn * error;        
        actionNext = mu + noise;
        quad = (-(actionNext - mu).^2) ./ (2.*error.^2);
        
    else
        actionNext = [];
        quad = (-(action - mu).^2) ./ (2.*error.^2);
        
    end
end

