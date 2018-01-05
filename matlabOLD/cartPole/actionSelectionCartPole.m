function [actionNext, quad] = actionSelectionCartPole(policy, state, action, ~)
    errorVariance = 1e-4;
       
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
        noise = randn * errorVariance;        
        actionNext = mu + noise;
        quad = (-(actionNext - mu).^2) ./ (2.*errorVariance);
        
    else
        actionNext = [];
        quad = (-(action - mu).^2) ./ (2.*errorVariance);
        
    end
end

