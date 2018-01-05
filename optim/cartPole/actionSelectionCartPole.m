function [actionNext, quad, mu] = actionSelectionCartPole(policy, state, action, opts)
    %position
    %velocity
    %acceleration
    %angle
    %angularVelocity
    feature = state(:,2:5);
    mu = feature * policy';
    mu(mu > 1) = 1;
    mu(mu < -1) = -1;
    
    if nargout < 3
        if isempty(action)        
            noise = randn * opts.errorVariance;        
            actionNext = mu + noise;
            quad = (-(actionNext - mu).^2) ./ (2.*opts.errorVariance);

        else
            actionNext = [];
            quad = (-(action - mu).^2) ./ (2.*opts.errorVariance);

        end
    else
        actionNext = [];
        quad = [];
    end
end

