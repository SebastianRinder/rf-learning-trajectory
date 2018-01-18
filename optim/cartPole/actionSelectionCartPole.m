function [actionNext, quad, mu] = actionSelectionCartPole(policy, state, action, opts)
    %position
    %velocity
    %acceleration
    %angle
    %angularVelocity
    errordeviation = 1e-4;
    feature = state(:,2:5);
    mu = feature * policy';
    mu(mu > 1) = 1;
    mu(mu < -1) = -1;
    
    if nargout < 3
        if isempty(action)        
            noise = randn * errordeviation;
            actionNext = mu + noise;
            quad = (-(actionNext - mu).^2) ./ (2.*errordeviation);

        else
            actionNext = [];
            quad = (-(action - mu).^2) ./ (2.*errordeviation);

        end
    else
        actionNext = [];
        quad = [];
    end
end

