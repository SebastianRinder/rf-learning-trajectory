function [actionNext, prob, mu] = actionSelectionContinuous(policy, state, action, errordeviation)
    feature = state; %[state, ones(size(state,1),1)];
    mu = feature * policy';
    mu(mu > 1) = 1;
    mu(mu < -1) = -1;
    
    if nargout < 3
        if isempty(action)        
            noise = randn * errordeviation;
            actionNext = mu + noise;
            prob = (-(actionNext - mu).^2);  %./ (2.*errordeviation.^2);

        else
            actionNext = [];
            prob = (-(action - mu).^2); %./ (2.*errordeviation.^2);

        end
    else
        actionNext = [];
        prob = [];
    end
end
