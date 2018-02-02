function [anext, prob, mu] = actionSelectionMountainCarContinuous(policy, state, action, errordeviation)
    %errordeviation = 1e-1;
    s1 = state(:,1);
    s2 = state(:,2);
    feature = [state, ones(size(state,1),1), s1.^2, s2.^2, s1.*s2, s1.^2.*s2, s2.^2.*s1, s1.^3, s2.^3];
    mu = feature * policy';
    mu(mu > 1) = 1;
    mu(mu < -1) = -1;
    
    if nargout < 3
        if isempty(action)        
            noise = randn * errordeviation;
            anext = mu + noise;
            prob = (-(anext - mu).^2); %./ (2.*errordeviation.^2);

        else
            anext = [];
            prob = (-(action - mu).^2); %./ (2.*errordeviation.^2);

        end
    else
        anext = [];
        prob = [];
    end
end