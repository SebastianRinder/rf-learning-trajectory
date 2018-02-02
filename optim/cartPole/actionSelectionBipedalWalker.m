function [actionNext, prob, mu] = actionSelectionBipedalWalker(policy, state, action, ~)
    errordeviation = 1e-4;
    feature = [state, ones(size(state,1),1)];
    feature = repmat(feature,1,4);
    mu = feature * policy';
    mu(mu > 1) = 1;
    mu(mu < -1) = -1;
    
    if nargout < 3
        if isempty(action)        
            noise = randn * errordeviation;
            actionNext = mu + noise;
            prob = (-(actionNext - mu).^2) ./ (2.*errordeviation.^2);

        else
            actionNext = [];
            prob = (-(action - mu).^2) ./ (2.*errordeviation.^2);

        end
    else
        actionNext = [];
        prob = [];
    end
end
