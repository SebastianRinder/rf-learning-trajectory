function [anext, prob, P] = actionSelectionMountainCarDiscrete(policy, state, action, actionList)
    s1 = state(:,1);
    s2 = state(:,2);
    feature = [state, ones(size(state,1),1), s1.^2, s2.^2, s1.*s2, s1.^2.*s2, s2.^2.*s1, s1.^3, s2.^3];
    policy = reshape(policy, [size(policy,2)/size(actionList,2), size(actionList,2)])';
        
    P = exp(feature * policy');
    P = P ./ sum(P,2);
    
    if isempty(action)
        action = actionList(1 + sum(rand(size(P,1),1) > cumsum(P,2),2))';
    end
    Pt = P';
    prob = Pt((actionList == action)');
    
    anext = [];
end