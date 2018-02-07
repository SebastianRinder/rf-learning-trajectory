function [anext, prob, P] = actionSelectionDiscrete(policy, state, action, actionList)
    feature = [state, ones(size(state,1),1)];
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