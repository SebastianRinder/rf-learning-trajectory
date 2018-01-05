function [aNext, prob] = actionSelectionAcroBot(policy, state, action, actionList)
    feature = state;
    policy = reshape(policy, [size(policy,2)/size(actionList,2), size(actionList,2)])';    
    
    P = exp(feature * policy');
    P = P ./ sum(P,2);
    
    if isempty(action)        
        aNext = actionList(1 + sum(rand > cumsum(P)));
        prob = P(actionList == aNext);
        
    else
        aNext = [];
        prob = sum((actionList == action) .* P, 2);
    end
end