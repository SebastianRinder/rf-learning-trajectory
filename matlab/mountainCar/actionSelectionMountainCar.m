function [aNext, prob] = actionSelectionMountainCar(policy, state, action)
    actionList = [-1,1];       %apply acceleration to the rear or forward
    
    pos = state(:,1);
    vel = state(:,2);
    feature = [pos, pos.^2, pos.^3, vel, vel.^2, vel.^3, pos.*vel, pos.*vel.^2, pos.^2.*vel]; 
    
    policy = reshape(policy, [size(policy,2)/size(actionList,2), size(actionList,2)])';
    
    P = exp(feature * policy');
    P = P ./ sum(P,2);
    
    %google: logsumexp matlab
    
    if isempty(action)        
        aNext = actionList(1 + sum(rand > cumsum(P)));
        prob = P(actionList == aNext);
        
    else
        aNext = [];
        prob = sum((actionList == action) .* P, 2);
        
    end
end