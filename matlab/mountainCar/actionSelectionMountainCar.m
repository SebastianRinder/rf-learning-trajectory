function [aNext, P] = actionSelectionMountainCar(theta, state, action)
    actionList = [-1,1];       %apply acceleration to the rear or forward
    
    theta = reshape(theta, [size(theta,2)/size(actionList,2), size(actionList,2)])';
    
    p = state.position;
    v = state.velocity;
    feature = [p; p^2; p^3; v; v^2; v^3; p*v; p*v^2; p^2*v]; 
       
    sumP = 0;
    probs = zeros(size(actionList,2),1);
    for i=1:size(actionList,2)
        probs(i,1) = exp(theta(i,:) * feature);
        sumP = sumP + probs(i,1);
    end
    
    %normalize
    for i=1:size(actionList,2)
        probs(i,1) = probs(i,1) / sumP;
    end

    %select action according to probabilities
    aNext = actionList(1+sum(rand>cumsum(probs)));

    %probability for given action
    if ~isempty(action)
        P = probs(actionList == action);
    end
end