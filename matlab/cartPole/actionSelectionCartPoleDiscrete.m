function [aNext, probs] = actionSelectionCartPole(policy, state, actionList)

%     f(1,1) = state.velocity;
%     f(1,2) = state.acceleration;
%     f(1,3) = state.angle;
%     f(1,4) = state.angleVelocity;

    feature(1,1) = state.position;
    feature(1,2) = state.velocity;
    feature(1,3) = state.angle;
    feature(1,4) = state.angleVelocity;
    
    %probs = zeros(size(actionList,2),1);
    
    probs = exp(policy * feature');
    
    %normalize
    probs = probs ./ sum(probs);    

    %select action according to probabilities
    aNext = actionList(1 + sum(rand > cumsum(probs)));

    %probability for given action
%     if ~isempty(action)
%         P = probs(actionList == action);
%     end
end