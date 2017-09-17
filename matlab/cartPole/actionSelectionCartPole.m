function [aNext, quad] = actionSelectionCartPole(theta, state, action)
    sigma = 0.1;
        
%     f(1,1) = state.velocity;
%     f(1,2) = state.acceleration;
%     f(1,3) = state.angle;
%     f(1,4) = state.angleVelocity;

    f(1,1) = state.position;
    f(1,2) = state.velocity;
    f(1,3) = state.angle;
    f(1,4) = state.angleVelocity;


    noise = randn * sigma;
    mu = theta * f';
    aNext = mu + noise;
    
    %probability for given action
    if ~isempty(action)
%         normalize = normpdf(0, 0, sigma);
%         P = normpdf(action, mu, sigma) / normalize;
        quad = -((action - mu)/(2*sigma))^2;
    end
end

