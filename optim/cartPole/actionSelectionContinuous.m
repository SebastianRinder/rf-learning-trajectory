function [actionNext, prob, mu] = actionSelectionContinuous(policy, state, action, errordeviation)
    feature = state; %[state, ones(size(state,1),1)];
    mu = feature * policy';
    mu(mu > 1) = 1;
    mu(mu < -1) = -1;
        
    if isempty(action)
        if size(state,1) == 1 %generate an action for simulation and save corresponding probability
            rr = randn;
            noise = rr * errordeviation;
            actionNext = mu + noise;
            prob = (-(actionNext - mu).^2)./ (2.*errordeviation.^2);

        else %only return mean of state feature for estimating distances between unknown trajectories
            actionNext = [];
            prob = [];
        end

    else %return probabilities of state feature for given actions 
        actionNext = [];
        prob = (-(action - mu).^2)./ (2.*errordeviation.^2);
    end
end

%         if isempty(action)
%             if size(state,1) == 1
%                 rr = randn;
%                 prob = normpdf(rr);
%                 noise = rr * errordeviation;
%                 actionNext = mu + noise;
%             else
%                 actionNext = [];
%                 prob = [];
%             end
%         
%         else
%             actionNext = [];
%             prob = normpdf(abs(action-mu)./errordeviation);
% 
%         end   