function [anext, prob, P] = actionSelectionMountainCarDiscrete(policy, state, action, actionList)
    s1 = state(:,1);
    s2 = state(:,2);
    feature = [state, ones(size(state,1),1), s1.^2, s2.^2, s1.*s2, s1.^2.*s2, s2.^2.*s1, s1.^3, s2.^3];
    policy = reshape(policy, [size(policy,2)/size(actionList,2), size(actionList,2)])';
        
    P = exp(feature * policy');
    P = P ./ sum(P,2);
        
    if nargout < 3
        Pt = P';
        prob = Pt((actionList == action)');
    else
        prob = [];
    end    
    anext = [];
end


% function [actionNext, prob, mu] = actionSelectionCartPole(policy, state, action, ~)
%     %position
%     %velocity
%     %acceleration
%     %angle
%     %angularVelocity
%     errordeviation = 1e-4;
%     feature = state(:,2:5);
%     mu = feature * policy';
%     mu(mu > 1) = 1;
%     mu(mu < -1) = -1;
%     
%     if nargout < 3
%         if isempty(action)        
%             noise = randn * errordeviation;
%             actionNext = mu + noise;
%             prob = (-(actionNext - mu).^2) ./ (2.*errordeviation);
% 
%         else
%             actionNext = [];
%             prob = (-(action - mu).^2) ./ (2.*errordeviation);
% 
%         end
%     else
%         actionNext = [];
%         prob = [];
%     end
% end
