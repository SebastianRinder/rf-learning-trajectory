function K = kernelKnownTrajectoryCartPole(thetai, thetaj, traji, trajj, hyper)    
%     [~, ~, sigma] = actionSelectionCartPole([], []);
%     normalize = normpdf(0, 0, sigma);
    D = 0;
    toBreak = false;
    
    for k = 1:size(traji,2)
        traj = traji{1,k};
        
        for t=1:size(traj,1)
            state = traj{t,1}.state;
            action = traj{t,1}.action;

            [~, mu, sigma] = actionSelectionCartPole(thetai,state);            
            Pi = normpdf(action, mu, sigma); % / normalize;

            [~, mu, ~] = actionSelectionCartPole(thetaj,state);
            Pj = normpdf(action, mu, sigma); % / normalize;
            
            D = D + log(Pi/Pj);
            
            if Pj == 0
                toBreak = true;
                break;
            end
        end
        if toBreak
            break;
        end
    end
    
    toBreak = false;
    
    for k = 1:size(trajj,2)
        traj = trajj{1,k};
        
        for t=1:size(traj,1)
            state = traj{t,1}.state;
            action = traj{t,1}.action;

            [~, mu, ~] = actionSelectionCartPole(thetai,state);
            Pi = normpdf(action, mu, sigma); % / normalize;

            [~, mu, ~] = actionSelectionCartPole(thetaj,state);
            Pj = normpdf(action, mu, sigma); % / normalize;

            D = D + log(Pj/Pi);
            
            if Pi == 0
                toBreak = true;
                break;
            end
        end
        if toBreak
            break;
        end
    end    
    
    K = exp(-hyper * D);
end

