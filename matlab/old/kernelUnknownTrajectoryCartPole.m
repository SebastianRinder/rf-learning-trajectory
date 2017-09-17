function Ks = kernelUnknownTrajectoryCartPole(thetaNew, thetaj, trajj, hyper)
    tempD = [1;0;0];
    D = 0;
    toBreak = false;
    
    for k = 1:size(trajj,2)
        traj = trajj{1,k};
    
        for t=1:size(traj,1)
            state = traj{t,1}.state;
            action = traj{t,1}.action;

            [~, mu, sigma] = actionSelectionCartPole(thetaNew,state);
            normalize = normpdf(0, 0, sigma);
            Pnew = normpdf(action, mu, sigma) / normalize;
            
            if Pnew == 0
                D = Inf;
                toBreak = true;
            end

            [~, mu, ~] = actionSelectionCartPole(thetaj,state);
            Pj = normpdf(action, mu, sigma) / normalize;            
            
            tempD(1,1) = tempD(1,1) * Pnew/Pj;
            tempD(2,1) = tempD(2,1) + log(Pnew/Pj);
            tempD(3,1) = tempD(3,1) + log(Pj/Pnew);
        end
        if toBreak
            break;
        else
            D = D + tempD(1,1) * tempD(2,1) + tempD(3,1);
        end
    end

    Ks = exp(-hyper * D);
end