function D = trajectoryDistance(xi, xj, traji, trajj, customOpts)    
    
    if isempty(traji) && ~isempty(trajj)
        D = importanceSampling(xi, trajj, customOpts);
        
	elseif ~isempty(traji) && ~isempty(trajj)
        D = monteCarloEst(xi, xj, traji, trajj, customOpts);
        
    else
        keyboard;
    end
end

function D = importanceSampling(xNew, knownTrajectory, customOpts)
    D = 0;

    for k = 1:size(knownTrajectory, 2)
        traj = knownTrajectory{1,k};

        [~,probNew] = customOpts.actionSelectionFcn(xNew, traj.state, traj.action);
        
        if isequal(customOpts.problem, 'cartPole')
            probDiff = sum(probNew - traj.prob);
            Dtemp1 = exp(probDiff);
            Dtemp2 = probDiff;
            Dtemp3 = -probDiff;
        else
            Dtemp1 = prod(probNew) / prod(traj.prob);
            Dtemp2 = log(Dtemp1);
            Dtemp3 = -Dtemp2;
        end

        D = D + (Dtemp1 .* Dtemp2 + Dtemp3) ./ length(traj.action);
    end

    D = D / size(knownTrajectory, 2);
end

function D = monteCarloEst(xi, xj, traji, trajj, customOpts)
    D1 = 0;
    for k = 1:size(trajj, 2)
        traj = traji{1,k};

        [~,probj] = customOpts.actionSelectionFcn(xj, traj.state, traj.action);
        if isequal(customOpts.problem, 'cartPole')
            D1 = D1 + sum(traj.prob - probj) ./ length(traj.action);
        else
            D1 = D1 + log(prod(traj.prob) / prod(probj)) ./ length(traj.action);
        end
    end
    
    D2 = 0;
    for k = 1:size(trajj, 2)
        traj = trajj{1,k};

        [~,probi] = customOpts.actionSelectionFcn(xi, traj.state, traj.action);
        if isequal(customOpts.problem, 'cartPole')
            D2 = D2 + sum(traj.prob - probi) ./ length(traj.action);
        else
            D2 = D2 + log(prod(traj.prob) / prod(probi)) ./ length(traj.action);
        end
    end

    D = (D1 + D2) / size(trajj, 2);
end

