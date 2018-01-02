function D = trajectoryDistance(Xm, Xn, trajectories, isCovMat, opts)
    m = size(Xm,1);
    n = size(Xn,1);
    
    D = zeros(m,n);    
    
    if m == n && all(all(Xm == Xn)) %symmetric = true
        if ~isCovMat
            for i = 1:m
                for j = i:n
                    if all(Xm(i,:) == Xn(j,:))
                        D(i,j) = 0;
                    else
                        D(i,j) = monteCarloEst(Xm(i,:), Xn(j,:), trajectories(i,:), trajectories(j,:), opts);
                    end
                end
            end
        else
            states = [];
            for i = 1:size(trajectories,1)
                for j = 1:size(trajectories{i,1}) %only 1 policy per trajectory
                    states = [states; trajectories{i,1}.state];
                end
            end
            if size(states,1) > 500
                states = states(randsample(size(states,1),500),:);
            end
            for i = 1:m
                for j = i:n
                    if all(Xm(i,:) == Xn(j,:))
                        D(i,j) = 0;
                    else
                        D(i,j) = closedForm(Xm(i,:), Xn(j,:), states, opts);
                    end
                end
            end
        end
        D = D + D';
        
    else
        for i = 1:m
            for j = 1:n
                if all(Xm(i,:) == Xn(j,:))
                    D(i,j) = 0;
                else
                    D(i,j) = importanceSampling(Xm(i,:), trajectories(j,:), opts);
                end
            end
        end
    end
end

function D = importanceSampling(xNew, knownTrajectory, opts)
    D = 0;

    for k = 1:size(knownTrajectory, 2)
        traj = knownTrajectory{1,k};

        [~,probNew] = opts.actionSelectionFcn(xNew, traj.state, traj.action, opts);
        
        if isequal(opts.environment, 'cartPole')
            probDiff = sum(probNew - traj.prob);
            Dtemp1 = exp(probDiff);
            Dtemp2 = probDiff;
            Dtemp3 = -Dtemp2;
        else
            Dtemp1 = prod(probNew) / prod(traj.prob);
            Dtemp2 = sum(log(probNew ./ traj.prob));
%             Dtemp2 = log(Dtemp1);
            Dtemp3 = -Dtemp2;
        end

        D = D + (Dtemp1 .* Dtemp2 + Dtemp3) ./ length(traj.state);
    end

    D = D / size(knownTrajectory, 2);
end

function D = monteCarloEst(xi, xj, traji, trajj, opts)
    D1 = 0;
    for k = 1:size(traji, 2)
        traj = traji{1,k};

        [~,probj] = opts.actionSelectionFcn(xj, traj.state, traj.action, opts);
        if isequal(opts.environment, 'cartPole')
            D1 = D1 + sum(traj.prob - probj) ./ length(traj.state);
        else
            D1 = D1 + sum(log(traj.prob ./ probj)) ./ length(traj.state);
        end
    end
    
    D2 = 0;
    for k = 1:size(trajj, 2)
        traj = trajj{1,k};
        try
            [~,probi] = opts.actionSelectionFcn(xi, traj.state, traj.action, opts);
        catch me
            disp('fail');
        end
        if isequal(opts.environment, 'cartPole')
            D2 = D2 + sum(traj.prob - probi) ./ length(traj.state);
        else
            D2 = D2 + sum(log(traj.prob ./ probi)) ./ length(traj.state);
        end
    end

    D = (D1 + D2) / size(trajj, 2);
end

function D = closedForm(xi, xj, states, opts)
    if size(states,1) > 500
        states = states(randsample(size(states,1),500),:);
    end
    
    [~,~,mu1] = opts.actionSelectionFcn(xi, states, [], opts);
    [~,~,mu2] = opts.actionSelectionFcn(xj, states, [], opts);
    muDiff = mu1 - mu2;
    D = (muDiff' * muDiff) ./ size(states,1);
end