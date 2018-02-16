function D = trajectoryDistance(Xm, Xn, trajectories, isThompsonsSample, opts)
    m = size(Xm,1);
    n = size(Xn,1);
    
    D = zeros(m,n);    
    
    if m == n && all(all(Xm == Xn)) %symmetric = true
        if ~isThompsonsSample
            for i = 1:m
                for j = i+1:n
                    D(i,j) = monteCarloEst(Xm(i,:), Xn(j,:), trajectories(i,:), trajectories(j,:), opts);
                end
            end
        else
            states = [];
            for i = 1:size(trajectories,1) %only 1 sample per trajectory
                states = [states; trajectories{i,1}.state];
            end
            if size(states,1) > 500
                states = states(randsample(size(states,1),500),:);
            end
            for i = 1:m
                for j = i+1:n
                    D(i,j) = closedForm(Xm(i,:), Xn(j,:), states, opts);
                end
            end
        end
        D = D + D';

    else        
        for i = 1:m
            for j = 1:n
                D(i,j) = importanceSampling(Xm(i,:), trajectories(j,:), opts);
            end
        end
    end
    
    D(D<0) = 0;
end

function D = importanceSampling(xNew, knownTrajectory, opts)
    D = 0;

    for k = 1:size(knownTrajectory, 2)
        traj = knownTrajectory{1,k};

        [~,probNew] = opts.actionSelectionFcn(xNew, traj.state, traj.action, opts.actionMisc);
        
        if isequal(opts.actionSpace, 'continuous')
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

        [~,probj] = opts.actionSelectionFcn(xj, traj.state, traj.action, opts.actionMisc);
        if isequal(opts.actionSpace, 'continuous')
            D1 = D1 + sum(traj.prob - probj) ./ length(traj.state);
        else
            D1 = D1 + sum(log(traj.prob ./ probj)) ./ length(traj.state);
        end
    end
    
    D2 = 0;
    for k = 1:size(trajj, 2)
        traj = trajj{1,k};
        
        [~,probi] = opts.actionSelectionFcn(xi, traj.state, traj.action, opts.actionMisc);
        if isequal(opts.actionSpace, 'continuous')
            D2 = D2 + sum(traj.prob - probi) ./ length(traj.state);
        else
            D2 = D2 + sum(log(traj.prob ./ probi)) ./ length(traj.state);
        end
    end

    D = (D1 + D2) / size(trajj, 2);
end

function D = closedForm(xi, xj, states, opts)
    if isequal(opts.actionSpace, 'continuous')
        [~,~,mu1] = opts.actionSelectionFcn(xi, states, [], opts.actionMisc);
        [~,~,mu2] = opts.actionSelectionFcn(xj, states, [], opts.actionMisc);
        muDiff = mu1 - mu2;
        D = (muDiff' * muDiff) ./ size(states,1);
    else
        [~,prob1,~] = opts.actionSelectionFcn(xi, states, [], opts.actionMisc);
        [~,prob2,~] = opts.actionSelectionFcn(xj, states, [], opts.actionMisc);
        Dtemp = log(prob1./prob2);
        D = sum(prob1.* Dtemp) + sum(-prob2.* Dtemp);
        D = D ./ size(states,1);
    end
end