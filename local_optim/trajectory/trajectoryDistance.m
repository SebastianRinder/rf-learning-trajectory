function D = trajectoryDistance(xi, xj, traji, trajj, opts)    
    
    if isempty(traji) && ~isempty(trajj)
        D = importanceSampling(xi, trajj, opts);
        
	elseif ~isempty(traji) && ~isempty(trajj)
        D = monteCarloEst(xi, xj, traji, trajj, opts);
        
    elseif ~isempty(traji) && isempty(trajj)
        D = closedForm(xi, xj, traji, opts);
        
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

function D = closedForm(xi, xj, trajectory, opts)
    D = 0;
    l = 0;
    for k = size(trajectory, 1):-1:1
        traj = trajectory{k,1};
        if l + length(traj.state) > 500
            traj.state = traj.state(end-500+l:end,:);
        end
        [~,~,mu1] = opts.actionSelectionFcn(xi, traj.state, [], opts);
        [~,~,mu2] = opts.actionSelectionFcn(xj, traj.state, [], opts);
        muDiff = mu1 - mu2;
        D = D + (muDiff' * muDiff) ./ length(traj.state);
        l = l + length(traj.state);
        if l > 500, break; end
    end
    
    %D = D / (opts.errorVariance ^ 2);
    
    if ~isequal(opts.environment, 'cartPole')       
        %TODO
        keyboard;
    end
end
