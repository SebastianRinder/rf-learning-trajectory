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
 %   trajectoriesPerPolicy = customOpts.trajectoriesPerPolicy;
    %actionSelectionFcn = customOpts.actionSelectionFcn;
    
%     Dtemp1 = 1;
%     Dtemp2 = 0;
%     Dtemp3 = 0;
    D = 0;

    for k = 1:size(knownTrajectory, 2)
        traj = knownTrajectory{1,k};

%         for t=1:length(traj.action)
%             state = traj.state(t,:);
%             action = traj.action(t,:);
%             [~,quadNew] = actionSelectionFcn(xNew, state, action);
%             %[~,quadj] = actionSelectionFcn(xj, state, action);
%             Dtemp1 = Dtemp1 * exp(quadNew - traj.probs(t));
%             Dtemp2 = Dtemp2 + quadNew - traj.probs(t);
%             Dtemp3 = Dtemp3 + traj.probs(t) - quadNew;
%         end
        
        [~,quadNew] = customOpts.actionSelectionFcn(xNew, traj.state, traj.action);
        Dtemp1 = prod(exp(quadNew - traj.prob));
        Dtemp2 = sum(quadNew - traj.prob);
        Dtemp3 = sum(traj.prob - quadNew);

        D = D + (Dtemp1 .* Dtemp2 + Dtemp3) ./ length(traj.action);
    end

    D = D / size(knownTrajectory, 2);
end

function D = monteCarloEst(xi, xj, traji, trajj, customOpts)
    %trajectoriesPerPolicy = customOpts.trajectoriesPerPolicy;
    %actionSelectionFcn = customOpts.actionSelectionFcn;
    
%    Dtemp = 0;
    D1 = 0;

    for k = 1:size(trajj, 2)
        traj = traji{1,k};

%         for t=1:length(traj.action)
%             state = traj.state(t,:);
%             action = traj.action(t,:);
%             %[~,quadi] = actionSelectionFcn(xi, state, action);
%             [~,quadj] = actionSelectionFcn(xj, state, action);
%             Dtemp = Dtemp + traj.prob(t) - quadj;
%         end

        [~,quadj] = customOpts.actionSelectionFcn(xj, traj.state, traj.action);
        D1 = D1 + sum(traj.prob - quadj) ./ length(traj.action);
    end
    
%    Dtemp = 0;
    D2 = 0;

    for k = 1:size(trajj, 2)
        traj = trajj{1,k};

%         for t=1:length(traj.action)
%             state = traj.state(t,:);
%             action = traj.action(t,:);
%             [~,quadi] = actionSelectionFcn(xi, state, action);
%             %[~,quadj] = actionSelectionFcn(xj, state, action);
%             Dtemp = Dtemp + traj.prob(t) - quadi;
%         end

        [~,quadi] = customOpts.actionSelectionFcn(xi, traj.state, traj.action);
        D2 = D2 + sum(traj.prob - quadi) ./ length(traj.action);
    end

    D = (D1 + D2) / size(trajj, 2);
end

