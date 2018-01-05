function [finalReward, traj] = objectiveFcn(policy, opts)
    traj.state = zeros(opts.timeSteps,size(opts.state0,2));
    traj.action = zeros(opts.timeSteps,1);
    traj.prob = zeros(opts.timeSteps,1);
    traj.cumReward = zeros(opts.timeSteps,1);
    
    cumReward = 0;    
    %probFac = 1;
    if isequal(opts.environment, 'cartPole')
        angle = randn * 0.01;
        s = sign(angle);
        if s == 0, s = 1; end
        opts.state0(4) = max(abs(angle),0.01) * s; %at least ~0.6 degree deviation from vertical
        
        %probFac = 1/sqrt(1e-4.*2.*pi);
    end
    state = opts.state0;    
    
    for i=1:opts.timeSteps
        [action, optProb] = opts.actionSelectionFcn(policy, state, [], opts.actionList);
        nextState = opts.simFcn(state, action, opts.bounds);        
        [reward, finished] = opts.rewardFcn(nextState, opts.bounds, i==opts.timeSteps);

        cumReward = cumReward + reward;
        traj.state(i,:) = state;
        traj.action(i,1) = action;
        traj.prob(i,1) = optProb;
        traj.cumReward(i,1) = cumReward;
        
        if finished
            traj.lastState = nextState;
            traj.state(i+1:end,:) = [];
            traj.action(i+1:end,:) = [];
            traj.prob(i+1:end,:) = [];
            traj.cumReward(i+1:end,:) = [];
            break;
        end
        state = nextState;
    end
    
    finalReward = cumReward; % / opts.timeSteps;
end