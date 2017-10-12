function [finalReward, traj] = execPolicyCartPole(policy, opts)
    traj.state = zeros(opts.timeSteps,5);
    traj.action = zeros(opts.timeSteps,1);
    traj.prob = zeros(opts.timeSteps,1);
    traj.cumReward = zeros(opts.timeSteps,1);
    
    cumReward = 0;    
    
    %position
    %velocity
    %acceleration
    %angle
    %angularVelocity
    angle = randn * 0.01;
    s = sign(angle);
    if s == 0
        s = 1;
    end
    opts.state0(4) = max(abs(angle),0.01) * s; %at least ~0.6 degree deviation from vertical
    state = opts.state0;    
    
    %policy = reshape(policy, [size(policy,2)/2, 2])';
    
    for i=1:opts.timeSteps
        [action, prob] = opts.actionSelectionFcn(policy,state,[]);
        [nextState, reward, finished] = simCartPole(state, action, opts.bounds);
        
        if i==opts.timeSteps
            if reward == 1
                reward = 1001;
            end
        end
        cumReward = cumReward + reward;
        traj.state(i,:) = state;
        traj.action(i) = action;
        traj.prob(i) = prob;
        traj.cumReward(i) = cumReward;
        
        if finished
            traj.state(i+1:end,:) = [];
            traj.action(i+1:end,:) = [];
            traj.prob(i+1:end,:) = [];
            traj.cumReward(i+1:end,:) = [];
            break;
        end
        state = nextState;
    end
    
    finalReward = cumReward / opts.timeSteps;
end