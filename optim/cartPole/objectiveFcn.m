function [finalReward, trajectory] = objectiveFcn(policy, opts)
    trajectory.state = zeros(opts.timeSteps,size(opts.state0,2));
    trajectory.action = zeros(opts.timeSteps,1);
    trajectory.prob = zeros(opts.timeSteps,1);
    trajectory.cumReward = zeros(opts.timeSteps,1);
    
    cumReward = 0;    
    %probFac = 1;
    if isequal(opts.environment, 'cartPole')
        opts.state0 = rand(1,5) * 0.1 - 0.05;
%         angle = randn * 0.01;
%         s = sign(angle);
%         if s == 0, s = 1; end
%         opts.state0(4) = max(abs(angle),0.01) * s; %at least ~0.6 degree deviation from vertical
        
        %probFac = 1/sqrt(1e-4.*2.*pi);
    end
    state = opts.state0;    
    
    for i=1:opts.timeSteps
        [action, optProb] = opts.actionSelectionFcn(policy, state, [], opts);
        nextState = opts.simFcn(state, action, opts.bounds);        
        [reward, finished] = opts.rewardFcn(nextState, opts.bounds, i==opts.timeSteps);

        cumReward = cumReward + reward;
        trajectory.state(i,:) = state;
        trajectory.action(i,1) = action;
        trajectory.prob(i,1) = optProb;
        trajectory.cumReward(i,1) = cumReward;
        
        if finished
            trajectory.lastState = nextState;
            trajectory.state(i+1:end,:) = [];
            trajectory.action(i+1:end,:) = [];
            trajectory.prob(i+1:end,:) = [];
            trajectory.cumReward(i+1:end,:) = [];
            break;
        end
        state = nextState;
    end
    
    finalReward = cumReward; %/ opts.timeSteps;
end