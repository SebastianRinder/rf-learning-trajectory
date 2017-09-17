function [eta,theta,ub,lb] = checkPoliciesCartPole()
    addpath('..');

    samples = 10000;
    executionTimeSteps = 1000;
    dim = 4;    
    
    worldBounds.position = [-10, 10];
    worldBounds.angle = [-pi/2, pi/2];
    worldBounds.rewardPosition = [-0.5, 0.5];
    worldBounds.rewardAngle = [-pi/9, pi/9];
    worldBounds.action = [-1,1];
    
    state0.position = 0;
    state0.velocity = 0;
    state0.acceleration = 0;
    state0.angle = 0;
    state0.angleVelocity = 0;

%      ub = [7.3 0.11 9.8 9.5];
%      lb = [-0.2 -0.087 -2.6 -0.0082];
     
    ub = ones(1, dim);
    lb = -ub;

    theta = zeros(samples,dim);
    eta = zeros(samples,1);
    %eps = zeros(steps,1);

    %theta(1,:) = randTheta(lb,ub);
    %theta = zeros(1,4);

    for i=1:samples
        %ub = randn(1,dim) * 10;
        %lb = randn(1,dim) * -10;
        theta(i,:) = randTheta(lb,ub);
        [eta(i,1), ~] = execPolicyCartPole(theta(i,:), state0, executionTimeSteps, worldBounds);
        if mod(i,floor(samples/10)) == 0
            disp([int2str(i*100/samples),'%']);
        end
    end
    
    histogram(eta);
end