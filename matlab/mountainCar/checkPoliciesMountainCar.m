function [eta,theta,xxx] = checkPoliciesMountainCar()
    addpath('..');

    samples = 1000;
    episodes = 1000;
    
    worldBounds.position = [-1.2, 0.5];
    worldBounds.velocity = [-0.07, 0.07];
    
    dim = 18;
        
    ub = ones(1,dim);
    lb = -ub;
    
    state0.position = -0.5;
    state0.velocity = 0;
    

    theta = zeros(samples,dim);
    eta = zeros(samples,1);

    maxEta = 0;
    for i=1:samples
        theta(i,:) = randTheta(lb,ub);
        [eta(i,1), traj] = execPolicyMountainCar(theta(i,:), state0, episodes, worldBounds);
        if eta(i,1) > maxEta
            maxEta = eta(i,1);
            %visMountainCar(traj, worldBounds);
        end
        if mod(i,floor(samples/10)) == 0
            disp([int2str(i*100/samples),'%']);
        end
    end
    
    figure;
    histogram(eta);
    figure;
    histogram(eta(eta~=-1000));
    
    x = [eta,theta];
    xx = x(eta>400,:);
    xxx = [min(xx);max(xx)];
end