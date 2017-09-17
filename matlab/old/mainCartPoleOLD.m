function ret = mainCartPoleOLD()
    addpath('..');
    
    steps = 100;
    episodes = 1000;
    dim = 4;
    kernel = @kernelSqExp;
    
    bounds = [      -10,         10;    %position bounds
                  -pi/2,	   pi/2;    %angular bounds
                   -0.5,        0.5;    %position bounds for positive reward
                 -pi/18,      pi/18];   %angular bounds for positive reward
    actions = 0;                        %continuous action mapping
    state0 = [0, 0, 0, 10e-4, 0];       %position, velocity, acceleration, angle, angular velocity
    
    ub = [7.3 0.11 9.8 9.5];            %policy boundaries from random sampling
    lb = [-0.2 -0.087 -2.6 -0.0082];
    
    theta = zeros(steps,dim);
    eta = zeros(steps,1);
    eps = zeros(steps,1);
    
    theta(1,:) = randTheta(lb,ub);
    [eta(1,1), ~, eps(1,1)] = execPolicy(theta(1,:), state0, actions, episodes, bounds, @simCartPole, @stateFeatureCartPole);
    
    for i=2:steps
        K = kernel(theta(1:i-1,:), theta(1:i-1,:), 1);
        fun = @(x) -acquisition(x, eta(1:i-1,1), theta(1:i-1,:), K, kernel);        
%         nextTheta = fmincon(fun, randTheta(lb,ub),[],[],[],[],lb,ub);
        nextTheta = globalMin(fun,lb,ub);
        theta(i,:) = nextTheta;
        
        [nextEta, traj, datEps] = execPolicy(nextTheta, state0, actions, episodes, bounds, @simCartPole, @stateFeatureCartPole);
        eta(i,1) = nextEta;
        eps(i,1) = datEps;
        
         if datEps == 1000 && nextEta < 3000
             visCartPole(traj, i, bounds);
         end
        i
    end
    ret = [eta,eps,theta];
end