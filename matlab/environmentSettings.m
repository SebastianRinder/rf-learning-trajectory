function [opts] = environmentSettings(env, visualize)
    if isequal(env, 'cartPole')
        opts.dim = 4;
        opts.actionList = [];	%continuous action selection
        opts.bounds.position = [-5, 5];
        opts.bounds.angle = [-90 * pi / 180, 90 * pi / 180];
        opts.bounds.rewardPosition = [-1, 1];
        opts.bounds.rewardAngle = [-12 * pi / 180, 12 * pi / 180];
        opts.state0 = zeros(1,5); %position, velocity, acceleration, angle, angularVelocity
        opts.timeSteps = 1000;
    
        opts.actionSelectionFcn = @actionSelectionCartPole;
        opts.simFcn = @simCartPole;
        opts.rewardFcn = @rewardCartPole;
        opts.visFcn = @visCartPole;
    elseif isequal(env, 'mountainCar')
        opts.dim = 9*2;
        opts.actionList = [-1,1];   %apply acceleration to the rear or forward
        opts.bounds.position = [-1.2, 0.5];
        opts.bounds.velocity = [-0.07, 0.07];
        opts.state0 = [-0.5, 0]; %position, velocity
        opts.timeSteps = 400;
                
        opts.actionSelectionFcn = @actionSelectionMountainCar;
        opts.simFcn = @simMountainCar;
        opts.rewardFcn = @rewardMountainCar;
        opts.visFcn = @visMountainCar;
    elseif isequal(env, 'acroBot')
        opts.dim = 4*3;
        opts.actionList = [-1,0,1];       %apply torque to hip
        opts.bounds.angle1 = [-pi, pi];
        opts.bounds.angle2 = [-pi, pi];
        opts.bounds.velocity1 = [-4*pi, 4*pi];
        opts.bounds.velocity2 = [-9*pi, 9*pi];
        opts.state0 = [0, 0, 0, 0]; %angle1, angle2, angularVelocity1, angularVelocity2
        opts.timeSteps = 400;
                
        opts.actionSelectionFcn = @actionSelectionAcroBot;
        opts.simFcn = @simAcroBot;
        opts.rewardFcn = @rewardAcroBot;
        opts.visFcn = @visAcroBot;        
    end
    
    opts.environment = env;
    opts.visualize = visualize;    
    addpath(env);
end