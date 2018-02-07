function [opts] = environmentSettings(env, platform, visualize)
    
    if isequal(platform, 'pygym')

        P = py.sys.path;
        if count(P,'pygym') == 0
            insert(P,int32(0),'pygym');
        end   

        if isequal(env, 'cartPole')
            pystr = 'CartPole-v0';        
            opts.dim = 10; %4 state values + 1 bias value per action
            opts.timeSteps = 200;
            opts.actionSelectionFcn = @actionSelectionDiscrete;
            mod = py.importlib.import_module('gymCartPole');
            py.reload(mod);            
            ret = py.gymCartPole.init_environment(pystr);
            
        elseif isequal(env, 'mountainCar')
            pystr = 'MountainCar-v0';
            opts.dim = 30;
            opts.timeSteps = 200;
            opts.actionSelectionFcn = @actionSelectionMountainCarDiscrete;
            mod = py.importlib.import_module('gymMountainCar');
            py.reload(mod);            
            ret = py.gymMountainCar.init_environment(pystr);
           
         elseif isequal(env, 'acroBot')
            pystr = 'Acrobot-v1';
            opts.dim = 21;
            opts.timeSteps = 500;
            opts.actionSelectionFcn = @actionSelectionDiscrete;
            mod = py.importlib.import_module('gymAcroBot');
            py.reload(mod);            
            ret = py.gymAcroBot.init_environment(pystr);
            
        elseif isequal(env, 'mountainCarContinuous')
            pystr = 'MountainCarContinuous-v0';
            opts.dim = 10;
            opts.timeSteps = 999;
            opts.actionSelectionFcn = @actionSelectionMountainCarContinuous;
            mod = py.importlib.import_module('gymMountainCarContinuous');
            py.reload(mod);            
            ret = py.gymMountainCarContinuous.init_environment(pystr);
        end
        
        pyObservationSpace = ret{1,2};
        pyObservationSpace = cell(pyObservationSpace.shape);
        opts.stateDim = double(pyObservationSpace{1,1});
        
        pyActionSpace = ret{1,1};
        if isa(pyActionSpace, 'py.gym.spaces.discrete.Discrete')            
            opts.actionDim = double(pyActionSpace.n);
            opts.actionSpace = 'discrete';
            opts.actionMisc = 0:double(pyActionSpace.n)-1;
        else
            opts.actionSpace = 'continuous';
            opts.actionMisc = 1e-3;
        end
    end

    if isequal(platform, 'matlab')
        if isequal(env, 'cartPole')
            opts.dim = 4; %5
            opts.timeSteps = 200;
            opts.actionSpace = 'continuous';
            opts.bounds.position = [-2.4, 2.4];
            opts.bounds.angle = [-12 * pi / 180, 12 * pi / 180];
            
            opts.actionSelectionFcn = @actionSelectionContinuous;
            opts.simFcn = @simCartPole;
            opts.rewardFcn = @rewardCartPole;
            opts.visFcn = @visCartPole;
            
            opts.actionMisc = 1e-3;
%         elseif isequal(env, 'bipedalWalker')
%             opts.dim = 100; %24 state values + 1 bias value per action
%             opts.actionSpace = 'continuous';
%             opts.actionList = [0,0,0,0];
%             opts.actionSelectionFcn = @actionSelectionBipedalWalker;
% 
%             mod = py.importlib.import_module('gymBipedalWalker');
%             py.reload(mod);
%             py.gymBipedalWalker.init_environment('BipedalWalker-v2');
%         elseif isequal(env, 'mountainCar')
%             opts.dim = 9*2;
%             opts.actionList = [-1,1];   %apply acceleration to the rear or forward
%             opts.bounds.position = [-1.2, 0.5];
%             opts.bounds.velocity = [-0.07, 0.07];
%             opts.state0 = [-0.5, 0]; %position, velocity
%             opts.timeSteps = 400;
% 
%             opts.actionSelectionFcn = @actionSelectionMountainCar;
%             opts.simFcn = @simMountainCar;
%             opts.rewardFcn = @rewardMountainCar;
%             opts.visFcn = @visMountainCar;
%         elseif isequal(env, 'acroBot')
%             opts.dim = 4*3;
%             opts.actionList = [-1,0,1];       %apply torque to hip
%             opts.bounds.angle1 = [-pi, pi];
%             opts.bounds.angle2 = [-pi, pi];
%             opts.bounds.velocity1 = [-4*pi, 4*pi];
%             opts.bounds.velocity2 = [-9*pi, 9*pi];
%             opts.state0 = [0, 0, 0, 0]; %angle1, angle2, angularVelocity1, angularVelocity2
%             opts.timeSteps = 400;
% 
%             opts.actionSelectionFcn = @actionSelectionAcroBot;
%             opts.simFcn = @simAcroBot;
%             opts.rewardFcn = @rewardAcroBot;
%             opts.visFcn = @visAcroBot;        
%         elseif isequal(env, 'randGauss')
%             %TODO
        else
            error(['no environment for ', env]);
        end
    end
    
    opts.environment = env;
    opts.visualize = visualize;    
    %addpath(env);
end