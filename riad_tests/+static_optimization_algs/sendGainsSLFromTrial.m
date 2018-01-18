% 	goto_pos       = 0,
% 	set_cannon_trg = 1,
% 	set_cannon_pos = 2,
% 	get_ball_state = 3,
% 	get_IK_hitPos  = 4,
% 	run_sim        = 5,
% 	set_ball_state = 6,
% 	execute_gains  = 7

% send gains to SL simulator
clc;
clear variables;
close all;

load('/home/akrour/postdoc/code/policysearchtoolbox/+Experiments/data/TableTennis_TimeDep2/TableTennis_TimeDependentREPS/settings001/eval001/trial009/trial.mat')
minHitPosition = [trial.contextSampler.MinBallPos trial.contextSampler.DefaultBallPos];
maxHitPosition = [trial.contextSampler.MaxBallPos trial.contextSampler.DefaultBallPos];
% minHitPosition = [-.60 .85 .15];
% maxHitPosition = [.30 .85 .15];

rangeHitPosition = maxHitPosition - minHitPosition;

biorob = Environments.SL.bioroblinaxis.BioroblinaxisCommunication();

n_dof = biorob.dimJoints;
waitTime    = 0;                % Waiting time between commands/states
maxSteps    = 2;                % Number of commands/states
trajLength = trial.GainProvider.maxTimeStep;
effectiveTrajLength = trial.GainProvider.trajLength;
effectiveStateDim = n_dof * 2 + 3;
dimControl = n_dof;
startSupraStep = 1;
startGains = ceil(effectiveTrajLength / n_dof) + 1;
usedStates = [ones(1, effectiveStateDim) zeros(1, 3)]; %we don't use ball velocity
controlledJoints = ones(1, n_dof);
firstTimeStep = trial.GainProvider.firstTimeStep;
reuseGainMode = trial.GainProvider.reuseGainMode;
additionalWaitingTime = 600;

applyNoise = 0;
%%

initPos = [-.15 .15 -.5 1 0 -1 -.8 0 -1.2];
nbEval = 150;
times = zeros(1, nbEval);
rewards = zeros(1, nbEval);
for i = 1:nbEval
    i
    tic 
    %go_to pos
    stateBuffer = [0 initPos];
    biorob.SLSendTrajectory(initPos, waitTime, 1, 1, stateBuffer);
    %set_cannon_target
    hitPos = rand(1, 3) .* rangeHitPosition + minHitPosition
    %%% hitPos = proTrials.trials{1}.cannonTarget;
    stateBuffer = [1 hitPos]; % Ballcannon serve target
    [reward, slbuff, flag] = biorob.SLSendTrajectory(initPos, waitTime, 1, maxSteps, stateBuffer);
    %execute_gains  
    %stateBuffer = [7 effectiveStateDim dimControl firstTimeStep trajLength startSupraStep startGains usedStates controlledJoints reuseGainMode applyNoise additionalWaitingTime];
    %%% first, how the supra step are structured
    trial.GainProvider.applyNoise = applyNoise;
    stateBuffer = trial.SLenvironment.getStateBufferForSL;
    traj = trial.SLenvironment.getGains;
    [reward, slbuff, flag] = biorob.SLSendTrajectory(traj, waitTime, 2, maxSteps, stateBuffer);
    %trial.SLenvironment.isInit = true;
%    [states, actions] = trial.SLenvironment.sampleEpisode();
    %%% checking if gains correctly applied
    [joints, jointsVel, jointsAcc, jointsDes, jointsDesVel, jointsDesAcc, torque, cart, SLstates, numCommand, stepIndex, torqueDes] = biorob.SLGetEpisode(trajLength + additionalWaitingTime, 13);
    stepClustering = trial.GainProvider.getStepClustering();
    sl_error_inGains = zeros(effectiveTrajLength, n_dof);
    for t = 0:(effectiveTrajLength - 1)
        K = trial.GainProvider.getGainsForTimeStep(stepClustering(t+1));
        sl_error_inGains(t + 1, :) = [1 joints(t + firstTimeStep, :) jointsVel(t + firstTimeStep, :)...
            SLstates(t + firstTimeStep, 1:3)] * K - torqueDes(t + firstTimeStep, :);
    end
    [~, rewardsSteps] = trial.SLenvironment.rewardFunc.getTrajectoryReward(SLstates, ...
                firstTimeStep, trial.GainProvider.stepLength, ...
                trial.GainProvider.numTimeSteps);
    reward = max(rewardsSteps)
    time = toc;
    times(i) = time;
    rewards(i) = reward;
end
plot (1:nbEval, times, '*b');
hold on;
plot (1:nbEval, rewards, '.r');
hold off;
