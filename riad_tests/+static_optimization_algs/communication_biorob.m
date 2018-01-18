clc;
clear all;
close all;

n_dof = Environments.SL.bioroblinaxis.BioroblinaxisCommunication.n_dof;
waitTime    = 0;                % Waiting time between commands/states
maxSteps    = 2;                % Number of commands/states
stateBuffer = [8 -0.3 0.85 0.15]; % Ballcannon serve target
%%

% WARNING: At execution, the state machine initial state is 0, so you
% must initially then next state, state 1. After the initial
% start, you must call SLSendTrajectory with state 0, then 1, 2...
%or you must use timeOut to add delay between state transitions

biorob = Environments.SL.bioroblinaxis.BioroblinaxisCommunication();

initPos = [-.8 0 -.25 .6 0 .4 -.4 0 0];


for i = 1:100000
    traj = rand(30,n_dof);
    biorob.SLSendTrajectory(traj, waitTime, 1, maxSteps, stateBuffer);
    stateBuffer = [8 -0.3 0.85 0.15]; % Ballcannon serve target
    [reward, slbuff, flag] = biorob.SLSendTrajectory(initPos, waitTime, 2, maxSteps, stateBuffer);
    [joints, jointsVel, jointsAcc, jointsDes, jointsDesVel, jointsDesAcc, torque, cart, SLstates, numCommand, stepIndex] = biorob.SLGetEpisode(1);
end


