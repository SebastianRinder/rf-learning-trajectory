%play arm and ball trajectory
load('/home/akrour/postdoc/code/policysearchtoolbox/+Experiments/data/test/TableTennis_TrajectoryBased_SimpleAlg/Table_Tennis_201507242016_01/eval001/trial001/iter_321491_321500.mat');
savedTraj = iter321500.data.steps;

biorob   = Environments.SL.bioroblinaxis.BioroblinaxisCommunication(false);

initPos  = [-.15 .15 -.5 1 0 -1 -.8 0 -1.2];

waitTime = 0;                % Waiting time between commands/states
maxSteps = 2;                % Number of commands/states

for i=1:10000
stateBuffer = [0 initPos]; % go to init pos
biorob.SLSendTrajectory(zeros(1, biorob.dimJoints), waitTime, 1, maxSteps, stateBuffer);

% building the trajectory
traj = savedTraj.states(:, 1:2:(2*biorob.dimJoints));
traj = [traj; [savedTraj.SLstates(:, 1:3) zeros(size(savedTraj.SLstates, 1), biorob.dimJoints-3)]];

stateBuffer(1) = 5; % run_sim
biorob.SLSendTrajectory(traj, waitTime, 2, maxSteps, stateBuffer);
end