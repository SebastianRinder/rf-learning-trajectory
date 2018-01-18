%play arm and ball trajectory
biorob   = Environments.SL.bioroblinaxis.BioroblinaxisCommunication(false);

%initPos  = [-.15 .15 -.5 1 0 -1 -.8 0 -1.2];

waitTime = 0;                % Waiting time between commands/states
maxSteps = 2;                % Number of commands/states

cannonTarget = [-0.8 4 2];   
stateBuffer = [1 cannonTarget]; % set cannon target
biorob.SLSendTrajectory(zeros(1, biorob.dimJoints), waitTime, 1, maxSteps, stateBuffer);

goalPos = [0 -3 -0.8];
stateBuffer = [4 goalPos]; % set cannon target
biorob.SLSendTrajectory(traj, waitTime, 2, maxSteps, stateBuffer);
