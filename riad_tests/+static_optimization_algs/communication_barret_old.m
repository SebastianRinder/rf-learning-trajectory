clc;
clear all;
close all;

barret = Environments.SL.barrett.BarrettCommunication();
parameters = rand(50000,7);

time = 0;
maxCommand = 2;
timeOut = 20;
stepLength = 10;

for i = 1:10
    ballPosition = rand(1,3);
    ballVelocity = rand(1,3);
    ballVelocity = ballVelocity ./ norm(ballVelocity, 2);
    
    barret.SLSendTrajectory(parameters, time, 1, maxCommand, [ballPosition ballVelocity stepLength], timeOut);
    [slReturn, slReturnState] = barret.SLSendTrajectory(parameters, time, 2, maxCommand, [], timeOut);
    
    [joints, jointsVel, jointsAcc, jointsDes, jointsDesVel, jointsDesAcc, torque, cart, SLstates, numCommand, stepIndex] = barret.SLGetEpisode(1);
end

