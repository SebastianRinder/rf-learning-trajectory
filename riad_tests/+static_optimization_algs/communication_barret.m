clc;
clear all;
close all;

barret = Environments.SL.barrett.BarrettCommunication();
parameters = rand(1,7);

for i = 1:10
    ballPosition = rand(1,3);
    ballVelocity = rand(1,3);
    ballVelocity = ballVelocity ./ norm(ballVelocity, 2);
    
    barret.SLSendTrajectory(parameters, 0, 1, 2, [ballPosition ballVelocity], 20);
    barret.SLSendTrajectory(parameters, 0, 2, 2, [], 20);
    [joints, jointsVel, jointsAcc, jointsDes, jointsDesVel, jointsDesAcc, torque, cart, SLstates, numCommand, stepIndex] = barret.SLGetEpisode(1);
end


