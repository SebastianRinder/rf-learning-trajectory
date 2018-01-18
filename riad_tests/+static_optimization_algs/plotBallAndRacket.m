function plotBallAndRacket(ballTraj, racketTraj, labelEveryXData)
plot3(ballTraj(:,1), ballTraj(:,2), ballTraj(:,3), 'b');
hold on;
plot3(racketTraj(:,1), racketTraj(:,2), racketTraj(:,3), 'm');
pause
traj1 = animatedline(ballTraj(1,1), ballTraj(1,2), ballTraj(1,3), 'Color','r');
traj2 = animatedline(racketTraj(1,1), racketTraj(1,2), racketTraj(1,3), 'Color','g');

if(~exist('labelEveryXData', 'var'))
    labelEveryXData = 100;
end

for k = 2:size(ballTraj, 1)
    % first line
    addpoints(traj1, ballTraj(k,1), ballTraj(k,2), ballTraj(k,3));
    addpoints(traj2, racketTraj(k,1), racketTraj(k,2), racketTraj(k,3));

    if(mod(k, labelEveryXData) == 0)
        text(ballTraj(k,1), ballTraj(k,2), ballTraj(k,3), int2str(k));
        text(racketTraj(k,1), racketTraj(k,2), racketTraj(k,3), int2str(k));
    end
    % update screen
    drawnow %limitrate
end
hold off;
end
