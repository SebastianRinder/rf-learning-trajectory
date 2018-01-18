function plotBallTrajectory(threeDTraj)
plot3(threeDTraj(:,1), threeDTraj(:,2), threeDTraj(:,3), 'b');
hold on;
traj = animatedline(threeDTraj(1,1), threeDTraj(1,2), threeDTraj(1,3), 'Color','r');
labelEveryXData = 50;

for k = 2:size(threeDTraj, 1)
    % first line
    addpoints(traj, threeDTraj(k,1), threeDTraj(k,2), threeDTraj(k,3));
    if(mod(k, labelEveryXData) == 0)
        text(threeDTraj(k,1), threeDTraj(k,2), threeDTraj(k,3), int2str(k));
    end
    % update screen
    drawnow %limitrate
end
hold off;
end

