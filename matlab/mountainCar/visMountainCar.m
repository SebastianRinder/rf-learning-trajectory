function visMountainCar(traj, bounds)
    close all;
    figure;
    hold on;
    
    x = linspace(bounds.position(1,1),bounds.position(1,2));
    % Draw y = sin(3*p) to representate a mountain
    plot(x, sin(3 .* x));

    % Text labels to show current input and simulation time.
    lblTime = uicontrol('style','text');
    lblAction = uicontrol('style','text');
    lblReward = uicontrol('style','text');
    set(lblTime,'Position', [1 20 70 20]);
    set(lblAction,'Position', [1 50 70 20]);
    set(lblReward,'Position', [1 80 70 20]);

    car = plot(0,0, 'or', 'LineWidth', 4);
    
    for t = 1:size(traj,1)
        p = traj{t, 1}.state.position(1,1);
        set(car, 'XData', p);
        set(car, 'YData', sin(3 * p)); 

        set(lblTime,'String', ['t: ',num2str(t - 1)]);
        set(lblAction,'String',  ['a: ',num2str(traj{t, 1}.action)]);
        set(lblReward,'String', ['cr: ',num2str(traj{t, 1}.cumReward)]);

        drawnow
        pause(0.02);
    end
end