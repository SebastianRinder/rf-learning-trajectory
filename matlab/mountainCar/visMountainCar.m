function visMountainCar(traj, bounds)
    skipVisValue = selectFigure('Visualize Mountain Car');
    clf;
    hold on;
    
    x = linspace(bounds.position(1,1),bounds.position(1,2));
    % Draw y = sin(3*p) to representate a mountain
    plot(x, sin(3 .* x));
    
    checkBoxSkip = uicontrol('Style', 'checkbox', 'String', 'Skip Vis',...
        'Position', [90 0 90 20],...
        'Value', skipVisValue);

    % Text labels to show current input and simulation time.
    lblTime = uicontrol('style','text');
    lblAction = uicontrol('style','text');
    lblReward = uicontrol('style','text');
    set(lblTime,'Position', [1 20 70 20]);
    set(lblAction,'Position', [1 50 70 20]);
    set(lblReward,'Position', [1 80 70 20]);

    car = plot(0,0, 'or', 'LineWidth', 4);
    
    T = length(traj{1,1}.action);
    
    for i = 1:T+1
        if get(checkBoxSkip, 'Value') == 1
            break;
        end
        
        if i == T+1
            p = traj.lastState(1,1);
        else
            p = traj.state(i,1);
        end
        
        set(car, 'XData', p);
        set(car, 'YData', sin(3 * p)); 

        if i < T+1
            set(lblTime,'String', ['t: ',num2str(i - 1)]);
            set(lblAction,'String',  ['a: ',num2str(traj.action(i))]);
            set(lblReward,'String', ['cr: ',num2str(traj.cumReward(i))]);
        end
        
        %drawnow
        pause(0.02);
    end
    
    hold off;
end