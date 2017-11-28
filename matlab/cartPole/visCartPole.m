function visCartPole(traj, bounds)
    skipVisValue = selectFigure('Visualize Cart Pole');
    clf;
    
    axis equal;
    axis([min(bounds.position), max(bounds.position) 0 5]);
    hold on;
    
    cart = rectangle;
    pole = line;
    l = 3;  %length of visualized pole
     
    checkBoxSkip = uicontrol('Style', 'checkbox', 'String', 'Skip Vis',...
        'Position', [360 60 90 20],...
        'Value', skipVisValue);
    
    lblTime = uicontrol('style','text');
    lblAction = uicontrol('style','text');
    lblReward = uicontrol('style','text');
    lblPosition = uicontrol('style','text');
    lblVelocity = uicontrol('style','text');
    lblAcceleration = uicontrol('style','text');
    lblAngle = uicontrol('style','text');
    lblAngularVelocity = uicontrol('style','text');
    
    set(lblTime,'Position', [0 60 90 20]);
    set(lblAction,'Position', [0 30 90 20]);
    set(lblReward,'Position', [0 0 90 20]);    
    set(lblPosition,'Position', [90 60 90 20]);
    set(lblVelocity,'Position', [90 30 90 20]);
    set(lblAcceleration,'Position', [90 0 90 20]);    
    set(lblAngle,'Position', [180 60 90 20]);
    set(lblAngularVelocity,'Position', [180 30 90 20]);
    
    T = length(traj.action);
    title(['time steps: ', num2str(T), '      final reward: ', num2str(traj.cumReward(end))]);
    
    for i = 1:T+1
        if get(checkBoxSkip, 'Value') == 1
            break;
        end
        if i == T+1
            x = traj.lastState(1,1);
            omega = traj.lastState(1,4);
        else
            x = traj.state(i,1);
            omega = traj.state(i,4);
        end
        
        pxp=[x x+l*sin(omega)];
        pyp=[1.25 1.25+l*cos(omega)];

        cart.Position = [x-1, 0.25, 2, 1];
        pole.XData = pxp;
        pole.YData = pyp;
        
        if i < T+1
            set(lblTime,'String', ['t: ',num2str(i)]);
            set(lblAction,'String',  ['a: ',num2str(traj.action(i))]);
            set(lblReward,'String', ['cr: ',num2str(traj.cumReward(i))]);

            set(lblPosition,'String', ['p: ',num2str(traj.state(i,1))]);
            set(lblVelocity,'String',  ['v: ',num2str(traj.state(i,2))]);
            set(lblAcceleration,'String', ['acc: ',num2str(traj.state(i,3))]);

            set(lblAngle,'String',  ['ang: ',num2str(traj.state(i,4))]);
            set(lblAngularVelocity,'String', ['angV: ',num2str(traj.state(i,5))]);
        end
        
        pause(0.02);
    end
    hold off;
end
