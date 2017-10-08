function visCartPole(traj, bounds)
    findFigure('Visualize Cart Pole');
    clf;
    
    axis equal;
    axis([bounds.position(1,1), bounds.position(1,2) 0 6]);
    hold on;
    
    cart = rectangle;
    pole = line;
    l = 3;  %length of visualized pole
    
    lblTime = uicontrol('style','text');
    lblAction = uicontrol('style','text');
    lblReward = uicontrol('style','text');
    lblPosition = uicontrol('style','text');
    lblVelocity = uicontrol('style','text');
    lblAcceleration = uicontrol('style','text');
    lblAngle = uicontrol('style','text');
    lblAngleVelocity = uicontrol('style','text');
    
    set(lblTime,'Position', [0 60 90 20]);
    set(lblAction,'Position', [0 30 90 20]);
    set(lblReward,'Position', [0 0 90 20]);    
    set(lblPosition,'Position', [90 60 90 20]);
    set(lblVelocity,'Position', [90 30 90 20]);
    set(lblAcceleration,'Position', [90 0 90 20]);    
    set(lblAngle,'Position', [180 60 90 20]);
    set(lblAngleVelocity,'Position', [180 30 90 20]);
    
    T = length(traj);
    title(['time steps: ', num2str(T)], ' --- final reward: ', num2str(traj{end}.cumReward));
    
    for i = 1:T
        x = traj{i}.state.position;
        omega = traj{i}.state.angle;
        
        pxp=[x x+l*sin(omega)];
        pyp=[1.25 1.25+l*cos(omega)];

        cart.Position = [x-1, 0.25, 2, 1];
        pole.XData = pxp;
        pole.YData = pyp;
        
        set(lblTime,'String', ['t: ',num2str(i)]);
        set(lblAction,'String',  ['a: ',num2str(traj{i}.action)]);
        set(lblReward,'String', ['cr: ',num2str(traj{i}.cumReward)]);

        set(lblPosition,'String', ['p: ',num2str(traj{i}.state.position)]);
        set(lblVelocity,'String',  ['v: ',num2str(traj{i}.state.velocity)]);
        set(lblAcceleration,'String', ['acc: ',num2str(traj{i}.state.acceleration)]);
        
        set(lblAngle,'String',  ['ang: ',num2str(traj{i}.state.angle)]);
        set(lblAngleVelocity,'String', ['angV: ',num2str(traj{i}.state.angleVelocity)]);
        
        pause(0.1); %0.02
    end
end
