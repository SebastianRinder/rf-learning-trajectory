function visCartPole(traj, bounds)
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
    set(lblTime,'Position', [1 20 90 20]);
    set(lblAction,'Position', [1 50 90 20]);
    set(lblReward,'Position', [1 80 90 20]);
    
    for i = 1:size(traj,1)
        x = traj{i,1}.state.position;
        omega = traj{i,1}.state.angle;
        
        pxp=[x x+l*sin(omega)];
        pyp=[1.25 1.25+l*cos(omega)];

        cart.Position = [x-1, 0.25, 2, 1];
        pole.XData = pxp;
        pole.YData = pyp;
        
        set(lblTime,'String', ['t: ',int2str(i)]);
        set(lblAction,'String',  ['a: ',num2str(traj{i,1}.action)]);
        set(lblReward,'String', ['cr: ',num2str(traj{i,1}.cumReward)]);

        pause(0.02);
    end
end
