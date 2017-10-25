function visAcroBot(traj, ~)
    skipVisValue = selectFigure('Visualize Acro Bot');
    clf;
    axis([-2.1 2.1 -2.1 2.1])
    axis square
    grid on
    hold on;
    
    T = length(traj.action);
    title(['time steps: ', num2str(T), '      final reward: ', num2str(traj.cumReward(end))]);
    
    checkBoxSkip = uicontrol('Style', 'checkbox', 'String', 'Skip Vis',...
        'Position', [90 0 90 20],...
        'Value', skipVisValue);
    
    lblTime = uicontrol('style','text');
    set(lblTime,'Position', [0 0 90 20]);
    
    pole1 = line;
    pole2 = line;
    
    x_acrobot(1)=0;
    y_acrobot(1)=0;
    
    for i = 1:T    
        if get(checkBoxSkip, 'Value') == 1
            break;
        end
        
        theta1 = traj.state(i,1);
        theta2 = traj.state(i,2);    

        x_acrobot(2) = x_acrobot(1) + sin(theta1); 
        y_acrobot(2) = y_acrobot(1) - cos(theta1);
        x_acrobot(3) = x_acrobot(2) + sin(theta2); 
        y_acrobot(3) = y_acrobot(2) - cos(theta2);
        
        pole1.XData = [x_acrobot(1) x_acrobot(2)];
        pole1.YData = [y_acrobot(1) y_acrobot(2)];
        
        pole2.XData = [x_acrobot(2) x_acrobot(3)];
        pole2.YData = [y_acrobot(2) y_acrobot(3)];

%         plot(x_acrobot ,y_acrobot,'ok-','LineWidth',1,'markersize',7,'MarkerFaceColor',[.7 .7 .7]);
%         
%         plot(x_acrobot(3),y_acrobot(3),'.r','markersize',20);
        
        set(lblTime,'String', ['t: ',num2str(i)]);
        
        pause(0.05);
    end
    
    %drawnow
    hold off
end