function [reward, finished] = rewardAcroBot(state, ~, ~)
    theta(1:2) = state(1:2);
    y_acrobot(1) = 0;
    y_acrobot(2) = y_acrobot(1) - cos(theta(1));
    y_acrobot(3) = y_acrobot(2) - cos(theta(2));
    goal = y_acrobot(1) + 1.0;

    reward = -1; %y_acrobot(3);
    finished = false;

    if(y_acrobot(3) >= goal) 
        reward = 100; %10*y_acrobot(3);        
        finished = true;
    end

end

