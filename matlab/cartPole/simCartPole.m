function [sNext, reward, finished] = simCartPole(state, action, simOpts)
    bounds = simOpts.bounds;
    
    x          = state.position;
    x_dot      = state.velocity;
    omega      = state.angle;
    omega_dot  = state.angleVelocity;

    g               = 9.8;      %Gravity
    Mass_Cart       = 1.0;      %Mass of the cart is assumed to be 1Kg
    Mass_Pole       = 0.1;      %Mass of the pole is assumed to be 0.1Kg
    Total_Mass      = Mass_Cart + Mass_Pole;
    Length          = 0.5;      %Half of the length of the pole 
    PoleMass_Length = Mass_Pole * Length;
    Force_Mag       = 10;
    Tau             = 0.02;     %Time interval for updating the values
    Fourthirds      = 4.0/3.0;

    action = bound(action, bounds.action);
    force = action * Force_Mag;
    
    sinOmega = sin(omega);
    cosOmega = cos(omega);
    
    temp     = (force + PoleMass_Length * omega_dot * omega_dot * sinOmega) / Total_Mass;
    omegaacc = (g * sinOmega - cosOmega * temp) / (Length * (Fourthirds - Mass_Pole * cosOmega * cosOmega / Total_Mass));
    xacc     = temp - PoleMass_Length * omegaacc * cosOmega / Total_Mass;
        
    % Update the four state variables, using Euler's method.
    x         = x + Tau * x_dot;
    x_dot     = x_dot + Tau * xacc;   
    omega     = omega + Tau * omega_dot;
    omega_dot = omega_dot+Tau*omegaacc;

    sNext.position = x;
    sNext.velocity = x_dot;
    sNext.acceleration = xacc;
    sNext.angle = omega;
    sNext.angleVelocity = omega_dot;

    reward = 1;
    finished = false;

    if ~inBounds(x, bounds.position) || ~inBounds(x, bounds.angle)
        finished = true;
%     else
%         if inBounds(x, bounds.rewardPosition)
%         	reward = reward + 1;
%         else
%             reward = reward - 1;
%         end
%         if inBounds(omega, bounds.rewardAngle)  %in angular reward bounds
%             reward = reward + 1;
%         else
%             reward = reward - 1;
%         end
    end
end

function ret = inBounds(x, bounds)
	minBound = min(bounds);
    maxBound = max(bounds);
    ret = x > minBound && x < maxBound;
end

function ret = bound(x, bounds)
    minBound = min(bounds);
    maxBound = max(bounds);
    ret = min(max(x, minBound), maxBound);
end