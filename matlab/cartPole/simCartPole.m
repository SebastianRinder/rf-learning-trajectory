function [sNext, reward, toBreak] = simCartPole(state, action, bounds)
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

    %action = min(max(action, bounds.action(1,1)), bounds.action(1,2));  %apply action bounds
    force = action * Force_Mag;

    temp     = (force + PoleMass_Length * omega_dot * omega_dot * sin(omega)) / Total_Mass;
    omegaacc = (g * sin(omega) - cos(omega) * temp) / (Length * (Fourthirds - Mass_Pole * cos(omega) * cos(omega) / Total_Mass));
    xacc     = temp - PoleMass_Length * omegaacc * cos(omega) / Total_Mass;

    % Update the four state variables, using Euler's method.
    x         = x + Tau * x_dot;
    x_dot     = x_dot + Tau * xacc;
    omega     = omega + Tau * omega_dot;
    omega_dot = omega_dot+Tau*omegaacc;

%    sNext = [x, x_dot, xacc, omega, omega_dot];
    sNext.position = x;
    sNext.velocity = x_dot;
    sNext.acceleration = xacc;
    sNext.angle = omega;
    sNext.angleVelocity = omega_dot;

    reward = 0;
    toBreak = false;

    if x < bounds.position(1,1) || x > bounds.position(1,2) ||...
            omega < bounds.angle(1,1) || omega > bounds.angle(1,2) %out of bounds
        toBreak = true;
    else
        if x > bounds.rewardPosition(1,1) && x < bounds.rewardPosition(1,2) %in position reward bounds
        	reward = reward + 1;
        else
            reward = reward - 1;
        end
        if omega > bounds.rewardAngle(1,1) && omega < bounds.rewardAngle(1,2)  %in angular reward bounds
            reward = reward + 1;
        else
            reward = reward - 1;
        end
    end
end