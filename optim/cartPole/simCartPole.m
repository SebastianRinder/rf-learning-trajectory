function sNext = simCartPole(state, action, ~)
    %position
    %velocity
    %acceleration
    %angle
    %angularVelocity
    x          = state(1);
    x_dot      = state(2);
    omega      = state(4);
    omega_dot  = state(5);

    g               = 9.8;      %Gravity
    Mass_Cart       = 1.0;      %Mass of the cart is assumed to be 1Kg
    Mass_Pole       = 0.1;      %Mass of the pole is assumed to be 0.1Kg
    Total_Mass      = Mass_Cart + Mass_Pole;
    Length          = 0.5;      %Half of the length of the pole 
    PoleMass_Length = Mass_Pole * Length;
    Force_Mag       = 10;
    Tau             = 0.02;     %Time interval for updating the values
    Fourthirds      = 4.0/3.0;

%     action = applyBound(action, bounds.action);
    force = action * Force_Mag;
    
    sinOmega = sin(omega);
    cosOmega = cos(omega);
    
    temp     = (force + PoleMass_Length * omega_dot * omega_dot * sinOmega) / Total_Mass;
    omega_dot_dot = (g * sinOmega - cosOmega * temp) / (Length * (Fourthirds - Mass_Pole * cosOmega * cosOmega / Total_Mass));
    x_dot_dot     = temp - PoleMass_Length * omega_dot_dot * cosOmega / Total_Mass;
        
    % Update the four state variables, using Euler's method.
    x         = x + Tau * x_dot;
    x_dot     = x_dot + Tau * x_dot_dot;   
    omega     = omega + Tau * omega_dot;
    omega_dot = omega_dot + Tau * omega_dot_dot;

    sNext(1) = x;
    sNext(2) = x_dot;
    sNext(3) = x_dot_dot;
    sNext(4) = omega;
    sNext(5) = omega_dot; 
end

