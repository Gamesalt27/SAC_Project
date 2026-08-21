function [b] = getMagVec_V2(r, t, params)
%getMagVec Calculate the magnetic field vector at a given point in ECI coordinates.
% Uses the simple dipole model used in the paper.
    arguments (Input)
        r (3,:) double  % [km]
        t (1,:) double  % [s]
        params.Me       (1,1) double = 7.8379e6   % [T km^3] earth magnetic dipole moment
        params.gamma_m  (1,1) double = 11.44      % [deg] tilt angle
        params.eta_m    (1,1) double = 0          % [deg] initial ascending node of the dipole
        params.we       (1,1) double = 7.292e-5   % [rad/s] angular rate of earth
        params.T        (1,1) double = 5855       % [s] orbital period
    end
    
    arguments (Output)
        b (3,:) double
    end

    Me = params.Me; gamma_m = params.gamma_m; eta_m = params.eta_m; we = params.we;

    bmag = Me./vecnorm(r,2,2).^3;
    phase = we*t + eta_m;
    

    
end