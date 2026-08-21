function [b] = getMagVec(r, t, params)
%getMagVec Calculate the magnetic field vector at a given point in ECI coordinates.
% Uses the simple dipole model used in the paper.
    arguments (Input)
        r (3,:) double  % [km]
        t (1,:) double  % [s]
        params.Me      (1,1) double = 7.8379e6   % [T km^3] earth magnetic dipole moment
        params.gamma_m (1,1) double = 11.44      % [deg] tilt angle
        params.beta_m  (1,1) double = 0          % [deg] initial plane angle of the dipole
        params.we      (1,1) double = 7.292e-5   % [rad/s] angular rate of earth
    end
    
    arguments (Output)
        b (3,:) double
    end

    Me = params.Me; gamma_m = params.gamma_m; beta_m = params.beta_m; we = params.we;

    phase = beta_m + rad2deg(we*t);   % Current phase of the dipole
    rmag = vecnorm(r,2,1);

    d = -Me*[sind(gamma_m).*sind(phase); -sind(gamma_m).*cosd(phase); cosd(gamma_m)];   % dipole vector
    b = 1./rmag.^3 .* (3*dot(d,r).*r./rmag.^2 - d);                                    % magnetic field vector  
end