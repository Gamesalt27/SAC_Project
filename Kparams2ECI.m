function [r_ECI, v_ECI] = Kparams2ECI(orbital_params)
%Kparam2ECI Generate the position and velocity assuming perfect spherical earth from
% Keplerian orbital parameters
%
%   Inputs:
%       orbital_params - Struct containing:
%           .a      : [km] Semi-major axis 
%           .e      : Eccentricity [0 <= e < 1]
%           .i      : [deg] Inclination angle
%           .omega  : [deg] Argument of periapsis
%           .Omega  : [deg] Longitude of ascending node
%           .f      : [deg] true anomaly 
%           .mu     : (Optional) Gravitational parameter. Default Earth.

arguments (Input)
    orbital_params (1,1) struct
end
arguments (Output)
    r_ECI (3,1) double
    v_ECI (3,1) double
end

% Check that all required fields are present
required_fields = ["a", "e", "i", "omega", "Omega", "f"];
for k = 1:numel(required_fields)
    if ~isfield(orbital_params, required_fields(k))
        error("Kparams2ECI:MissingParam", ...
            "The input struct must contain the field '%s'.", required_fields(k));
    end
end

% Assign defaults for optional fields if they are missing
if ~isfield(orbital_params, 'mu')
    orbital_params.mu = 3.986e5;
end

% Extract variables and convert degrees to radians
a = orbital_params.a;
e = orbital_params.e;
inc = deg2rad(orbital_params.i);
f = deg2rad(orbital_params.f);
if e ~= 0
    w = deg2rad(orbital_params.omega);
else
    w = 0;
    if orbital_params.omega ~= 0
        warning("In circular orbits argument of perigee is undefined. Ignoring provided value.")
    end
end
if inc ~= 0
    W = deg2rad(orbital_params.Omega);
else
    W = 0;
    if orbital_params.Omega ~= 0
        warning("In equatorial orbits longitude of the ascending node is undefined. Ignoring provided value.")
    end
end
mu = orbital_params.mu;

p = a * (1 - e^2);

% Calculate position and velocity in Perifocal Frame
r_Pf = zeros(3, 1);
r_Pf(1, 1) = p / (1 + e * cos(f)) * cos(f);
r_Pf(2, 1) = p / (1 + e * cos(f)) * sin(f);

v_Pf = zeros(3, 1);
v_Pf(1, :) = sqrt(mu / p) * (-sin(f));
v_Pf(2, :) = sqrt(mu / p) * (e + cos(f));

% Transform vectors to the ECI frame
R_PF2ECI = eul2rotm([W,inc,w],"ZXZ");         
r_ECI = R_PF2ECI * r_Pf;
v_ECI = R_PF2ECI * v_Pf;

end


