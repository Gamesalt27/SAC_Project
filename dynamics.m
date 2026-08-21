function [x] = dynamics(t, eta, eps, omega, r0, v0, params)
%DYNAMICS Non linear dyanimcs of the controller
    arguments (Input)
        t     (1,1) double  % [s] time
        eta   (1,1) double  % [-] quaternion scalar B to LVLH
        eps   (3,1) double  % [-] quaternion vector B to LVLH
        omega (3,1) double  % [rad/s] angular rate B to I in body
        r0    (3,1) double  % [km] initial position in ECI
        v0    (3,1) double  % [km/s] initial velocity in ECI
        params.kw   (1,1) double      % [-] control gain
        params.I    (3,3) double      % [kg m^2] inertia tensor
        params.mlim (1,1) double = 2  % [A m^2] actuator saturation limit
    end
    
    arguments (Output)
        x (7,1) double % State vector
    end

    % Unpacking
    I = params.I; kw = params.kw; mlim = params.mlim;

    % Normlize the quaternion to 1
    Qnorm = norm([eta; eps]);   
    eta = eta/Qnorm; eps = eps/Qnorm;

    % Find magnetic field
    [r, v] = propOrbit(r0, v0, t);
    b_ECI = getMagVec(r, t);
    R_I2O = getECI2O(r, v);             % DCM from inertial frame to orbit
    R_O2B = quat2dcm([eta; eps].');     % DCM from orbit to body    
    b = R_O2B*R_I2O*b_ECI;              % body frame

    % Control law
    m = -kw/norm(b)^2 * cross(b,omega);                 % actuators strength
    m(abs(m) > mlim) = sign(m(abs(m) > mlim))*mlim;     % saturation limit
    M = cross(m, b);

    % Angular velocity dynamics 
    omegadot = I \ (M - cross(omega, I*omega));
    
    % Quaternion dynamics
    omega_OII = cross(r,v)/norm(r)^2;    % orbit angular rate in inertial frame
    omega_OIB = R_O2B*R_I2O*omega_OII;   % in body frame
    omega_BOB = omega - omega_OIB;       % angular rate needed for quaternion propagation
    epsdot = 0.5*(eta*omega_BOB + cross(eps,omega_BOB));
    etadot = -0.5*dot(eps,omega_BOB);
    qdot = normqdot([eta; eps], [etadot; epsdot]);

    x = [qdot; omegadot];
end

function [qdot] = normqdot(q, qdot)
%NORMQDOT remove the residual part of qdot that changes the length of q. 
arguments (Input)
    q    (4,1) double   % quaternion, assumed already normalized
    qdot (4,1) double   % quaternion derivative
end

arguments (Output)
    qdot (4,1) double   % quaternion derivative with size preservation
end

    res = dot(qdot, q)*q / (norm(q)^2);   % projection of qdot onto q
    qdot = qdot-res;

end