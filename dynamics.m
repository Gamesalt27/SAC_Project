function [x, b, m, M] = dynamics(t, eta, eps, omega, r0, v0, params)
%DYNAMICS Non linear dyanimcs of the controller
    arguments (Input)
        t     (1,:) double  % [s] time
        eta   (1,:) double  % [-] quaternion scalar B to LVLH
        eps   (3,:) double  % [-] quaternion vector B to LVLH
        omega (3,:) double  % [rad/s] angular rate B to I in body
        r0    (3,1) double  % [km] initial position in ECI
        v0    (3,1) double  % [km/s] initial velocity in ECI
        params.kw   (1,1) double      % [-] control gain
        params.I    (3,3) double      % [kg m^2] inertia tensor
        params.mlim (1,1) double = 2  % [A m^2] actuator saturation limit
    end
    
    arguments (Output)
        x (7,:) double   % State vector
        b (3,:) double   % Body frame field
        m (3,:) double   % Actuator fields
        M (3,:) double   % Total torque
    end

    % Unpacking
    I = params.I; kw = params.kw; mlim = params.mlim;

    % Normlize the quaternion to 1
    Qnorm = vecnorm([eta; eps], 2, 1);   
    eta = eta./Qnorm; eps = eps./Qnorm;

    % Find magnetic field
    [r, v] = propOrbit(r0, v0, t);
    b_ECI = reshape(getMagVec(r, t),3,1,[]);
    R_I2O = getECI2O(r, v);             % DCM from inertial frame to orbit
    R_O2B = quat2dcm([eta; eps].');     % DCM from orbit to body    
    b = pagemtimes(R_O2B,pagemtimes(R_I2O,b_ECI));   % body frame
    b = reshape(b,3,[]);

    % Control law
    m = -kw./vecnorm(b,2,1).^2 .* cross(b,omega,1);     % actuators strength
    m(abs(m) > mlim) = sign(m(abs(m) > mlim))*mlim;     % saturation limit
    M = cross(m,b,1);

    % Angular velocity dynamics 
    omegadot = I \ (M - cross(omega,I*omega,1));
    
    % Quaternion dynamics
    omega_OII = reshape(cross(r,v,1)./vecnorm(r,2,1).^2,3,1,[]);    % orbit angular rate in inertial frame
    omega_OIB = pagemtimes(R_O2B,pagemtimes(R_I2O,omega_OII));     % in body frame
    omega_OIB = reshape(omega_OIB,3,[]);
    omega_BOB = omega - omega_OIB;       % angular rate needed for quaternion propagation
    epsdot = 0.5*(eta.*omega_BOB + cross(eps,omega_BOB,1));
    etadot = -0.5*dot(eps,omega_BOB,1);
    qdot = normqdot([eta; eps], [etadot; epsdot]);

    x = [qdot; omegadot];
end

function [qdot] = normqdot(q, qdot)
%NORMQDOT remove the residual part of qdot that changes the length of q. 
arguments (Input)
    q    (4,:) double   % quaternion, assumed already normalized
    qdot (4,:) double   % quaternion derivative
end

arguments (Output)
    qdot (4,:) double   % quaternion derivative with size preservation
end

    res = dot(qdot,q,1).*q;   % projection of qdot onto q
    qdot = qdot-res;

end