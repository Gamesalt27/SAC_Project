function [R_I2O] = getECI2O(r,v)
%getECI2O Calculate the rotation matrix between ECI and orbital frame.
arguments (Input)
    r (3,:) double
    v (3,:) double
end

arguments (Output)
    R_I2O (3,3,:) double
end
    z_OI = -r./vecnorm(r,2,1);         % towards earth
    h = cross(r, v, 1); hmag = vecnorm(h, 2, 1);
    y_OI = -h ./ hmag;                  % opposite of angular momentum
    x_OI = cross(y_OI, z_OI, 1);       % completes right handed system

    R_I2O = [reshape(x_OI,1,3,[]);reshape(y_OI,1,3,[]);reshape(z_OI,1,3,[])];    % DCM from inertial frame to orbit (as in the paper)
end