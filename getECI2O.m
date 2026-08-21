function [R_I2O] = getECI2O(r,v)
%getECI2O Calculate the rotation matrix between ECI and orbital frame.
arguments (Input)
    r (3,1) double
    v (3,1) double
end

arguments (Output)
    R_I2O (3,3) double
end
    z_OI = -r/norm(r);                                    % towards earth 
    y_OI = -cross(r, v) / norm(cross(r, v));              % opposite of angular momentum
    x_OI = cross(y_OI, z_OI) / norm(cross(y_OI, z_OI));   % normalized to reduce errors

    R_I2O = [x_OI.';y_OI.';z_OI.'];                       % DCM from inertial frame to orbit (as in the paper)
end