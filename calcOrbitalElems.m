function [elems] = calcOrbitalElems(r, v, params)
%calcOrbitalElems Calculate Keplerian orbital parameters from ECI state vectors.
% Doesn't currently support parabolic or hyperbolic orbits, will output
% garbage values.
arguments (Input)
    r (3,:) double
    v (3,:) double
    params.mu (1,1) double = 3.986e+5
    params.tol (1,1) double = 12
end
arguments (Output)
    elems (1,1) struct
end


% Helper function to prevent values outside range of acos due to numerical errors.
clamp = @(x) max(min(x, 1), -1);

h = cross(r, v);
elems.h = vecnorm(h, 2, 1);
n = zeros(size(r));
n(1:2,:) = [-h(2,:); h(1,:)];

r_mag = vecnorm(r, 2, 1);
h_mag = vecnorm(h, 2, 1);
n_mag = vecnorm(n, 2, 1);

e_vec = cross(v, h)./params.mu - r./r_mag; 
elems.e = vecnorm(e_vec, 2, 1);

% Initialize angles
elems.i = acosd(clamp(h(3,:)./h_mag));
elems.Omega = zeros(1, size(r, 2));
elems.omega = zeros(1, size(r, 2));
elems.f = zeros(1, size(r, 2));
elems.E = zeros(1, size(r, 2));

has_inc = abs(h_mag - abs(h(3,:))) > 10^(-params.tol);     % Non-equatorial orbit
has_ecc = elems.e > 10^(-params.tol);    % Eliptic orbit

% For eliptic non-equatorial orbits
elip_non_eq = has_ecc & has_inc;
elems.Omega(has_inc) = rad2deg(atan2(n(2,has_inc), n(1,has_inc)));  % Doesn't depend on eccentricity 
elems.Omega(has_inc & n(2,:) < 0) = 360 + elems.Omega(has_inc & n(2,:) < 0);  % Correct to desired angle behaviour.

elems.omega(elip_non_eq) = acosd(clamp( dot(n(:,elip_non_eq), e_vec(:,elip_non_eq)) ...
    ./ (n_mag(elip_non_eq) .* elems.e(elip_non_eq)) ));
elems.omega(elip_non_eq & e_vec(3,:) < 0) = 360 - elems.omega(elip_non_eq & e_vec(3,:) < 0);    % Correct to desired angle behaviour.

rdotv = dot(r, v);
elems.f(has_ecc) = acosd(clamp(dot(e_vec(:,has_ecc), r(:,has_ecc)) ./ (elems.e(has_ecc) .* r_mag(has_ecc)))); % Doesn't depend on inclination
elems.f(has_ecc & rdotv < 0) = 360 - elems.f(has_ecc & rdotv < 0);  % Correct to desired angle behaviour.

elems.E(has_ecc) = acosd( (elems.e+cosd(elems.f))./(1+elems.e.*cosd(elems.f)) );
elems.E(has_ecc & rdotv < 0) = 360 - elems.E(has_ecc & rdotv < 0);     % Correct to desired angle behaviour.

% For eliptic equatorial orbits
elip_eq = has_ecc & ~has_inc;
elems.omega(elip_eq) = rad2deg(atan2(e_vec(2,elip_eq), e_vec(1,elip_eq)));
elems.omega(elip_eq & e_vec(2,:) < 0) = 360 + elems.omega(elip_eq & e_vec(2,:) < 0);   % Correct to desired angle behaviour.

% For circular non-equatorial orbits
circ_non_eq = ~has_ecc & has_inc;
elems.f(circ_non_eq) = acosd(clamp(dot(n(:,circ_non_eq), r(:,circ_non_eq)) ./ ...
    (n_mag(circ_non_eq) .* r_mag(circ_non_eq))));

% For circular and equatorial orbit
circ_eq = ~(has_ecc | has_inc);
elems.f(circ_eq) = acosd(clamp(r(1,circ_eq) ./ r_mag(circ_eq)));
elems.f(~has_ecc & r(2,:) < 0) = 360 - elems.f(~has_ecc & r(2,:) < 0);    % Correct to desired angle behaviour.

elems.E(~has_ecc) = elems.f(~has_ecc);
elems.M = deg2rad(elems.E)-elems.e.*sind(elems.E);

elems.a = h_mag.^2 ./ (params.mu .* (1 - elems.e.^2));
elems.eps = -params.mu./(2*elems.a);
elems.T = 2*pi.*sqrt(elems.a.^3 ./ params.mu);

% Round everything up to tolarence to clean up floating point errors
fields = fieldnames(elems);
for i = 1:numel(fields)
    val = elems.(fields{i});
    if isnumeric(val)
        elems.(fields{i}) = round(val, params.tol);
    end
end

end
