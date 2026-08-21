function [r1, v1] = propOrbit(r0, v0, dt, params)
%propOrbit Propagate a satellite orbit using the Lagrange functions method.
arguments (Input)
    r0 (3,1) double       
    v0 (3,1) double
    dt (1,:) double
    params.mu  (1,1) double = 3.986e+5
    params.tol (1,1) double = 12
end

arguments (Output)
    r1 (3,:) double
    v1 (3,:) double
end
    
    % Find eccentric anomly
    elems = calcOrbitalElems(r0, v0, mu=params.mu, tol=params.tol);
    M1    = mod(sqrt(params.mu/elems.a^3).*dt + elems.M, 2*pi);
    E1    = rad2deg(getE(M1, elems.e, params.tol));

    % Lagrange functions
    F  = 1 - elems.a/norm(r0) .* (1 - cosd(E1-elems.E));
    G  = elems.a/params.mu * dot(r0, v0) .* (1 - cosd(E1-elems.E)) +... 
        norm(r0).*sqrt(elems.a/params.mu).*sind(E1-elems.E);
    r1 = F.*r0 + G.*v0;

    Ft = -sqrt(params.mu*elems.a)./(vecnorm(r1, 2, 1).*norm(r0)).*sind(E1-elems.E);
    Gt = 1-elems.a./vecnorm(r1, 2, 1) .* (1 - cosd(E1-elems.E));
    v1 = Ft.*r0 + Gt.*v0;

    
end

function [E] = getE(M, e, tol)
%getE Finds the eccentric anomaly from the mean anomaly using
%Newton-Raphson. Works for vector inputs.
    if nargin < 3, tol = 12; end
    E = M + e.*sin(M) ./ (1 - sin(M + e) + sin(M)); % Cubic approximation by Mikkola
    for i=1:50
        dE = (M + e.*sin(E) - E) ./ (e.*cos(E) - 1);
        E = E - dE;
        if all(abs(dE) < 10^(-tol)), break; end
    end

end
