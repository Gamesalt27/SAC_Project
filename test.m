clc; clearvars; close all;

%%

% case C
% Omega = 272:0.001:273; 
% u = 35:0.001:37;
% i = 65;
% Q = [0.646  0.525 -0.514  0.206];
% omega = [1 0 0].';

% case A
Omega = 0:1:360; 
u = 0:1:360;
i = 11;
Q = [0.375 -0.062  0.925 -0.007];
omega = [0.604 -0.76 -0.384];
r_c = 7021;
Q = Q/norm(Q);
mu = 3.986e5;

[OmegaGrid, uGrid] = meshgrid(Omega, u);
Omega = OmegaGrid(:).';
u = uGrid(:).';

N = numel(Omega);

cO = cosd(Omega);
sO = sind(Omega);
cu = cosd(u);
su = sind(u);
ci = cosd(i);
si = sind(i);

r = r_c*[ ...
     cO.*cu - sO.*su.*ci;
     sO.*cu + cO.*su.*ci;
     su.*si ];
v = sqrt(mu/r_c)*[ ...
    -cO.*su - sO.*cu.*ci;
    -sO.*su + cO.*cu.*ci;
     cu.*si ];

b_ECI = reshape(getMagVec(r, zeros(1,length(r))),3,1,[]);
R_I2O = getECI2O(r, v);
R_O2B = quat2dcm(Q);
% b = pagemtimes(R_O2B+zeros(size(R_I2O)),pagemtimes(R_I2O,b_ECI));
b = pagemtimes(R_I2O,b_ECI);
b = reshape(b,3,[]);

[alignment, idx] = min(vecnorm(cross([0;-1;0]+zeros(size(b)),b),2,1));
disp(Omega(idx))
disp(u(idx))
disp(alignment)
