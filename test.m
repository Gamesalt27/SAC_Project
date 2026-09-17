clc; clearvars; close all;

%%

% case C
% Omega = 272:0.001:273; 
% u = 35:0.001:37;
% i = 65;
% Q = [0.646  0.525 -0.514  0.206];
% omega = [1 0 0].';

% case A
Omega = 0:2:360; 
u = 0:2:360;
i = 11;
Q = [0.375 -0.062  0.925 -0.007];
omega = [0.604 -0.76 -0.384];
omega_hat = omega/norm(omega);
r_c = 7021;
Q = Q/norm(Q);
mu = 3.986e5;

[OmegaGrid, uGrid] = meshgrid(Omega, u);
OmegaGrid = OmegaGrid(:).';
uGrid = uGrid(:).';

N = numel(OmegaGrid);

cO = cosd(OmegaGrid);
sO = sind(OmegaGrid);
cu = cosd(uGrid);
su = sind(uGrid);
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

b_ECI = getMagVec(r, zeros(1,length(r)));
bhat_ECI = reshape(b_ECI./vecnorm(b_ECI,2,1),3,1,[]);
R_I2O = getECI2O(r, v);
R_O2B = quat2dcm(Q);
bhat = pagemtimes(R_O2B+zeros(size(R_I2O)),pagemtimes(R_I2O,bhat_ECI));
% bhat = pagemtimes(R_I2O,b_ECI);
bhat= reshape(bhat,3,[]);

score = abs(dot(omega.'+zeros(size(bhat)),bhat) - cos(5/8*pi));
[alignment, idx] = min(score);
disp(OmegaGrid(idx))
disp(uGrid(idx))
disp(alignment)

score = reshape(score,length(Omega),length(u));
