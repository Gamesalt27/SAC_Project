% clearvars; close all; clc;

%% PAPER SIMULATIONS
disp("============== Simulations From The Paper ==============")

%% data from table 1

J1 = 0.33; J2 = 0.37; J3 = 0.35;    % [kg m^2] Principal moments of inertia
I = diag([J1, J2, J3]);

mlim    = 2;        % [A m^2] acuator saturation limit
r_c     = 7021;     % [km] orbit radius
To      = 5855;     % [s] orbit period
gamma_m = 11.44;    % [deg] geomagnetic plane tilt
beta_m  = 0;        % [deg] initial geomagnetic phase
mu      = 3.986e5;  % [km^3/s^2] earth gravitational parameter
Re      = 6378;     % [km] earth radius
N       = 5;        % [-] number of orbits
n       = 1e3;      % [-] points per orbit

% Cases A, B, and C
i      = [11, 65, 65];                     % [deg] inclination
omega0 = [0.604 -0.76 -0.384; 
          0.604 -0.76 -0.384; 
          1      0     0].';               % [rad/s] initial angular velocity b to I (body frame)
Q0     = [0.375 -0.062  0.925 -0.007; 
          0.375 -0.062  0.925 -0.007; 
          0.646  0.525 -0.514  0.206].';   % [-] initial quaternions (LVLH to body)
Q0 = Q0./vecnorm(Q0, 2, 1);

Omega0  = [30, 272, 272];   % inital ascending node
arglat0 = [0, 36, 36];    % initial argument of latitude



%% Simulation

% results = cell(1,length(i));
opts = odeset(Refine=1, Stats='on');
tspan = [0, N*To];

for j=1:length(i)
    fprintf("Trial number %d\n",j);
    R_O2I = angle2dcm(deg2rad(Omega0(j)),deg2rad(i(j)),deg2rad(arglat0(j)), "ZXZ").';   % LVLH to ECI
    r0 = R_O2I*[r_c;0;0];
    v0 = R_O2I*[0;sqrt(mu/r_c);0];
    xi_m = acosd( cosd(i(j))*cosd(gamma_m) + sind(i(j))*sind(gamma_m)*cosd(beta_m-Omega0(j)));
    kw = 4*pi/To * (1+sind(xi_m)) * min(diag(I));

    fun = @(t,x) dynamics(t, x(1), x(2:4), x(5:7), r0, v0, I=I, kw=kw);
    x0 = [Q0(:,j); omega0(:,j)].';
    [Time, x] = ode89(fun, tspan, x0, opts);
    Q = x(:,1:4); omega = x(:,5:end); 
    [rates,b,m,M] = dynamics(Time, Q(:,1), Q(:,2:4).', omega.', r0, v0, I=I, kw=kw);
    b = b.'; m = m.'; M = M.';
    results(j) = {table(Time, Q, omega, b, m, M)};
end


%% Plots

handle = figure('Position', [100, 100, 600, 400]);
hold on; grid on;
for j=1:length(results)
    plot(results{j}.Time/To, vecnorm(results{j}.omega,2,2), LineWidth=1.5);
end
hold off; title("Angular rate over time"); xlabel('t/T [-]'); xlim([0,5])
ylabel('||\omega|| [rad/s]'); legend('||\omega_A||', '||\omega_B||', '||\omega_C||')

% handle = figure('Position', [100, 100, 600, 400]);
% hold on; grid on;
% for j=1:length(results)
%     plot(results{j}.Time/To, vecnorm(results{j}.Q,2,2), LineWidth=1.5);
% end
% hold off; title("Quaternion over time"); xlabel('t/T [-]');
% ylabel('||q||'); legend('||q_A||', '||q_B||', '||q_C||')
% 
% handle = figure('Position', [100, 100, 600, 400]);
% tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
% nexttile;
% hold on; grid on;
% plot(results{1}.Time/To, results{1}.omega(:,1), 'b', 'LineWidth', 1.5);
% title('\omega_1');
% xlabel('t/T [-]');
% ylabel('\omega_1 [rad/s]');
% hold off;
% 
% nexttile;
% hold on; grid on;
% plot(results{1}.Time/To, results{1}.omega(:,2), 'b', 'LineWidth', 1.5);
% title('\omega_2');
% xlabel('t/T [-]');
% ylabel('\omega_2 [rad/s]');
% hold off;
% 
% nexttile;
% hold on; grid on;
% plot(results{1}.Time/To, results{1}.omega(:,3), 'b', 'LineWidth', 1.5);
% title('\omega_3');
% xlabel('t/T [-]');
% ylabel('\omega_3 [rad/s]');
% hold off;
% 
