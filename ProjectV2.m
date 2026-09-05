clearvars; close all; clc;

%% Setup

% Which simulation to run:
% paper - replicate paper sim
% ecc - test different eccentricities
% optKw - find optimal kw for scenario A
% NearA - optimal value vs analytic at scenarios like A
% attitudeSweep - optimal vs analytic kw over a grid of initial attitudes
sim_type = "NearA";
parallels = false;  % Run on multiple cores
conv_tol = 0.01; % Convergence tolarence.

J1 = 0.33; J2 = 0.37; J3 = 0.35; % [kg m^2] Principal moments of inertia
I = diag([J1, J2, J3]);

mlim = 2; % [A m^2] acuator saturation limit
r_c = 7021; % [km] orbit radius
To = 5855; % [s] orbit period
gamma_m = 11.44; % [deg] geomagnetic plane tilt
beta_m = 0; % [deg] initial geomagnetic phase
mu = 3.986e5; % [km^3/s^2] earth gravitational parameter
Re = 6378; % [km] earth radius
N = 8; % [-] number of orbits
kwopt_A = 8.340986872382716e-04; % Optimal kw for case A found numerically.


a = (mu*To^2/(4*pi^2))^(1/3); % [km] semi major axis

switch sim_type
case 'paper'
% Cases A, B, and C
i = [11, 65, 65]; % [deg] inclination
trials = length(i);
omega0 = [0.604 -0.76 -0.384;
          0.604 -0.76 -0.384;
          1      0     0].'; % [rad/s] initial angular velocity b to I (body frame)
Q0 = [0.375 -0.062  0.925 -0.007;
      0.375 -0.062  0.925 -0.007;
      0.646  0.525 -0.514  0.206].'; % [-] initial quaternions (LVLH to body)
Q0 = Q0./vecnorm(Q0, 2, 1);

Omega0 = [30, 272, 272]; % [deg] inital ascending node
arglat0 = [0, 36, 36]; % [deg] initial argument of latitude
w0 = zeros(size(i)); % [deg] initial argument of perigee
ecc = zeros(size(i)); % [-] eccentricity
leg_text = struct("omega", ['||\omega_A||';'||\omega_B||';'||\omega_C||'],...
"Q", ['||q_A||';'||q_B||';'||q_C||'],...
"alpha", ['\alpha_A';'\alpha_B';'\alpha_C']);

kw = analyticKw(i, Omega0, gamma_m, beta_m, To, I);
IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Q0, omega0, kw);

    case 'ecc'  % Not a good simulation, should try to stay at a similar average altitude instead of keeping period same
ecc = [0 0.1 0.5 0.85 0.99]; % [-] eccentricity
trials = length(ecc);
i = ones(size(ecc))*65; % [deg] inclination
omega0 = [0.604;-0.76;-0.384]*ones(size(ecc)); % [rad/s] initial angular velocity b to I (body frame)
Q0 = [0.375;-0.062;0.925;-0.007]*ones(size(ecc)); % [-] initial quaternions (LVLH to body)
Q0 = Q0./vecnorm(Q0, 2, 1);

Omega0 = zeros(size(ecc)); % [deg] inital ascending node
arglat0 = zeros(size(ecc)); % [deg] initial argument of latitude
w0 = zeros(size(ecc)); % [deg] initial argument of perigee
leg_text = struct("omega", compose('e=%0.2f',ecc),...
"Q", compose('e=%0.2f',ecc),...
"alpha", compose('e=%0.2f',ecc));

kw = analyticKw(i, Omega0, gamma_m, beta_m, To, I);
IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Q0, omega0, kw);

case 'optKw'
trials = 200;
parallels = true;
i = 11; % [deg] inclination
omega0 = [0.604 -0.76 -0.384].'; % [rad/s] initial angular velocity b to I (body frame)
Q0 = [0.375 -0.062 0.925 -0.007].'; % [-] initial quaternions (LVLH to body)
Q0 = Q0./vecnorm(Q0, 2, 1);

Omega0 = 0; % [deg] inital ascending node
arglat0 = 0; % [deg] initial argument of latitude
w0 = 0; % [deg] initial argument of perigee
ecc = 0; % [-] eccentricity

kw_ref = analyticKw(i, Omega0, gamma_m, beta_m, To, I);
kw = [kw_ref, linspace(0.9*kw_ref, 1.3*kw_ref, trials-1)];
IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Q0, omega0, kw);
leg_text = struct("omega", compose('k_\\omega=%e',kw),...
"Q", compose('k_\\omega=%e',kw),...
"alpha", compose('k_\\omega=%e',kw));

case "NearA"
trials = 8;
i = 11; % [deg] inclination
omega0 = [0.604 -0.76 -0.384;
          0.604 -0.76 -0.384;
         -0.76   0.604 -0.384;
          0.604 -0.76 -0.384].'; % [rad/s] initial angular velocity b to I (body frame)
omega0 = [omega0,omega0];
Q0 = [0.375 -0.062  0.925 -0.007;
      0.646  0.525 -0.514  0.206;
      0.375 -0.062  0.925 -0.007;
      0.375 -0.062  0.925 -0.007].'; % [-] initial quaternions (LVLH to body)
Q0 = Q0./vecnorm(Q0, 2, 1);
Q0 = [Q0, Q0];

Omega0 = 0; % [deg] inital ascending node
arglat0 = [0 0 0 -10]; % [deg] initial argument of latitude
arglat0 = [arglat0,arglat0];
w0 = 0; % [deg] initial argument of perigee
ecc = 0; % [-] eccentricity

kwAnalytic = analyticKw(i, Omega0, gamma_m, beta_m, To, I);
kw = [kwAnalytic*ones(1,trials/2),...
    kwopt_A*ones(1,trials/2)]; 
IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Q0, omega0, kw);
leg_text = struct("omega", compose('k_\\omega=%e',kw),...
"Q", compose('k_\\omega=%e',kw),...
"alpha", compose('k_\\omega=%e',kw));

case "attitudeSweep"
qA = [0.375; -0.062; 0.925; -0.007];
qA = qA / norm(qA);

% Generate equally spaced quaternions relative to qA
axisColat = 15:15:165; % Angle between rotation vectors
axisLon = 0:30:330; % Complementry angle
rotAngle = 15:15:165; % Rotation angle

[axisColatGrid,axisLonGrid,rotAngleGrid] = meshgrid(axisColat,axisLon,rotAngle);
axisColatGrid = reshape(axisColatGrid,1,[]);
axisLonGrid = reshape(axisLonGrid,1,[]);
rotAngleGrid = reshape(rotAngleGrid,1,[]);

% Add pure rotations at the beginning
axisColatGrid = [zeros(size(rotAngle)), axisColatGrid];
axisLonGrid = [zeros(size(rotAngle)), axisLonGrid];
rotAngleGrid = [rotAngle, rotAngleGrid];
axisColat = [0, axisColat];

u = [sind(axisColatGrid).*cosd(axisLonGrid).*sind(rotAngleGrid/2);...
     sind(axisColatGrid).*sind(axisLonGrid).*sind(rotAngleGrid/2);...
     cosd(axisColatGrid).*sind(rotAngleGrid/2)];
dQ = [cosd(rotAngleGrid/2);u]; % Difference Quaternion
Qsample = quatmultiply(dQ.',qA.').';    % Rotate qA by the difference quaternion
Qsample = [Qsample, Qsample];

trials = length(Qsample);
parallels = true;

i = 11;
Omega0 = 0;
arglat0 = 0;
w0 = 0;
ecc = 0;

omega0v = [0.604; -0.76; -0.384];

kwAnalytic = analyticKw(i, Omega0, gamma_m, beta_m, To, I);

kw = [kwAnalytic*ones(1,trials/2), kwopt_A*ones(1,trials/2)];
IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Qsample, omega0v, kw);

leg_text = struct("omega", [], "Q", [], "alpha", []);
end

%% Simulation

results = cell(1,trials);
convIdx = zeros(1, trials);
convTime = zeros(1, trials);

tspan = [0, N*To];
if parallels
    D = parallel.pool.DataQueue;
    startTime = tic;
    afterEach(D, @(~) reportProgress(trials, startTime));
    disp("Simulation Began")
    
    parfor j=1:trials
        ICj = IC(j); 
        fun = @(t,x) dynamics(t, x(1), x(2:4), x(5:7), ICj.r0, ICj.v0, I=I, kw=ICj.kw);
        x0 = [ICj.Q0; ICj.omega0];
        
        [Time, x] = ode89(fun, tspan, x0);
        omega = x(:,5:end);
        lastAbove = find(vecnorm(omega,2,2) >= conv_tol*norm(omega(1,:)), 1, 'last');
        if lastAbove == length(Time)
            fprintf("Didn't converge in %d orbits", N)
        end
        convIdx(j) = min(lastAbove+1,length(Time));
        convTime(j) = Time(convIdx(j))/To;
        send(D, 1);
    end

else

    opts = odeset(Refine=1, Stats='on');
    for j=1:trials
        fprintf("Trial number %d\n",j);
        ICj = IC(j);
        
        fun = @(t,x) dynamics(t, x(1), x(2:4), x(5:7), ICj.r0, ICj.v0, I=I, kw=ICj.kw);
        x0 = [ICj.Q0; ICj.omega0];
        
        [Time, x] = ode89(fun, tspan, x0, opts);
        Q = x(:,1:4); omega = x(:,5:end);
        [rates,b,m,M] = dynamics(Time, Q(:,1), Q(:,2:4).', omega.', ICj.r0, ICj.v0, I=I, kw=ICj.kw);
        b = b.'; m = m.'; M = M.';
        results(j) = {table(Time, Q, omega, b, m, M)};
        lastAbove = find(vecnorm(omega,2,2) >= conv_tol*norm(omega(1,:)), 1, 'last');
        if lastAbove == length(Time)
            fprintf("Didn't converge in %d orbits", N)
        end
        convIdx(j) = min(lastAbove+1,length(Time));
        convTime(j) = Time(convIdx(j))/To;
    end

end

%% Plots

plots = plotting();

plots.angularRate(results, To, leg_text);
plots.quaternionNorm(results, To, leg_text);

l = min(629, trials); % Which trial to plot the actuator dipoles for
plots.actuatorDipoles(results, To, mlim, l);

plots.controlAuthority(results, To, convIdx, leg_text);

if strcmp(sim_type, "optKw")
    plots.settlingTimeVsKw(kw, convTime);
end

if strcmp(sim_type, "NearA")
    plots.nearAComparison(kw, convTime, trials);
end

if strcmp(sim_type, "attitudeSweep")
    [~, improvementPct] = plots.attitudeSweepByTrial(convTime, trials);
    plots.attitudeSweepByAngles(axisColat, axisLon, rotAngle, ...
        axisColatGrid, axisLonGrid, rotAngleGrid, improvementPct, ...
        kwOptimal, kwAnalytic);
end

%% Local functions

function IC = buildIC(a, ecc, incl, Omega0, arglat0, w0, Q0, omega0, kw)
%BUILDIC Generate IC struct

arguments
    a       (1,1) double {mustBePositive}
    ecc     (1,:) double {mustBeNonnegative}
    incl    (1,:) double
    Omega0  (1,:) double
    arglat0 (1,:) double
    w0      (1,:) double
    Q0      (4,:) double
    omega0  (3,:) double
    kw      (1,:) double {mustBePositive}
end

    trials = max([numel(ecc), numel(incl), numel(Omega0), numel(arglat0), ...
        numel(w0), size(Q0,2), size(omega0,2), numel(kw)]);
    
    ecc     = broadcastCols(ecc, trials);
    incl    = broadcastCols(incl, trials);
    Omega0  = broadcastCols(Omega0, trials);
    arglat0 = broadcastCols(arglat0, trials);
    w0      = broadcastCols(w0, trials);
    Q0      = broadcastCols(Q0, trials);
    omega0  = broadcastCols(omega0, trials);
    kw      = broadcastCols(kw, trials);
    
    IC(trials) = struct('r0',[],'v0',[],'Q0',[],'omega0',[],'kw',[]);
    for j = 1:trials
        [r0, v0] = Kparams2ECI(struct("a",a, "e",ecc(j), "i",incl(j), ...
            "omega",w0(j), "Omega",Omega0(j), "f",arglat0(j)-w0(j)));
        IC(j).r0 = r0;
        IC(j).v0 = v0;
        IC(j).Q0 = Q0(:,j);
        IC(j).omega0 = omega0(:,j);
        IC(j).kw = kw(j);
    end
end

function out = broadcastCols(val, trials)
%BROADCASTCOLS like z fill, fills the list with copies of initial element
%   if needed.
arguments
    val    double
    trials (1,1) double {mustBePositive, mustBeInteger}
end
    if size(val, 2) == 1 && trials > 1
        out = repmat(val, 1, trials);
    else
        out = val;
    end
end

function kw = analyticKw(incl, Omega0, gamma_m, beta_m, To, I)
%ANALYTICKW vectorized implementation of the analytical Kw from paper.

arguments
    incl    (1,:) double
    Omega0  (1,:) double
    gamma_m (1,1) double
    beta_m  (1,1) double
    To      (1,1) double {mustBePositive}
    I       (3,3) double
end

    xi_m = acosd( cosd(incl).*cosd(gamma_m) + sind(incl).*sind(gamma_m).*cosd(beta_m-Omega0) );
    kw = 4*pi/To * (1+sind(xi_m)) * min(diag(I));
end

function reportProgress(total, startTime)
%REPORTPROGRESS Helper function for parfor loop progress tracking.
arguments
    total     (1,1) double {mustBePositive, mustBeInteger}
    startTime (1,1) uint64
end
    persistent done
    if isempty(done); done = 0; end
    done = done + 1;
    elapsed = toc(startTime);
    pctDone = done/total;
    eta = elapsed/pctDone - elapsed;
    fprintf("Progress: %d/%d (%.1f%%), ETA %.0f s\n", done, total, 100*pctDone, eta);
end
