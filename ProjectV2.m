clearvars; close all; clc;

%% Setup

% Which simulation to run:
% paper - replicate paper sim
% ecc - test different eccentricities
% Mol - simulate Molniya orbit
% optKw - find optimal kw for scenario A
% NearA - optimal value vs analytic at scenarios like A
% attitudeSweep - optimal vs analytic kw over a grid of initial attitudes
% alphaSweep - optimal vs analytic over varying initial alpha
sim_type = "cmp";
parallels = false;  % Run on multiple cores
conv_tol = 0.01; % Convergence tolarence.
enabledPlots = struct( ...      % initilize which plots to plot
    'angularRate', false, ...
    'controlTorque', false, ...
    'quaternionNorm', false, ...
    'actuatorDipoles', false, ...
    'controlAuthority', false, ...
    'settlingTimeVsKw', false, ...
    'nearAComparison', false, ...
    'attitudeSweepByTrial', false, ...
    'attitudeSweepByAngles', false, ...
    'controlAngleSweep', false);
% saveName = struct( ...      % initilize which plots to save
%     'angularRate', "", ...
%     'quaternionNorm', "", ...
%     'actuatorDipoles', "", ...
%     'controlAuthority', "", ...
%     'settlingTimeVsKw', "", ...
%     'nearAComparison', "", ...
%     'attitudeSweepByTrial', "", ...
%     'attitudeSweepByAngles', "", ...
%     'controlAngleSweep', "");
saveName = struct( ...      % initilize which plots to save
    'angularRate', "pt3_MolARcmp.png", ...
    'controlTorque', "pt3_MolCTcmp.png", ...
    'quaternionNorm', "", ...
    'actuatorDipoles', "", ...
    'controlAuthority', "pt3_MolCAcmp.png", ...
    'settlingTimeVsKw', "", ...
    'nearAComparison', "", ...
    'attitudeSweepByTrial', "", ...
    'attitudeSweepByAngles', "", ...
    'controlAngleSweep', "");

J1 = 0.33; J2 = 0.37; J3 = 0.35; % [kg m^2] Principal moments of inertia
I = diag([J1, J2, J3]);

mlim = 2; % [A m^2] acuator saturation limit
r_c = 7021; % [km] orbit radius
To = 5855; % [s] orbit period
gamma_m = 11.44; % [deg] geomagnetic plane tilt
beta_m = 0; % [deg] initial geomagnetic phase
mu = 3.986e5; % [km^3/s^2] earth gravitational parameter
Re = 6378; % [km] earth radius
N = 4; % [-] number of orbits
l = 1; % which trial to use for actuator dipoles
kwopt_A = 8.340986872382716e-04; % Optimal kw for case A found numerically.
kwopt_B = 0.001145973154362;
kwopt_Mol = 0.001164433622563;


a = (mu*To^2/(4*pi^2))^(1/3); % [km] semi major axis

switch sim_type
    case 'paper'
        % Cases A, B, and C
        i = [11, 65, 65]; % [deg] inclination
        trials = 3;
        omega0 = [0.604 -0.76 -0.384;
                  0.604 -0.76 -0.384;
                  1      0     0].'; % [rad/s] initial angular velocity b to I (body frame)
        Q0 = [0.375 -0.062  0.925 -0.007;
              0.375 -0.062  0.925 -0.007;
              0.646  0.525 -0.514  0.206].'; % [-] initial quaternions (LVLH to body)
        Q0 = Q0./vecnorm(Q0, 2, 1);
        
        Omega0 = [0, 272, 272]; % [deg] inital ascending node
        arglat0 = [0, 0, 36]; % [deg] initial argument of latitude
        w0 = zeros(size(i)); % [deg] initial argument of perigee
        ecc = zeros(size(i)); % [-] eccentricity
        leg_text = struct("omega", ['||\omega_A||';'||\omega_B||';'||\omega_C||'],...
        "Q", ['||q_A||';'||q_B||';'||q_C||'],...
        "alpha", ['\alpha_A';'\alpha_B';'\alpha_C']);
        
        kw = analyticKw(i, Omega0, gamma_m, beta_m, To, I);
        IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Q0, omega0, kw, trials, To);
        enabledPlots.angularRate = true;
        enabledPlots.quaternionNorm = true;
        % enabledPlots.actuatorDipoles = true; l = 1;
        enabledPlots.controlAuthority = true;

     case 'cmp'
        % Compare two cases 
        i = 62.8; % [deg] inclination
        trials = 2;
        omega0 = [0.604 -0.76 -0.384].'; % [rad/s] initial angular velocity b to I (body frame)
        Q0 = [1 0 0 0].'; % [-] initial quaternions (LVLH to body)
        Q0 = Q0./vecnorm(Q0, 2, 1);
        
        Omega0 = 0; % [deg] inital ascending node
        arglat0 = 270; % [deg] initial argument of latitude
        w0 = 270; % [deg] initial argument of perigee
        ecc = 0.7411; % [-] eccentricity
        leg_text = struct("omega", ["||\omega_{analytic}||";"||\omega_{optimal}||"],...
        "Q", ["||q_{analytic}||";"||q_{optimal}||"],...
        "alpha", ["\alpha_{analytic}";"\alpha_{optimal}"],...
        "M", ["||M_{analytic}||";"||M_{optimal}||"]);
        
        a = 27015; Torb = 2*pi*sqrt(a.^3/mu); To = Torb;
        kw = analyticKw(i, Omega0, gamma_m, beta_m, To, I);
        IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Q0, omega0, [kw kwopt_Mol], trials, To);
        enabledPlots.angularRate = true;
        enabledPlots.controlTorque = true;
        % enabledPlots.quaternionNorm = true;
        % enabledPlots.actuatorDipoles = true; l = 1;
        enabledPlots.controlAuthority = true;

    case 'ecc'  
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
        
        
        as = r_c./sqrt(1-ecc.^2); Torbs = 2*pi*sqrt(as.^3/mu);  % Normalize average magnetic field strength
        kw = analyticKw(i, Omega0, gamma_m, beta_m, Torbs, I);
        IC = buildIC(as, ecc, i, Omega0, arglat0, w0, Q0, omega0, kw, trials, Torbs);
        
        enabledPlots.angularRate = true;
        enabledPlots.quaternionNorm = true;
        % enabledPlots.actuatorDipoles = true; l = 1;
        enabledPlots.controlAuthority = true;

    case 'Mol'  
        ecc = 0.7411; % [-] eccentricity
        trials = 1;
        i = 62.8; % [deg] inclination
        omega0 = [0.604;-0.76;-0.384]*ones(size(ecc)); % [rad/s] initial angular velocity b to I (body frame)
        Q0 = [1;0;0;0]*ones(size(ecc)); % [-] initial quaternions (LVLH to body)
        Q0 = Q0./vecnorm(Q0, 2, 1);
        
        Omega0 = 0; % [deg] inital ascending node
        arglat0 = 270; % [deg] initial argument of latitude
        w0 = 270; % [deg] initial argument of perigee
        
        a = 27015; Torb = 2*pi*sqrt(a.^3/mu); To = Torb;
        kw = analyticKw(i, Omega0, gamma_m, beta_m, Torb, I);
        % kw = kwopt_Mol;
        IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Q0, omega0, kw, trials, Torb);
        
        enabledPlots.angularRate = true;
        enabledPlots.quaternionNorm = true;
        enabledPlots.actuatorDipoles = true; l = 1;
        enabledPlots.controlAuthority = true;

        leg_text = struct();

    case 'optKw'
        trials = 80;
        parallels = true;
        i = 62.8; % [deg] inclination
        omega0 = [0.604 -0.76 -0.384].'; % [rad/s] initial angular velocity b to I (body frame)
        Q0 = [1 0 0 0].'; % [-] initial quaternions (LVLH to body)
        Q0 = Q0./vecnorm(Q0, 2, 1);
        
        Omega0 = 0; % [deg] inital ascending node
        arglat0 = 270; % [deg] initial argument of latitude
        w0 = 270; % [deg] initial argument of perigee
        ecc = 0.7411; % [-] eccentricity
        
        a = 27015; Torb = 2*pi*sqrt(a.^3/mu); To = Torb;
        kw_ref = analyticKw(i, Omega0, gamma_m, beta_m, To, I);
        kw = [kw_ref, linspace(0.9*kw_ref, 10*kw_ref, trials-1)];
        IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Q0, omega0, kw, trials, Torb);
        leg_text = struct();

        enabledPlots.settlingTimeVsKw = true;

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
        IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Q0, omega0, kw, trials, To);
        leg_text = struct("omega", compose('k_\\omega=%e',kw),...
        "Q", compose('k_\\omega=%e',kw),...
        "alpha", compose('k_\\omega=%e',kw));

        enabledPlots.angularRate = true;
        enabledPlots.quaternionNorm = true;
        % enabledPlots.actuatorDipoles = true; l = 1;
        enabledPlots.controlAuthority = true;
        enabledPlots.nearAComparison = true;

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
        
        i = 65;
        Omega0 = 0;
        arglat0 = 0;
        w0 = 0;
        ecc = 0;
        
        omega0v = [0.604; -0.76; -0.384];
        
        kwAnalytic = analyticKw(i, Omega0, gamma_m, beta_m, To, I);
        
        kw = [kwAnalytic*ones(1,trials/2), kwopt_B*ones(1,trials/2)];
        IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Qsample, omega0v, kw, trials, To);
        
        leg_text = struct("omega", [], "Q", [], "alpha", []);

        enabledPlots.attitudeSweepByTrial = true;
        enabledPlots.attitudeSweepByAngles = true;
        enabledPlots.controlAngleSweep = true;

    case "alphaSweep"
        
        q = [0.375; -0.062; 0.925; -0.007];
        q = q / norm(q);
        omega0v = [0.604; -0.76; -0.384]; 
        
        ecc     = 0.7411;   % [-]
        i       = 62.8;     % [deg]
        Omega0  = 0;        % [deg]
        arglat0 = 270;      % [deg]
        w0      = 270;      % [deg]
        a       = 27015;    % [km]
        Torb    = 2*pi*sqrt(a^3/mu);
        To      = Torb;
        
        parallels = true;
        
        [r0ref, v0ref] = Kparams2ECI(struct("a",a, "e",ecc, "i",i, ...
            "omega",w0, "Omega",Omega0, "f",arglat0-w0));
        [~, bB_ref] = dynamics(0, q(1), q(2:4), omega0v, r0ref, v0ref, I=eye(3), kw=0);
        alpha_ref = acosd( dot(omega0v, bB_ref) / (norm(omega0v)*norm(bB_ref)) );
        n_hat = cross(omega0v, bB_ref); n_hat = n_hat / norm(n_hat);
        
        trialsPerKw = 7;
        alpha0 = linspace(1, 179, trialsPerKw); % [deg]
        
        dtheta = alpha0 - alpha_ref;                   
        dQ = [cosd(dtheta/2); n_hat .* sind(dtheta/2)]; 
        Qalpha = quatmultiply(dQ.', q.').';
        Qalpha = Qalpha ./ vecnorm(Qalpha, 2, 1);
        

        kwAnalytic = analyticKw(i, Omega0, gamma_m, beta_m, To, I); 
        
        Q0     = [q, Qalpha, q Qalpha];
        omega0 = [omega0v, omega0v*ones(1,trialsPerKw), omega0v, omega0v*ones(1,trialsPerKw)];
        kw     = [kwAnalytic, kwAnalytic*ones(1,trialsPerKw), kwopt_Mol, kwopt_Mol*ones(1,trialsPerKw)];
        trials = length(kw);
        
        IC = buildIC(a, ecc, i, Omega0, arglat0, w0, Q0, omega0, kw, trials, To);
        
        leg_text = struct();
        
        enabledPlots.attitudeSweepByTrial = true;
        enabledPlots.controlAngleSweep = true;

end

%% Simulation

convTime = zeros(1, trials);
results = cell(1,trials);
convIdx = zeros(1, trials);

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
        convTime(j) = interp1(vecnorm(omega,2,2), Time, conv_tol*norm(omega(1,:)))/To;
        send(D, 1);
    end

else
    opts = odeset(Refine=1, Stats='on');
    for j=1:trials
        fprintf("Trial number %d\n",j);
        ICj = IC(j);
        
        fun = @(t,x) dynamics(t, x(1), x(2:4), x(5:7), ICj.r0, ICj.v0, I=I, kw=ICj.kw);
        x0 = [ICj.Q0; ICj.omega0];
        
        [Time, x] = ode89(fun, [0 N*ICj.Torb], x0, opts);
        Q = x(:,1:4); omega = x(:,5:end);
        [rates,b,m,M] = dynamics(Time, Q(:,1), Q(:,2:4).', omega.', ICj.r0, ICj.v0, I=I, kw=ICj.kw);
        b = b.'; m = m.'; M = M.';
        results(j) = {table(Time, Q, omega, b, m, M)};
        lastAbove = find(vecnorm(omega,2,2) >= conv_tol*norm(omega(1,:)), 1, 'last');
        if lastAbove == length(Time)
            fprintf("Didn't converge in %d orbits\n", N)
        end
        convIdx(j) = min(lastAbove+1,length(Time));
        convTime(j) = interp1(vecnorm(omega,2,2), Time, conv_tol*norm(omega(1,:)))/ICj.Torb;
    end

end

%% Plots

plots = plotting();

ctx = struct('results',{results}, 'To',To, 'leg_text',leg_text, ...
    'mlim',mlim, 'l',l, 'convIdx',convIdx, 'kw',kw, ...
    'convTime',convTime, 'trials',trials);

if strcmp(sim_type, "attitudeSweep")
    ctx.axisColat = axisColat; ctx.axisLon = axisLon; ctx.rotAngle = rotAngle;
    ctx.axisColatGrid = axisColatGrid; ctx.axisLonGrid = axisLonGrid;
    ctx.rotAngleGrid = rotAngleGrid;
    ctx.kwOptimal = kwopt_B; ctx.kwAnalytic = kwAnalytic;
    ctx.improvementPct = 100*(convTime(1:trials/2)-convTime(trials/2+1:end)) ...
                          ./ convTime(1:trials/2);
    ctx.alpha0 = initialControlAngle(IC(1:trials/2));
end

if strcmp(sim_type, "alphaSweep")
    ctx.alpha0 = [alpha_ref, alpha_ref, alpha0, alpha0];
    ctx.kwOptimal  = kwopt_Mol;
    ctx.kwAnalytic = kwAnalytic;
    ctx.baselineIdx = [1, 2];
    ctx.improvementPct = 100*(convTime(1:trials/2)-convTime(trials/2+1:end)) ...
                          ./ convTime(1:trials/2);
end

runPlots(plots, enabledPlots, ctx, saveName);


%% Local functions

function IC = buildIC(a, ecc, incl, Omega0, arglat0, w0, Q0, omega0, kw, trials, Torb)
%BUILDIC Generate IC struct

    arguments
        a       (1,:) double {mustBePositive}
        ecc     (1,:) double {mustBeNonnegative}
        incl    (1,:) double
        Omega0  (1,:) double
        arglat0 (1,:) double
        w0      (1,:) double
        Q0      (4,:) double
        omega0  (3,:) double
        kw      (1,:) double {mustBePositive}
        trials  (1,1) {mustBeInteger, mustBePositive}
        Torb    (1,:) double {mustBePositive}
    end
    
    a       = broadcastCols(a, trials);
    ecc     = broadcastCols(ecc, trials);
    incl    = broadcastCols(incl, trials);
    Omega0  = broadcastCols(Omega0, trials);
    arglat0 = broadcastCols(arglat0, trials);
    w0      = broadcastCols(w0, trials);
    Q0      = broadcastCols(Q0, trials);
    omega0  = broadcastCols(omega0, trials);
    kw      = broadcastCols(kw, trials);
    Torb    = broadcastCols(Torb, trials);
    
    IC(trials) = struct('r0',[],'v0',[],'Q0',[],'omega0',[],'kw',[], 'Torb', []);
    for j = 1:trials
        [r0, v0] = Kparams2ECI(struct("a",a(j), "e",ecc(j), "i",incl(j), ...
            "omega",w0(j), "Omega",Omega0(j), "f",arglat0(j)-w0(j)));
        IC(j).r0 = r0;
        IC(j).v0 = v0;
        IC(j).Q0 = Q0(:,j);
        IC(j).omega0 = omega0(:,j);
        IC(j).kw = kw(j);
        IC(j).Torb = Torb(j);
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

function kw = analyticKw(i, Omega0, gamma_m, beta_m, To, I)
%ANALYTICKW vectorized implementation of the analytical Kw from paper.

    arguments
        i       (1,:) double
        Omega0  (1,:) double
        gamma_m (1,1) double
        beta_m  (1,1) double
        To      (1,:) double {mustBePositive}
        I       (3,3) double
    end

    xi_m = acosd( cosd(i).*cosd(gamma_m) + sind(i).*sind(gamma_m).*cosd(beta_m-Omega0) );
    kw = 4*pi./To .* (1+sind(xi_m)) * min(diag(I));
end

function alpha0 = initialControlAngle(IC)
%INITIALCONTROLANGLE Calculate the inital angle between angular velocity and
% magnetic field.
    arguments
        IC (1,:) struct
    end
    alpha0 = zeros(size(IC));
    for j=1:length(IC)
        ICj = IC(j);
        [~, b_B] = dynamics(0, ICj.Q0(1), ICj.Q0(2:4), ICj.omega0, ICj.r0, ICj.v0, ...
                             I=eye(3), kw=0);   % Reusing the function only for b so many parameters are irrelevant
        cosT = dot(ICj.omega0, b_B) / (norm(ICj.omega0)*norm(b_B));
        alpha0(j) = acosd(min(max(cosT,-1),1));
    end
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

function reportProgressBis(data, total, startTime)
%REPORTPROGRESSBis Helper function for parfor loop progress tracking in bisection loop.
    arguments
        data      (1,1) struct
        total     (1,1) double {mustBePositive, mustBeInteger}
        startTime (1,1) uint64
    end
    persistent done
    if isempty(done); done = 0; end
    
    if data.iter < 0
        done = done + 1;
        elapsed = toc(startTime);
        pctDone = done/total;
        eta = elapsed/pctDone - elapsed;
        fprintf("Orbit %d Converged\n", data.id)
        fprintf("Progress: %d/%d (%.1f%%), ETA %.0f s\n", done, total, 100*pctDone, eta);
    else
        fprintf("Orbit %d finished iteration %d\n", data.id, data.iter)
    end
    
end
