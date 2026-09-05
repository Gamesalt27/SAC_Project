clearvars; close all; clc;

%% Setup

% Which simulation to run: 
% paper - replicate paper sim 
% ecc - test different eccentricities
% optKw - find optimal kw for scenario A
% NearA - optimal value vs analytic at scenarios like A
sim_type = "attitudeSweep";
parallels = false;       % Log entire trial
conv_tol = 0.01;  % Convergence tolarence.

J1 = 0.33; J2 = 0.37; J3 = 0.35;    % [kg m^2] Principal moments of inertia
I = diag([J1, J2, J3]);

mlim    = 2;        % [A m^2] acuator saturation limit
r_c     = 7021;     % [km] orbit radius
To      = 5855;     % [s] orbit period
gamma_m = 11.44;    % [deg] geomagnetic plane tilt
beta_m  = 0;        % [deg] initial geomagnetic phase
mu      = 3.986e5;  % [km^3/s^2] earth gravitational parameter
Re      = 6378;     % [km] earth radius
N       = 8;        % [-] number of orbits
n       = 1e3;      % [-] points per orbit

a = (mu*To^2/(4*pi^2))^(1/3);  % [km] semi major axis

switch sim_type
    case 'paper'
        % Cases A, B, and C
        i      = [11, 65, 65];                     % [deg] inclination
        trials = length(i);
        omega0 = [0.604 -0.76 -0.384; 
                  0.604 -0.76 -0.384; 
                  1      0     0].';               % [rad/s] initial angular velocity b to I (body frame)
        Q0     = [0.375 -0.062  0.925 -0.007; 
                  0.375 -0.062  0.925 -0.007; 
                  0.646  0.525 -0.514  0.206].';   % [-] initial quaternions (LVLH to body)
        Q0 = Q0./vecnorm(Q0, 2, 1);
        
        Omega0  = [30, 272, 272];   % [deg] inital ascending node
        arglat0 = [0, 36, 36];      % [deg] initial argument of latitude
        w0      = zeros(size(i));   % [deg] initial argument of perigee
        ecc     = zeros(size(i));   % [-] eccentricity
        leg_text = struct("omega", ['||\omega_A||';'||\omega_B||';'||\omega_C||'],...
            "Q", ['||q_A||';'||q_B||';'||q_C||'],...
            "alpha", ['\theta_A';'\theta_B';'\theta_C']);

        IC(trials) = struct('r0',[],'v0',[],'Q0',[],'omega0',[],'kw',[]);
        for j = 1:trials
            [r0,v0] = Kparams2ECI(struct("a",a,"e",ecc(j),"i",i(j),"omega",w0(j), ...
                                            "Omega",Omega0(j),"f",arglat0(j)-w0(j)));
            xi_m = acosd( cosd(i(j))*cosd(gamma_m) + sind(i(j))*sind(gamma_m)*cosd(beta_m-Omega0(j)) );
            IC(j).r0 = r0; IC(j).v0 = v0;
            IC(j).Q0 = Q0(:,j); IC(j).omega0 = omega0(:,j);
            IC(j).kw = 4*pi/To * (1+sind(xi_m)) * min(diag(I));
        end

    case 'ecc'
        ecc    = [0 0.1 0.5 0.85 0.99];                        % [-] eccentricity 
        trials = length(ecc);
        i      = ones(size(ecc))*65;                           % [deg] inclination
        omega0 = [0.604;-0.76;-0.384]*ones(size(ecc));         % [rad/s] initial angular velocity b to I (body frame)
        Q0     = [0.375;-0.062;0.925;-0.007]*ones(size(ecc));  % [-] initial quaternions (LVLH to body)
        Q0 = Q0./vecnorm(Q0, 2, 1);
        
        Omega0  = zeros(size(ecc));   % [deg] inital ascending node
        arglat0 = zeros(size(ecc));   % [deg] initial argument of latitude
        w0      = zeros(size(ecc));   % [deg] initial argument of perigee 
        leg_text = struct("omega", compose('e=%0.2f',ecc),...
            "Q", compose('e=%0.2f',ecc),...
            "alpha", compose('e=%0.2f',ecc));
        IC(trials) = struct('r0',[],'v0',[],'Q0',[],'omega0',[],'kw',[]);
        for j = 1:trials
            [r0,v0] = Kparams2ECI(struct("a",a,"e",ecc(j),"i",i(j),"omega",w0(j), ...
                                            "Omega",Omega0(j),"f",arglat0(j)-w0(j)));
            xi_m = acosd( cosd(i(j))*cosd(gamma_m) + sind(i(j))*sind(gamma_m)*cosd(beta_m-Omega0(j)) );
            IC(j).r0 = r0; IC(j).v0 = v0;
            IC(j).Q0 = Q0(:,j); IC(j).omega0 = omega0(:,j);
            IC(j).kw = 4*pi/To * (1+sind(xi_m)) * min(diag(I));
        end

    case 'optKw'
        trials = 200;
        parallels = true;
        i      = 11;                               % [deg] inclination
        omega0 = [0.604 -0.76 -0.384].';           % [rad/s] initial angular velocity b to I (body frame)
        Q0     = [0.375 -0.062  0.925 -0.007].';   % [-] initial quaternions (LVLH to body)
        Q0 = Q0./vecnorm(Q0, 2, 1);
        
        Omega0  = 0;   % [deg] inital ascending node
        arglat0 = 0;   % [deg] initial argument of latitude
        w0      = 0;   % [deg] initial argument of perigee
        ecc     = 0;   % [-] eccentricity
        
   
        [r0,v0] = Kparams2ECI(struct("a",a,"e",ecc,"i",i,"omega",w0,"Omega", Omega0, "f", arglat0-w0));
        xi_m = acosd( cosd(i)*cosd(gamma_m) + sind(i)*sind(gamma_m)*cosd(beta_m-Omega0));
        kw_ref = 4*pi/To * (1+sind(xi_m)) * min(diag(I));
        x0 = [Q0; omega0].';
        kws = [kw_ref, linspace(0.9*kw_ref,1.3*kw_ref,trials-1)];
        IC = struct('r0',ones(1,trials).*r0,'v0',ones(1,trials).*v0,...
            'Q0',ones(1,trials).*Q0,'omega0',ones(1,trials).*omega0,'kw',kws);
        leg_text = struct("omega", compose('k_\\omega=%e',kws),...
            "Q", compose('k_\\omega=%e',kws),...
            "alpha", compose('k_\\omega=%e',kws));

    case "NearA"
        trials = 8;
        i      = 11;                               % [deg] inclination
        omega0 = [0.604 -0.76  -0.384;
                  0.604 -0.76  -0.384;
                 -0.76   0.604 -0.384;
                  0.604 -0.76  -0.384].';           % [rad/s] initial angular velocity b to I (body frame)
        omega0 = [omega0,omega0];
        Q0     = [0.375 -0.062  0.925 -0.007;
                  0.646  0.525 -0.514  0.206;
                  0.375 -0.062  0.925 -0.007;
                  0.375 -0.062  0.925 -0.007].';   % [-] initial quaternions (LVLH to body)
        Q0 = Q0./vecnorm(Q0, 2, 1);
        Q0 = [Q0, Q0];
        
        Omega0  = 0;   % [deg] inital ascending node
        arglat0 = [0 0 0 -10];   % [deg] initial argument of latitude
        arglat0 = [arglat0,arglat0];
        w0      = 0;   % [deg] initial argument of perigee
        ecc     = 0;   % [-] eccentricity

        xi_m = acosd( cosd(i)*cosd(gamma_m) + sind(i)*sind(gamma_m)*cosd(beta_m-Omega0) );
        kw = [4*pi/To * (1+sind(xi_m)) * min(diag(I))*ones(1,trials/2),...
            8.340986872382716e-04*ones(1,trials/2)];    % Optimal Kw for scenario A        
   
        IC = struct('r0',zeros(size(omega0)),'v0',zeros(size(omega0)),'Q0',Q0,'omega0',omega0,'kw',kw);
        for j = 1:trials
            [r0,v0] = Kparams2ECI(struct("a",a,"e",ecc,"i",i,"omega",w0, ...
                                            "Omega",Omega0,"f",arglat0(j)-w0));
            IC.r0(:,j) = r0; IC.v0(:,j) = v0;
        end
        leg_text = struct("omega", compose('k_\\omega=%e',kw),...
            "Q", compose('k_\\omega=%e',kw),...
            "alpha", compose('k_\\omega=%e',kw));

    case "attitudeSweep"
        qA = [0.375; -0.062; 0.925; -0.007];
        qA = qA / norm(qA);

        % Generate equally spaced quaternions 
        axisColat = 15:15:165;  % Angle between rotation vectors
        axisLon   = 0:30:330;  % Complementry angle
        rotAngle  = 15:15:165;  % Rotation angle

        [axisColatGrid,axisLonGrid,rotAngleGrid] = meshgrid(axisColat,axisLon,rotAngle);
        axisColatGrid = reshape(axisColatGrid,1,[]);
        axisLonGrid   = reshape(axisLonGrid,1,[]);
        rotAngleGrid  = reshape(rotAngleGrid,1,[]);

        % Add pure rotations at the beginning
        axisColatGrid = [zeros(size(rotAngle)), axisColatGrid];     
        axisLonGrid   = [zeros(size(rotAngle)), axisLonGrid]; 
        rotAngleGrid  = [rotAngle, rotAngleGrid];
        axisColat     = [0, axisColat];


        u = [sind(axisColatGrid).*cosd(axisLonGrid).*sind(rotAngleGrid/2);...
             sind(axisColatGrid).*sind(axisLonGrid).*sind(rotAngleGrid/2);...
             cosd(axisColatGrid).*sind(rotAngleGrid/2)];
        dQ = [cosd(rotAngleGrid/2);u];   % Difference Quaternion
        Qsample = quatmultiply(dQ.',qA.').';
        Qsample = [Qsample, Qsample];

        trials = length(Qsample);             
        parallels = true;      
    
        i       = 11;   
        Omega0  = 0;
        arglat0 = 0;
        w0      = 0;
        ecc     = 0;
    
        omega0v = [0.604; -0.76; -0.384];
    
        [r0,v0] = Kparams2ECI(struct("a",a, "e",ecc, "i",i, "omega",w0, ...
            "Omega",Omega0, "f",arglat0-w0));
    
        xi_m = acosd(cosd(i)*cosd(gamma_m) + sind(i)*sind(gamma_m)*cosd(beta_m-Omega0));
        kwAnalytic = 4*pi/To * (1+sind(xi_m)) * min(diag(I));
        kwOptimal  = 8.340986872382716e-4;  % Found numerically
    
        IC = struct( ...
            'r0',     repmat(r0, 1, trials), ...
            'v0',     repmat(v0, 1, trials), ...
            'Q0',     Qsample, ...
            'omega0', repmat(omega0v, 1, trials), ...
            'kw',     [kwAnalytic*ones(1,trials/2), ...
                        kwOptimal *ones(1,trials/2)]);
    
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
        fun = @(t,x) dynamics(t, x(1), x(2:4), x(5:7), IC.r0(:,j), IC.v0(:,j), I=I, kw=IC.kw(j));
        x0 = [IC.Q0(:,j); IC.omega0(:,j)];
    
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
    
        fun = @(t,x) dynamics(t, x(1), x(2:4), x(5:7), IC.r0(:,j), IC.v0(:,j), I=I, kw=IC.kw(j));
        x0 = [IC.Q0(:,j); IC.omega0(:,j)];
    
        [Time, x] = ode89(fun, tspan, x0, opts);
        Q = x(:,1:4); omega = x(:,5:end); 
        [rates,b,m,M] = dynamics(Time, Q(:,1), Q(:,2:4).', omega.', IC.r0(:,j), IC.v0(:,j), I=I, kw=IC.kw(j));
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

function reportProgress(total, startTime)
    persistent done
    if isempty(done); done = 0; end
    done = done + 1;
    elapsed = toc(startTime);
    pctDone = done/total;
    eta = elapsed/pctDone - elapsed;
    fprintf("Progress: %d/%d (%.1f%%), ETA %.0f s\n", done, total, 100*pctDone, eta);
end


%% angular rate
handle = figure('Position', [100, 100, 600, 400]);
hold on; grid on;
for j=1:trials
    plot(results{j}.Time/To, vecnorm(results{j}.omega,2,2), LineWidth=1.5);
end
hold off; title("Angular rate over time"); xlabel('t/T [-]');
ylabel('||\omega|| [rad/s]'); legend(leg_text.omega, Interpreter="tex", Location="best")
% exportgraphics(handle, "graphs/pt2_omega.png", "Resolution", 300);

%% Quaternion size (verify sim accuracy)
handle = figure('Position', [100, 100, 600, 400]);
hold on; grid on;
for j=1:trials
    plot(results{j}.Time/To, vecnorm(results{j}.Q,2,2), LineWidth=1.5);
end
hold off; title("Quaternion over time"); xlabel('t/T [-]');
ylabel('||q||'); legend(leg_text.Q, Interpreter="tex", Location="best")
% exportgraphics(handle, "graphs/pt2_Q.png", "Resolution", 300);

%% Actuators
l = 629;  % Which sim to plot
handle = figure('Position', [100, 100, 600, 400]);
tl = tiledlayout(3, 1, TileSpacing='compact', Padding='compact');
axLabels = {'m_1', 'm_2', 'm_3'};

for k=1:3
    nexttile;
    hold on; grid on;
    plot(results{l}.Time/To, results{l}.m(:,k), 'b', 'LineWidth', 1.5);
    yline( mlim, 'r--', 'LineWidth', 1.2);
    yline(-mlim, 'r--', 'LineWidth', 1.2);
    title(axLabels{k});
    ylabel([axLabels{k} ' [A m^2]']);
    ylim([-1.15*mlim, 1.15*mlim]);
    hold off;
end
xlabel('t/T [-]');

title(tl, 'Magnetic Actuator Dipole Moments');
lgd = legend({'m_i', 'Saturation limit'}, Orientation='horizontal', Interpreter='tex');
lgd.Layout.Tile = 'north';

% exportgraphics(handle, "graphs/pt2_actu.png", "Resolution", 300);

%% Control authority

handle = figure('Position', [100, 100, 600, 400]);
hold on; grid on;

for j=i:trials
    omega_j = results{j}.omega;  
    b_j     = results{j}.b;       

    cosAlpha = dot(omega_j, b_j, 2) ./ (vecnorm(omega_j,2,2).*vecnorm(b_j,2,2));
    cosAlpha = min(max(cosAlpha, -1), 1);   % prevent numerical errors leading to complex result
    alpha = acosd(cosAlpha);               

    plot(results{j}.Time(1:convIdx(j))/To, alpha(1:convIdx(j)), 'LineWidth', 1.5);
end

yline(90, 'k--', 'LineWidth', 1);
hold off;
title('Angle Between \omega and b');
xlabel('t/T [-]');
ylabel('\alpha(\omega, b) [deg]');
ylim([0, 180]);
legend(leg_text.alpha, Location='best', Interpreter='tex');

%% settling time to kw

if strcmp(sim_type, "optKw")
    [Tmin, idx] = min(convTime);
    
    handle = figure('Position', [100, 100, 600, 400]);
    hold on; grid on;
    plot(IC.kw(2:end), convTime(2:end), 'b', LineWidth=1.5, HandleVisibility='off');
    plot(IC.kw(idx), Tmin, 'g*', DisplayName="Optimal");
    plot(IC.kw(1), convTime(1), 'k*', DisplayName="Analytical")
    hold off; title('Settling time vs k_\omega'); ylabel('t/T');
    xlabel('k_\omega'); legend("show", Interpreter="tex", Location="best")
    % exportgraphics(handle, "graphs/pt2_omega.png", "Resolution", 300);
end

%% Near A comparison

nScenarios = trials/2;

analyticalTime = convTime(1:nScenarios);
optimalTime    = convTime(nScenarios+1:end);

analyticalKw = IC.kw(1);
optimalKw    = IC.kw(nScenarios+1);

pctImprovement = 100*(analyticalTime - optimalTime)./analyticalTime;   
scenarioLabels = {'Control', 'Initial Q', 'Initial $\omega$', 'Initial \theta'};

handle = figure('Position', [100, 100, 700, 450]);
hold on; grid on;

barColors = repmat([0.2 0.7 0.3], nScenarios, 1);  
barColors(pctImprovement < 0, :) = repmat([0.8 0.2 0.2], sum(pctImprovement<0), 1);

b = bar(categorical(scenarioLabels, scenarioLabels), pctImprovement, FaceColor='flat', Interpreter="tex");
b.CData = barColors;
yline(0, 'k-', 'LineWidth', 1);

% annotate each bar with its exact numeric value, since bars are short
for k = 1:nScenarios
    va = 'bottom'; offset = 0.05*max(abs(pctImprovement));
    if pctImprovement(k) < 0, va = 'top'; offset = -offset; end
    text(k, pctImprovement(k) + offset, sprintf('%.3f%%', pctImprovement(k)), ...
         'HorizontalAlignment','center', 'FontSize', 9);
end

title(sprintf('Settling-Time Improvement: Optimal vs. Analytical k_\\omega (k_\\omega=%.3e vs %.3e)', ...
               optimalKw, analyticalKw));
ylabel('Improvement [%]');
ylim([min(0,1.4*min(pctImprovement)), 1.4*max(pctImprovement)]);

% exportgraphics(handle, "graphs/pt2_NearA_improvement.png", "Resolution", 300);

%% attitude sweep by trial number

tAnalytic = convTime(1:trials/2);
tOptimal  = convTime(trials/2+1:end);

improvementPct = 100*(tAnalytic - tOptimal)./tAnalytic;

figure('Position',[100 100 750 400]);
hold on; grid on;

stem(1:trials/2, improvementPct, 'filled', LineWidth=1.1, MarkerSize=4, Color='g');

yline(0, 'k-', 'LineWidth', 1.2);
xlabel('Initial-attitude sample');
ylabel('Settling-time change [%]');
title('Numerical k_\omega Improvement over Analytical k_\omega');

%% attitude sweep by angles
colatMed = zeros(size(axisColat));
colatStd = zeros(size(axisColat));
lonMed = zeros(size(axisLon));
lonStd = zeros(size(axisLon));
rotMed = zeros(size(rotAngle));
rotStd = zeros(size(rotAngle));

for j=1:length(axisColat)
    colatMed(j) = median(improvementPct(axisColatGrid == axisColat(j)));
    colatStd(j) = std(improvementPct(axisColatGrid == axisColat(j)));
end

for j=1:length(axisLon)
    lonMed(j) = median(improvementPct(axisLonGrid == axisLon(j)));
    lonStd(j) = std(improvementPct(axisLonGrid == axisLon(j)));
end

for j=1:length(rotAngle)
    rotMed(j) = median(improvementPct(rotAngleGrid == rotAngle(j)));
    rotStd(j) = std(improvementPct(rotAngleGrid == rotAngle(j)));
end

kwInfo = sprintf(' (k_\\omega=%.3e vs %.3e)', kwOptimal, kwAnalytic);

plotAngleMarginalBar(axisColat, colatMed, colatStd, 'deg', ...
                     'Axis Colatitude \chi', kwInfo);

plotAngleMarginalBar(axisLon, lonMed, lonStd, 'deg', ...
                     'Axis Longitude \lambda', kwInfo);

plotAngleMarginalBar(rotAngle, rotMed, rotStd, 'deg', ...
                     'Rotation Angle \alpha', kwInfo);

% exportgraphics(gcf, "graphs/pt2_attitudeSweep_rotAngle.png", "Resolution", 300);

%% Marginal median/spread bar plots for the attitude sweep
% One figure per swept angle (axis colatitude, axis longitude, rotation
% angle). Each bar shows the median improvement for that angle value,
% with an error bar showing +/- 3*std.

function plotAngleMarginalBar(angleValues, medVals, stdVals, angleUnit, angleName, kwInfo)
%PLOTANGLEMARGINALBAR One median+spread bar chart for a swept angle.
% angleValues : the grid values for this angle (e.g. axisColat)
% medVals     : median(improvementPct) at each angleValues(k)
% varVals     : var(improvementPct) at each angleValues(k)
% angleUnit   : string, e.g. 'deg'
% angleName   : string used in title/labels, e.g. 'Axis Colatitude \chi'
% kwInfo      : string appended to the title (e.g. kw values used)

    nBars = numel(angleValues);
    spread3sigma = 3*stdVals;

    labels = compose('%g %s', angleValues(:), angleUnit);

    barColors = repmat([0.2 0.7 0.3], nBars, 1);
    barColors(medVals < 0, :) = repmat([0.8 0.2 0.2], sum(medVals<0), 1);

    figure('Position', [100, 100, 700, 450]);
    hold on; grid on;

    b = bar(categorical(labels, labels), medVals, FaceColor='flat');
    b.CData = barColors;

    errorbar(1:nBars, medVals, spread3sigma, 'k', LineStyle='none', ...
             LineWidth=1.2, CapSize=8, HandleVisibility='off');

    yline(0, 'k-', 'LineWidth', 1);

    % annotate each bar with its exact median value
    yPad = 0.08*max(abs(medVals) + spread3sigma);
    for k = 1:nBars
        offset = yPad;
        if medVals(k) < 0, offset = -offset; end
        text(k, medVals(k) + spread3sigma(k)*sign(offset) + offset, ...
             sprintf('%.3f%%', medVals(k)), ...
             'HorizontalAlignment', 'center', 'FontSize', 9);
    end

    title(sprintf('Settling-Time Change vs. %s%s', angleName, kwInfo));
    xlabel(angleName, 'Interpreter', 'tex');
    ylabel('Median Change [%]', 'Interpreter', 'tex');

    yTop = max(medVals + spread3sigma);
    yBot = min(medVals - spread3sigma);
    ylim([min(0, 1.3*yBot), 1.3*yTop]);

    hold off;
end



