clearvars; close all; clc;

%% PAPER SIMULATIONS
disp("============== Simulations From The Paper ==============")

%% data from table 1

% Which simulation to run: 
% paper - replicate paper sim 
% ecc - test different eccentricities
% optKw - find optimal kw for scenario A
% NearA - optimal value vs analytic at scenarios like A
sim_type = "NearA"; 
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
N       = 6;        % [-] number of orbits
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
            "theta", ['\theta_A';'\theta_B';'\theta_C']);

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
            "theta", compose('e=%0.2f',ecc));
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
            "theta", compose('k_\\omega=%e',kws));

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
            "theta", compose('k_\\omega=%e',kw));
end



%% Simulation

results = cell(1,trials);
convIdx = zeros(1, trials);
convTime = zeros(1, trials);
opts = odeset(Refine=1, Stats='on');
tspan = [0, N*To];

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


%% angular rate
handle = figure('Position', [100, 100, 600, 400]);
hold on; grid on;
for j=1:length(results)
    plot(results{j}.Time/To, vecnorm(results{j}.omega,2,2), LineWidth=1.5);
end
hold off; title("Angular rate over time"); xlabel('t/T [-]');
ylabel('||\omega|| [rad/s]'); legend(leg_text.omega, Interpreter="tex", Location="best")
% exportgraphics(handle, "graphs/pt2_omega.png", "Resolution", 300);

%% Quaternion size (verify sim accuracy)
handle = figure('Position', [100, 100, 600, 400]);
hold on; grid on;
for j=1:length(results)
    plot(results{j}.Time/To, vecnorm(results{j}.Q,2,2), LineWidth=1.5);
end
hold off; title("Quaternion over time"); xlabel('t/T [-]');
ylabel('||q||'); legend(leg_text.Q, Interpreter="tex", Location="best")
% exportgraphics(handle, "graphs/pt2_Q.png", "Resolution", 300);

%% Actuators
l = 1;  % Which sim to plot
handle = figure('Position', [100, 100, 600, 400]);
tl = tiledlayout(3, 1, TileSpacing='compact', Padding='compact');
axLabels = {'m_1', 'm_2', 'm_3'};

for k = 1:3
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

for j = 1:length(results)
    omega_j = results{j}.omega;  
    b_j     = results{j}.b;       

    cosTheta = dot(omega_j, b_j, 2) ./ (vecnorm(omega_j,2,2).*vecnorm(b_j,2,2));
    cosTheta = min(max(cosTheta, -1), 1);   % prevent numerical errors leading to complex result
    theta = acosd(cosTheta);               

    plot(results{j}.Time(1:convIdx(j))/To, theta(1:convIdx(j)), 'LineWidth', 1.5);
end

yline(90, 'k--', 'LineWidth', 1);
hold off;
title('Angle Between \omega and b');
xlabel('t/T [-]');
ylabel('\theta(\omega, b) [deg]');
ylim([0, 180]);
legend(leg_text.theta, Location='best', Interpreter='tex');

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

if strcmp(sim_type, "NearA")
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
end

