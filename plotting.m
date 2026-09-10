function plots = plotting(graphsDir)
%PLOTTING Returns struct of all plotting functions
%
%   Usage in Project.m:
%       plots = plotting();             
%       plots.angularRate(results, To, leg_text);
%       plots.actuatorDipoles(results, To, mlim, l, SaveName="pt2_actu.png");
%       [~, improvementPct] = plots.attitudeSweepByTrial(convTime, trials);
%       plots.attitudeSweepByAngles(axisColat, axisLon, rotAngle, ...
%           axisColatGrid, axisLonGrid, rotAngleGrid, improvementPct, ...
%           kwOptimal, kwAnalytic, SaveNamePrefix="pt2_attitudeSweep");

    arguments
        graphsDir (1,1) string = "graphs"
    end
    
    if ~isfolder(graphsDir)
        mkdir(graphsDir);
    end
    
    graphsDirRegistry(graphsDir);
    
    plots.graphsDir = graphsDir;
    
    plots.angularRate           = @angularRate;
    plots.quaternionNorm        = @quaternionNorm;
    plots.actuatorDipoles       = @actuatorDipoles;
    plots.controlAuthority      = @controlAuthority;
    plots.settlingTimeVsKw      = @settlingTimeVsKw;
    plots.nearAComparison       = @nearAComparison;
    plots.attitudeSweepByTrial  = @attitudeSweepByTrial;
    plots.attitudeSweepByAngles = @attitudeSweepByAngles;
    plots.angleMarginalBar      = @angleMarginalBar;
    plots.controlAngleSweep     = @controlAngleSweep;

end

%% ------------------------------------------------------------------ %%

function fig = angularRate(results, To, leg_text, trialIdx, options)
%ANGULARRATE ||omega|| vs t/T.

    arguments
        results   cell
        To        (1,1) double {mustBePositive}
        leg_text  struct
        trialIdx  (1,:) double {mustBePositive, mustBeInteger} = 1:numel(results)
        options.SaveName (1,1) string = ""
    end
    
    fig = figure('Position', [100, 100, 600, 400]);
    hold on; grid on;
    for j = trialIdx
        plot(results{j}.Time/To, vecnorm(results{j}.omega, 2, 2), 'LineWidth', 1.5);
    end
    hold off;
    title('Angular rate over time'); xlabel('t/T [-]');
    ylabel('||\omega|| [rad/s]');
    if isfield(leg_text, 'omega') && ~isempty({leg_text.omega})
        legend(leg_text.omega, 'Interpreter', 'tex', 'Location', 'best');
    end
    
    saveFigure(fig, options.SaveName);
end

%% ------------------------------------------------------------------ %%

function fig = quaternionNorm(results, To, leg_text, trialIdx, options)
%QUATERNIONNORM ||q|| vs t/T, confirm integrator accuracy.

    arguments
        results   cell
        To        (1,1) double {mustBePositive}
        leg_text  struct
        trialIdx  (1,:) double {mustBePositive, mustBeInteger} = 1:numel(results)
        options.SaveName (1,1) string = ""
    end
    
    fig = figure('Position', [100, 100, 600, 400]);
    hold on; grid on;
    for j = trialIdx
        plot(results{j}.Time/To, vecnorm(results{j}.Q, 2, 2), 'LineWidth', 1.5);
    end
    hold off;
    title('Quaternion over time'); xlabel('t/T [-]');
    ylabel('||q||');
    if isfield(leg_text, 'Q') && ~isempty({leg_text.Q})
        legend(leg_text.Q, 'Interpreter', 'tex', 'Location', 'best');
    end
    
    saveFigure(fig, options.SaveName);
end

%% ------------------------------------------------------------------ %%

function fig = actuatorDipoles(results, To, mlim, trialIdx, options)
%ACTUATORDIPOLES figure of dipole moments with saturation limits, for a single trial.

    arguments
        results  cell
        To       (1,1) double {mustBePositive}
        mlim     (1,1) double {mustBePositive}
        trialIdx (1,1) double {mustBePositive, mustBeInteger}
        options.SaveName (1,1) string = ""
    end
    
    fig = figure('Position', [100, 100, 600, 400]);
    tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    axLabels = {'m_1', 'm_2', 'm_3'};
    
    for k = 1:3
        nexttile;
        hold on; grid on;
        plot(results{trialIdx}.Time/To, results{trialIdx}.m(:, k), 'b', 'LineWidth', 1.5);
        yline( mlim, 'r--', 'LineWidth', 1.2);
        yline(-mlim, 'r--', 'LineWidth', 1.2);
        title(axLabels{k});
        ylabel([axLabels{k} ' [A m^2]']);
        ylim([-1.15*mlim, 1.15*mlim]);
        hold off;
    end
    xlabel('t/T [-]');
    
    title(tl, 'Magnetic Actuator Dipole Moments');
    lgd = legend({'m_i', 'Saturation limit'}, 'Orientation', 'horizontal', 'Interpreter', 'tex');
    lgd.Layout.Tile = 'north';
    
    saveFigure(fig, options.SaveName);
end

%% ------------------------------------------------------------------ %%

function fig = controlAuthority(results, To, convIdx, leg_text, trialIdx, options)
%CONTROLAUTHORITY Angle alpha between omega and the geomagnetic field b,
%   plotted up to each trial's convergence index.

    arguments
        results  cell
        To       (1,1) double {mustBePositive}
        convIdx  (1,:) double {mustBePositive, mustBeInteger}
        leg_text struct
        trialIdx (1,:) double {mustBePositive, mustBeInteger} = 1:numel(results)
        options.SaveName (1,1) string = ""
    end
    
    fig = figure('Position', [100, 100, 600, 400]);
    hold on; grid on;
    
    for j = trialIdx
        omega_j = results{j}.omega;
        b_j = results{j}.b;
    
        cosAlpha = dot(omega_j, b_j, 2) ./ (vecnorm(omega_j, 2, 2) .* vecnorm(b_j, 2, 2));
        cosAlpha = min(max(cosAlpha, -1), 1); % prevent numerical errors leading to complex result
        alpha = acosd(cosAlpha);
    
        plot(results{j}.Time(1:convIdx(j))/To, alpha(1:convIdx(j)), 'LineWidth', 1.5);
    end
    
    yline(90, 'k--', 'LineWidth', 1);
    hold off;
    title('Angle Between \omega and b');
    xlabel('t/T [-]');
    ylabel('\alpha(\omega, b) [deg]');
    ylim([0, 180]);
    if isfield(leg_text, 'alpha') && ~isempty({leg_text.alpha})
        legend(leg_text.alpha, 'Location', 'best', 'Interpreter', 'tex');
    end
    
    saveFigure(fig, options.SaveName);
end

%% ------------------------------------------------------------------ %%

function fig = settlingTimeVsKw(kw, convTime, options)
%SETTLINGTIMEVSKW Settling time (t/T) vs kw for the 'optKw' sweep.
%   kw(1) should be the analytical reference gain, kw(2:end)
%   should be the numerically swept gains.

    arguments
        kw       (1,:) double
        convTime (1,:) double
        options.SaveName (1,1) string = ""
    end
    
    [Tmin, idx] = min(convTime);
    
    fig = figure('Position', [100, 100, 600, 400]);
    hold on; grid on;
    plot(kw(2:end), convTime(2:end), 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(kw(idx), Tmin, 'g*', 'DisplayName', 'Optimal');
    plot(kw(1), convTime(1), 'k*', 'DisplayName', 'Analytical');
    hold off;
    title('Settling time vs k_\omega'); ylabel('t/T');
    xlabel('k_\omega'); legend('show', 'Interpreter', 'tex', 'Location', 'best');
    
    saveFigure(fig, options.SaveName);
end

%% ------------------------------------------------------------------ %%

function [fig, pctImprovement] = nearAComparison(kw, convTime, trials, options)
%NEARACOMPARISON Bar chart comparing settling time for the analytical
%   vs. optimal kw for the 'NearA' scenarios (control +
%   perturbed Q0/omega0/theta0 cases).

    arguments
        kw       (1,:) double
        convTime (1,:) double
        trials   (1,1) double {mustBePositive, mustBeInteger}
        options.ScenarioLabels cell = {'Control', 'Initial Q', 'Initial \omega', 'Initial \theta'}    % I don't know why \omega requires $ but \theta doesn't
        options.SaveName (1,1) string = ""
    end
    
    nScenarios = trials/2;
    
    analyticalTime = convTime(1:nScenarios);
    optimalTime    = convTime(nScenarios+1:end);
    
    analyticalKw = kw(1);
    optimalKw    = kw(nScenarios+1);
    
    pctImprovement = 100*(analyticalTime - optimalTime)./analyticalTime;
    scenarioLabels = options.ScenarioLabels;
    
    fig = figure('Position', [100, 100, 600, 400]);
    hold on; grid on;
    
    barColors = repmat([0.2 0.7 0.3], nScenarios, 1);
    barColors(pctImprovement < 0, :) = repmat([0.8 0.2 0.2], sum(pctImprovement < 0), 1);
    
    b = bar(categorical(scenarioLabels, scenarioLabels), pctImprovement, 'FaceColor', 'flat');
    b.CData = barColors;
    yline(0, 'k-', 'LineWidth', 1);
    
    for k = 1:nScenarios
        offset = 0.05*max(abs(pctImprovement));
        if pctImprovement(k) < 0, offset = -offset; end
        text(k, pctImprovement(k) + offset, sprintf('%.3f%%', pctImprovement(k)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    
    title(sprintf('Settling-Time Improvement: Optimal vs. Analytical k_\\omega (k_\\omega=%.3e vs %.3e)', ...
        optimalKw, analyticalKw));
    ylabel('Improvement [%]');
    ylim([min(0, 1.4*min(pctImprovement)), 1.4*max(pctImprovement)]);
    
    saveFigure(fig, options.SaveName);
end

%% ------------------------------------------------------------------ %%

function [fig, improvementPct] = attitudeSweepByTrial(convTime, trials, options)
%ATTITUDESWEEPBYTRIAL Stem plot of settling time change (assumes first half 
%   uses analytic kw) for every sampled initial attitude, indexed by
%   trial number. Also returns improvementPct so it can be reused.

    arguments
        convTime (1,:) double
        trials   (1,1) double {mustBePositive, mustBeInteger}
        options.SaveName (1,1) string = ""
    end
    
    tAnalytic = convTime(1:trials/2);
    tOptimal  = convTime(trials/2+1:end);
    
    improvementPct = 100*(tAnalytic - tOptimal)./tAnalytic;
    
    fig = figure('Position', [100, 100, 600, 400]);
    hold on; grid on;
    
    stem(1:trials/2, improvementPct, 'filled', 'LineWidth', 1.1, 'MarkerSize', 4, 'Color', 'g');
    
    yline(0, 'k-', 'LineWidth', 1.2);
    xlabel('Initial-attitude sample');
    ylabel('Settling-time change [%]');
    title('Numerical k_\omega Improvement over Analytical k_\omega');
    
    saveFigure(fig, options.SaveName);
end

%% ------------------------------------------------------------------ %%

function figs = attitudeSweepByAngles(axisColat, axisLon, rotAngle, ...
    axisColatGrid, axisLonGrid, rotAngleGrid, improvementPct, ...
    kwOptimal, kwAnalytic, options)
%ATTITUDESWEEPBYANGLES Median +/- 3*std change bar charts for each swept angle.

    arguments
        axisColat       (1,:) double
        axisLon         (1,:) double
        rotAngle        (1,:) double
        axisColatGrid   (1,:) double
        axisLonGrid     (1,:) double
        rotAngleGrid    (1,:) double
        improvementPct  (1,:) double
        kwOptimal       (1,1) double {mustBePositive}
        kwAnalytic      (1,1) double {mustBePositive}
        options.SaveNamePrefix (1,1) string = ""
    end
    
    colatMed = zeros(size(axisColat)); colatStd = zeros(size(axisColat));
    lonMed   = zeros(size(axisLon));   lonStd   = zeros(size(axisLon));
    rotMed   = zeros(size(rotAngle));  rotStd   = zeros(size(rotAngle));
    
    for j = 1:length(axisColat)
        colatMed(j) = median(improvementPct(axisColatGrid == axisColat(j)));
        colatStd(j) = std(improvementPct(axisColatGrid == axisColat(j)));
    end
    for j = 1:length(axisLon)
        lonMed(j) = median(improvementPct(axisLonGrid == axisLon(j)));
        lonStd(j) = std(improvementPct(axisLonGrid == axisLon(j)));
    end
    for j = 1:length(rotAngle)
        rotMed(j) = median(improvementPct(rotAngleGrid == rotAngle(j)));
        rotStd(j) = std(improvementPct(rotAngleGrid == rotAngle(j)));
    end
    
    kwInfo = sprintf(' (k_\\omega=%.3e vs %.3e)', kwOptimal, kwAnalytic);
    prefix = options.SaveNamePrefix;
    
    figs(1) = angleMarginalBar(axisColat, colatMed, colatStd, 'deg', ...
        'Axis Colatitude', kwInfo, SaveName=nameOrEmpty(prefix, "_colat.png"));
    figs(2) = angleMarginalBar(axisLon, lonMed, lonStd, 'deg', ...
        'Axis Longitude', kwInfo, SaveName=nameOrEmpty(prefix, "_lon.png"));
    figs(3) = angleMarginalBar(rotAngle, rotMed, rotStd, 'deg', ...
        'Rotation Angle', kwInfo, SaveName=nameOrEmpty(prefix, "_rotAngle.png"));

end

%% ------------------------------------------------------------------ %%

function fig = angleMarginalBar(angleValues, medVals, stdVals, angleUnit, angleName, kwInfo, options)
%ANGLEMARGINALBAR Median change bar chart with error bars for a swept angle.

    arguments
        angleValues (1,:) double    % the grid values for this angle
        medVals     (1,:) double    % median change for each angle
        stdVals     (1,:) double {mustBeNonnegative}    % standard deviation for each angle     
        angleUnit   (1,1) string    
        angleName   (1,1) string
        kwInfo      (1,1) string
        options.SaveName (1,1) string = ""
    end
    
    nBars = numel(angleValues);
    spread3sigma = 3*stdVals;
    
    labels = compose('%g %s', angleValues(:), angleUnit);
    
    barColors = repmat([0.2 0.7 0.3], nBars, 1);
    barColors(medVals < 0, :) = repmat([0.8 0.2 0.2], sum(medVals < 0), 1);
    
    fig = figure('Position', [100, 100, 600, 400]);
    hold on; grid on;
    
    b = bar(categorical(labels, labels), medVals, 'FaceColor', 'flat');
    b.CData = barColors;
    
    errorbar(1:nBars, medVals, spread3sigma, 'k', 'LineStyle', 'none', ...
        'LineWidth', 1.2, 'CapSize', 8, 'HandleVisibility', 'off');
    
    yline(0, 'k-', 'LineWidth', 1);
    
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
    
    saveFigure(fig, options.SaveName);
end

%% ------------------------------------------------------------------ %%

function [fig, change, alpha0Sorted] = controlAngleSweep(alpha0, convTime, trials, options)
%CONTROLANGLESWEEP Settling-time change vs. initial angle between omega
%   and b.

    arguments
        alpha0 (1,:) double
        convTime (1,:) double
        trials (1,1) double {mustBePositive, mustBeInteger}
        options.SaveName (1,1) string = ""
    end
    
    tAnalytic = convTime(1:trials/2);
    tOptimal  = convTime(trials/2+1:end);
    
    [alpha0Sorted, sortIdx] = sort(alpha0);
    tAnalytic = tAnalytic(sortIdx);
    tOptimal  = tOptimal(sortIdx);
    
    % improvementPct = (tAnalytic - tOptimal) ./ tAnalytic * 100;
    change = tAnalytic - tOptimal;
    disp(mean(change))
    
    fig = figure('Position', [100, 100, 600, 400]);
    hold on; grid on;
    
    plot(alpha0Sorted, change, 'b', 'LineWidth', 1.5);
    yline(0, 'k-', 'LineWidth', 1);
    xline(90, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    hold off;
    title('Settling-Time Change vs. \alpha_0');
    xlabel('\alpha_0 [deg]', 'Interpreter', 'tex');
    ylabel('Change in orbits [-]');
    xlim([0, 180]);
    
    saveFigure(fig, options.SaveName);
end

%% ------------------------------------------------------------------ %%

function n = nameOrEmpty(prefix, suffix)
%NAMEOREMPTY Small helper so attitudeSweepByAngles can skip saving
%   (empty SaveName) when no prefix was requested.

    arguments
        prefix (1,1) string
        suffix (1,1) string
    end
    
    if prefix == ""
        n = "";
    else
        n = prefix + suffix;
    end
end

%% ------------------------------------------------------------------ %%

function saveFigure(fig, saveName)
%SAVEFIGURE Shared export helper.

    arguments
        fig      matlab.ui.Figure
        saveName (1,1) string = ""
    end
    
    if saveName ~= ""
        exportgraphics(fig, fullfile(graphsDirRegistry(), saveName), 'Resolution', 300);
    end
end

%% ------------------------------------------------------------------ %%

function out = graphsDirRegistry(newDir)
%GRAPHSDIRREGISTRY Get/set the shared graphs directory.

    persistent storedDir
    if nargin == 1
        storedDir = newDir;
    end
    if isempty(storedDir)
        storedDir = "graphs";
    end
    out = storedDir;
end
