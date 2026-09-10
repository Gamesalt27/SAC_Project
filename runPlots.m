function runPlots(plots, enabled, ctx)
%RUNPLOTS Calls all the plotting function asked for.

    arguments
        plots (1,1) struct
        enabled (1,1) struct
        ctx (1,1) struct
    end

    % Function mapping struct
    dispatch = struct();
    dispatch.angularRate          = @(p,c) p.angularRate(c.results, c.To, c.leg_text);
    dispatch.quaternionNorm       = @(p,c) p.quaternionNorm(c.results, c.To, c.leg_text);
    dispatch.actuatorDipoles      = @(p,c) p.actuatorDipoles(c.results, c.To, c.mlim, c.l);
    dispatch.controlAuthority     = @(p,c) p.controlAuthority(c.results, c.To, c.convIdx, c.leg_text);
    dispatch.settlingTimeVsKw     = @(p,c) p.settlingTimeVsKw(c.kw, c.convTime);
    dispatch.nearAComparison      = @(p,c) p.nearAComparison(c.kw, c.convTime, c.trials);
    dispatch.attitudeSweepByTrial = @(p,c) p.attitudeSweepByTrial(c.convTime, c.trials);
    dispatch.attitudeSweepByAngles = @(p,c) p.attitudeSweepByAngles( ...
        c.axisColat, c.axisLon, c.rotAngle, ...
        c.axisColatGrid, c.axisLonGrid, c.rotAngleGrid, ...
        c.improvementPct, c.kwOptimal, c.kwAnalytic);
    dispatch.controlAngleSweep    = @(p,c) p.controlAngleSweep(c.alpha0, c.convTime, c.trials);

    names = fieldnames(enabled);
    for k = 1:numel(names)
        name = names{k};
    
        if ~islogical(enabled.(name)) && ~isnumeric(enabled.(name))
            warning("runPlots:badFlag", ...
                "enabled.%s is not a logical/numeric flag -- skipping.", name);
            continue
        end
    
        if ~enabled.(name)
            continue
        end
    
        if ~isfield(dispatch, name)
            warning("runPlots:unknownPlot", ...
                "'%s' has no entry in the dispatch table - skipping...", name);
            continue
        end
    
        if ~isfield(plots, name)
            warning("runPlots:missingHandle", ...
                "plots.%s does not exist - skipping...", name);
            continue
        end
    
        try
            dispatch.(name)(plots, ctx);
        catch ME
            warning("runPlots:plotFailed", ...
                "Plot '%s' failed (%s) - skipping...", name, ME.message);
        end
    end

end
