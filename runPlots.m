function runPlots(plots, enabled, ctx, saveName)
%RUNPLOTS Calls all the plotting function asked for.

    arguments
        plots (1,1) struct
        enabled (1,1) struct
        ctx (1,1) struct
        saveName (1,1) struct = struct()
    end

    % Function mapping struct
    dispatch = struct();
    dispatch.angularRate           = @(p,c,s) p.angularRate(c.results, c.To, c.leg_text, SaveName=s);
    dispatch.controlTorque         = @(p,c,s) p.controlTorque(c.results, c.To, c.leg_text, SaveName=s);
    dispatch.quaternionNorm        = @(p,c,s) p.quaternionNorm(c.results, c.To, c.leg_text, SaveName=s);
    dispatch.actuatorDipoles       = @(p,c,s) p.actuatorDipoles(c.results, c.To, c.mlim, c.l, SaveName=s);
    dispatch.controlAuthority      = @(p,c,s) p.controlAuthority(c.results, c.To, c.convIdx, c.leg_text, SaveName=s);
    dispatch.settlingTimeVsKw      = @(p,c,s) p.settlingTimeVsKw(c.kw, c.convTime, SaveName=s);
    dispatch.nearAComparison       = @(p,c,s) p.nearAComparison(c.kw, c.convTime, c.trials, SaveName=s);
    dispatch.attitudeSweepByTrial  = @(p,c,s) p.attitudeSweepByTrial(c.convTime, c.trials, SaveName=s);
    dispatch.attitudeSweepByAngles = @(p,c,s) p.attitudeSweepByAngles( ...
        c.axisColat, c.axisLon, c.rotAngle, c.axisColatGrid, c.axisLonGrid, c.rotAngleGrid, ...
        c.improvementPct, c.kwOptimal, c.kwAnalytic, SaveName=s);
    dispatch.controlAngleSweep     = @(p,c,s) p.controlAngleSweep(c.alpha0, c.convTime, c.trials, SaveName=s);

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
            if isfield(saveName, name)
                dispatch.(name)(plots, ctx, saveName.(name));
            else
                dispatch.(name)(plots, ctx, "");
            end
        catch ME
            warning("runPlots:plotFailed", ...
                "Plot '%s' failed (%s) - skipping...", name, ME.message);
        end
    end

end
