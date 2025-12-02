function crypt_to_villus_gradients_from_IntTau
%CRYPT_TO_VILLUS_GRADIENTS_FROM_INTTAU Analyze crypt->villus gradients for Int/TauPhase/TauModulation.
%   Interactively selects CSV files (one row per cell) with required columns:
%   Int, TauPhase, TauModulation, Xcoord, Ycoord.
%   The user defines the crypt-to-villus axis by drawing a line on the
%   scatter plot (crypt end first, villus end second). The function projects all cells
%   onto this axis, bins the normalized positions, computes pooled mean and
%   SEM per bin for each metric, and saves a single Prism-like figure with
%   the three pooled mean ± SEM curves plus a CSV summary of bin statistics.
%
%   Output directory: ./output_IntTau/
%
%   This function relies only on base MATLAB (R2021b+).

%% --------------------- USER SETTINGS ---------------------
numBins    = 50;         % number of bins between crypt (0) and villus (1)
outFolder  = './output_IntTau';
smoothSpan = 0.12;       % fraction of bins used for smoothing plotted curves

if ~exist(outFolder,'dir'), mkdir(outFolder); end

%% --------------------- FILE SELECTION ---------------------
[files, pth] = uigetfile('*.csv','Select CSV files','MultiSelect','on');
if isequal(files,0)
    % user cancelled
    return;
end
if ischar(files)
    files = {files};
end
fileList = fullfile(pth, files);

%% --------------------- LOAD, VALIDATE, AND GET AXIS PER FILE ---------------------
requiredVars = {'Int','TauPhase','TauModulation','Xcoord','Ycoord'};
allRows = table();
for k = 1:numel(fileList)
    fpath = fileList{k};
    T = readtable(fpath);
    assert(all(ismember(requiredVars, T.Properties.VariableNames)), ...
        'File %s must contain columns: %s', fpath, strjoin(requiredVars, ', '));
    [~,fname,ext] = fileparts(fpath);
    sampleName = string([fname ext]);

    % Ask the user to draw the crypt-to-villus axis for this sample.
    [baseXY, tipXY] = getAxisFromLine(T.Xcoord, T.Ycoord, sampleName);
    v = tipXY - baseXY;
    vlen2 = sum(v.^2);
    assert(vlen2 > 0, 'Crypt and villus clicks are identical for sample %s.', sampleName);

    pts = [T.Xcoord, T.Ycoord];
    proj = ((pts - baseXY) * v.') / vlen2; % projection scalar
    T.Pos01 = min(max(proj,0),1);          % clip to [0,1]

    T.Sample = repmat(sampleName, height(T), 1);
    allRows = [allRows; T]; %#ok<AGROW>
end

if isempty(allRows)
    disp('No data loaded. Exiting.');
    return;
end

%% --------------------- BINNING ---------------------
edges = linspace(0,1,numBins+1);
binCenters = edges(1:end-1) + diff(edges)/2;
[~,~,binIdx] = histcounts(allRows.Pos01, edges);

metrics = {'Int','TauPhase','TauModulation'};
stats = struct();
for i = 1:numel(metrics)
    stats.(metrics{i}) = computeBinStats(allRows.(metrics{i}), binIdx, numBins);
end

%% --------------------- CSV SUMMARY ---------------------
summaryTbl = table((1:numBins).', binCenters(:), ...
    stats.Int.N(:), stats.Int.Mean(:), stats.Int.SEM(:), ...
    stats.TauPhase.N(:), stats.TauPhase.Mean(:), stats.TauPhase.SEM(:), ...
    stats.TauModulation.N(:), stats.TauModulation.Mean(:), stats.TauModulation.SEM(:), ...
    'VariableNames', {'Bin','BinCenter', ...
    'N_Int','Mean_Int','SEM_Int', ...
    'N_TauPhase','Mean_TauPhase','SEM_TauPhase', ...
    'N_TauModulation','Mean_TauModulation','SEM_TauModulation'});
writetable(summaryTbl, fullfile(outFolder, 'pooled_bin_summary.csv'));

%% --------------------- PLOT MEAN ± SEM ---------------------
fFig = figure('Color','w','Position',[200 200 1200 500]);
set(fFig,'DefaultAxesFontName','Arial', ...
         'DefaultAxesFontSize',12, ...
         'DefaultAxesFontWeight','bold', ...
         'DefaultAxesLineWidth',2, ...
         'DefaultAxesTickDir','out', ...
         'DefaultAxesBox','off');

ylabs = {'Intensity (a.u.)','TauPhase (ns)','TauModulation (ns)'};
for i = 1:numel(metrics)
    ax = subplot(1,3,i,'Parent',fFig); hold(ax,'on'); box(ax,'off');
    mRaw = stats.(metrics{i}).Mean;
    sRaw = stats.(metrics{i}).SEM;
    % Smooth only for visualization; CSV retains the raw bin statistics.
    m = smoothForPlot(mRaw, smoothSpan);
    s = smoothForPlot(sRaw, smoothSpan);
    fill(ax, [binCenters fliplr(binCenters)], [m-s; flipud(m+s)].', ...
        [0 0 0], 'FaceAlpha',0.15, 'EdgeColor','none');
    plot(ax, binCenters, m, 'k-', 'LineWidth',2);
    xlabel(ax, 'Normalized position (0 = crypt, 1 = villus)');
    ylabel(ax, ylabs{i});
    title(ax, sprintf('Pooled %s mean \pm SEM', metrics{i}), 'Interpreter','none');
    xlim(ax,[0 1]);
end

try
    exportgraphics(fFig, fullfile(outFolder,'pooled_mean_sem.png'), 'Resolution', 600);
catch
    print(fFig, fullfile(outFolder,'pooled_mean_sem.png'), '-dpng','-r600');
end

end

%% --------------------- HELPERS ---------------------
function S = computeBinStats(values, binIdx, numBins)
%COMPUTEBINSTATS Compute N, mean, and SEM per bin for a vector of values.
    Mean = nan(numBins,1);
    SEM  = nan(numBins,1);
    N    = zeros(numBins,1);
    for b = 1:numBins
        mask = (binIdx == b);
        x = values(mask);
        x = x(~isnan(x));
        N(b) = numel(x);
        if N(b) > 0
            Mean(b) = mean(x);
            SEM(b)  = std(x) / sqrt(N(b));
        end
    end
    S = struct('Mean',Mean,'SEM',SEM,'N',N);
end

function [baseXY, tipXY] = getAxisFromLine(x, y, sampleName)
%GETAXISFROMLINE Ask the user to draw the crypt-to-villus axis on a scatter plot.
    fScatter = figure('Color','w');
    scatter(x, y, 8, 'k', 'filled');
    axis equal;
    xlabel('Xcoord'); ylabel('Ycoord');
    title(sprintf('%s: draw axis (crypt base to villus tip)', sampleName), ...
        'Interpreter','none');
    set(gca,'FontName','Arial','FontWeight','bold','LineWidth',2,'TickDir','out');

    roiLine = [];
    try
        roiLine = drawline('Color',[0.2 0.6 1],'LineWidth',2); %#ok<DRAWLINE>
        wait(roiLine); % Wait for double-click or right-click
        pos = roiLine.Position;
    catch
        % Fallback to simple two-click input if drawline is unavailable
        warning('drawline unavailable; falling back to two-click axis selection.');
        try
            [xClick,yClick] = ginput(2);
        catch
            close(fScatter);
            error('Axis selection aborted for sample %s.', sampleName);
        end
        if numel(xClick) < 2
            close(fScatter);
            error('Two points are required to define the axis for sample %s.', sampleName);
        end
        pos = [xClick(:), yClick(:)];
    end

    baseXY = [pos(1,1), pos(1,2)]; % crypt end
    tipXY  = [pos(end,1), pos(end,2)]; % villus end
    hold on; plot([baseXY(1) tipXY(1)], [baseXY(2) tipXY(2)], 'b-', 'LineWidth',2);
    drawnow;
    close(fScatter);
end

function y = smoothForPlot(x, span)
%SMOOTHFORPLOT Light smoothing helper for plotting jagged curves.
%   Uses LOWESS via SMOOTHDATA when available; falls back to a moving mean.
    if nargin < 2 || isempty(span)
        span = 0.12;
    end
    y = x;
    if isempty(x) || all(isnan(x))
        return;
    end
    win = max(1, round(span * numel(x)));
    if win < 2
        return;
    end
    try
        y = smoothdata(x, 'lowess', win, 'omitnan');
    catch
        y = movmean(x, win, 'omitnan');
    end
end

end
