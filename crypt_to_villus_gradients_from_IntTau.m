function crypt_to_villus_gradients_from_IntTau
%CRYPT_TO_VILLUS_GRADIENTS_FROM_INTTAU Analyze crypt->villus gradients for Int/TauPhase/TauModulation across groups.
%   Each "group" corresponds to a folder containing CSV files (one row per cell)
%   with required columns: Int, TauPhase, TauModulation, Xcoord, Ycoord.
%
%   Workflow overview (copy/paste friendly):
%     1) When prompted, select one or more folders. Each folder = one group.
%        If you pick a parent folder without CSVs, each CSV-containing subfolder
%        under it becomes a group automatically.
%     2) For every group, a scatter plot of all cells opens. Click crypt base
%        first, villus tip second to define that group's axis (same logic as the
%        original single-group script).
%     3) Cells are projected to the axis, normalized to [0,1], binned, and mean
%        ± SEM are computed per bin for Int, TauPhase, and TauModulation.
%     4) A combined Prism-like figure overlays mean ± SEM curves for all groups
%        per metric, with group names taken from folder names.
%     5) Grouped bin statistics are written to ./output_IntTau/grouped_bin_summary.csv.
%
%   Output directory: ./output_IntTau/
%
%   This function relies only on base MATLAB (R2021b+).

%% --------------------- USER SETTINGS ---------------------
numBins    = 50;         % number of bins between crypt (0) and villus (1)
outFolder  = './output_IntTau';
smoothSpan = 0.12;       % fraction of bins used for smoothing plotted curves

if ~exist(outFolder,'dir'), mkdir(outFolder); end

%% --------------------- GROUP (FOLDER) SELECTION ---------------------
% Allow the user to select one or more folders. Each folder = one group.
% If a selected folder contains no CSV files, each of its immediate
% subfolders is treated as a separate group (facilitating "choose a
% parent folder to analyze all subfolders" workflows).
groupDirsInput = selectGroupFolders();
groupDirs = resolveGroupFolders(groupDirsInput);
if isempty(groupDirs)
    disp('No group folders selected. Exiting.');
    return;
end

groupStruct = struct('Name',{},'Table',{},'BinIdx',{},'Stats',{});
requiredVars = {'Int','TauPhase','TauModulation','Xcoord','Ycoord'};

%% --------------------- LOAD DATA & DEFINE AXIS PER GROUP ---------------------
for g = 1:numel(groupDirs)
    grpDir = groupDirs{g};
    [~, grpName] = fileparts(grpDir);

    csvList = dir(fullfile(grpDir, '*.csv'));
    if isempty(csvList)
        warning('No CSV files found in folder %s. Skipping.', grpDir);
        continue;
    end

    grpTable = table();
    for k = 1:numel(csvList)
        fpath = fullfile(grpDir, csvList(k).name);
        T = readtable(fpath);
        assert(all(ismember(requiredVars, T.Properties.VariableNames)), ...
            'File %s must contain columns: %s', fpath, strjoin(requiredVars, ', '));
        T.SourceFile = repmat(string(csvList(k).name), height(T), 1);
        T.Group = repmat(string(grpName), height(T), 1);
        grpTable = [grpTable; T]; %#ok<AGROW>
    end

    if isempty(grpTable)
        warning('No data loaded for group %s. Skipping.', grpName);
        continue;
    end

    % Axis definition per group (same logic as original single-group code).
    [baseXY, tipXY] = getAxisFromLine(grpTable.Xcoord, grpTable.Ycoord, grpName);
    v = tipXY - baseXY;
    vlen2 = sum(v.^2);
    assert(vlen2 > 0, 'Crypt and villus clicks are identical for group %s.', grpName);

    % Project and clip to [0,1].
    pts = [grpTable.Xcoord, grpTable.Ycoord];
    proj = ((pts - baseXY) * v.') / vlen2;
    grpTable.Pos01 = min(max(proj,0),1);

    groupStruct(end+1).Name = grpName; %#ok<AGROW>
    groupStruct(end).Table = grpTable;
end

if isempty(groupStruct)
    disp('No valid group data found. Exiting.');
    return;
end

%% --------------------- BINNING & GROUP-WISE STATISTICS ---------------------
edges = linspace(0,1,numBins+1);
binCenters = edges(1:end-1) + diff(edges)/2;
metrics = {'Int','TauPhase','TauModulation'};

summaryRows = table();
for g = 1:numel(groupStruct)
    grpTable = groupStruct(g).Table;
    [~,~,binIdx] = histcounts(grpTable.Pos01, edges);
    groupStruct(g).BinIdx = binIdx;

    stats = struct();
    for i = 1:numel(metrics)
        stats.(metrics{i}) = computeBinStats(grpTable.(metrics{i}), binIdx, numBins);
    end
    groupStruct(g).Stats = stats;

    % Assemble per-bin summary for CSV output.
    grpCol = repmat(string(groupStruct(g).Name), numBins, 1);
    summaryRows = [summaryRows; table(grpCol, (1:numBins).', binCenters(:), ...
        stats.Int.N(:), stats.Int.Mean(:), stats.Int.SEM(:), ...
        stats.TauPhase.N(:), stats.TauPhase.Mean(:), stats.TauPhase.SEM(:), ...
        stats.TauModulation.N(:), stats.TauModulation.Mean(:), stats.TauModulation.SEM(:), ...
        'VariableNames', {'Group','Bin','BinCenter', ...
        'N_Int','Mean_Int','SEM_Int', ...
        'N_TauPhase','Mean_TauPhase','SEM_TauPhase', ...
        'N_TauModulation','Mean_TauModulation','SEM_TauModulation'})]; %#ok<AGROW>
end

writetable(summaryRows, fullfile(outFolder, 'grouped_bin_summary.csv'));

%% --------------------- PLOT MULTI-GROUP MEAN ± SEM ---------------------
fFig = figure('Color','w','Position',[200 200 1200 500]);
set(fFig,'DefaultAxesFontName','Arial', ...
         'DefaultAxesFontSize',12, ...
         'DefaultAxesFontWeight','bold', ...
         'DefaultAxesLineWidth',2, ...
         'DefaultAxesTickDir','out', ...
         'DefaultAxesBox','off');

colors = lines(numel(groupStruct));
ylabels = {'Intensity (a.u.)','TauPhase (ns)','TauModulation (ns)'};

for i = 1:numel(metrics)
    ax = subplot(1,3,i,'Parent',fFig); hold(ax,'on'); box(ax,'off');
    for g = 1:numel(groupStruct)
        mRaw = groupStruct(g).Stats.(metrics{i}).Mean;
        sRaw = groupStruct(g).Stats.(metrics{i}).SEM;
        % Smooth only for visualization; CSV retains raw bin statistics.
        m = smoothForPlot(mRaw, smoothSpan);
        s = smoothForPlot(sRaw, smoothSpan);
        c = colors(g,:);
        fill(ax, [binCenters fliplr(binCenters)], [m-s; flipud(m+s)].', ...
            c, 'FaceAlpha',0.15, 'EdgeColor','none');
        plot(ax, binCenters, m, 'Color', c, 'LineWidth', 2, 'DisplayName', groupStruct(g).Name);
    end
    xlabel(ax, 'Normalized position (0 = crypt, 1 = villus)');
    ylabel(ax, ylabels{i});
    title(ax, sprintf('%s mean \pm SEM', metrics{i}), 'Interpreter','none');
    xlim(ax,[0 1]);
    legend(ax,'Location','best','Box','off');
end

try
    exportgraphics(fFig, fullfile(outFolder,'grouped_mean_sem.png'), 'Resolution', 600);
catch
    print(fFig, fullfile(outFolder,'grouped_mean_sem.png'), '-dpng','-r600');
end

end

%% --------------------- HELPERS ---------------------
function folders = selectGroupFolders()
%SELECTGROUPFOLDERS Let the user choose one or more folders as groups.
%   Uses repeated UIGETDIR calls so users can pick any set of folders. The
%   loop ends when the user cancels the dialog.
    folders = {};
    while true
        pth = uigetdir(pwd, 'Select a folder containing CSV files for one group');
        if isequal(pth,0)
            break;
        end
        folders{end+1} = pth; %#ok<AGROW>
        choice = questdlg('Add another group folder?', 'Select groups', ...
            'Yes','No','No');
        if isempty(choice) || strcmpi(choice,'No')
            break;
        end
    end
    folders = unique(folders, 'stable');
end

function folders = resolveGroupFolders(inputDirs)
%RESOLVEGROUPFOLDERS Expand selected folders into actual group folders.
%   If a selected folder contains CSV files, it is treated as one group.
%   Otherwise, each immediate subfolder that contains CSV files becomes a
%   group. This mirrors "select parent folder to analyze all subfolders".
    folders = {};
    if nargin == 0 || isempty(inputDirs)
        return;
    end
    for i = 1:numel(inputDirs)
        p = inputDirs{i};
        if ~isfolder(p)
            warning('Path %s is not a folder. Skipping.', p);
            continue;
        end

        csvHere = dir(fullfile(p, '*.csv'));
        if ~isempty(csvHere)
            folders{end+1} = p; %#ok<AGROW>
            continue;
        end

        % No CSVs in the selected folder; look one level down.
        subDirs = dir(p);
        subDirs = subDirs([subDirs.isdir]);
        subDirs = subDirs(~ismember({subDirs.name}, {'.','..'}));
        added = false;
        for s = 1:numel(subDirs)
            subPath = fullfile(p, subDirs(s).name);
            if ~isempty(dir(fullfile(subPath, '*.csv')))
                folders{end+1} = subPath; %#ok<AGROW>
                added = true;
            end
        end
        if ~added
            warning('No CSV files found in %s or its immediate subfolders.', p);
        end
    end
    folders = unique(folders, 'stable');
end

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

function [baseXY, tipXY] = getAxisFromLine(x, y, groupName)
%GETAXISFROMLINE Ask the user to draw the crypt-to-villus axis on a scatter plot.
    fScatter = figure('Color','w');
    scatter(x, y, 8, 'k', 'filled');
    axis equal;
    xlabel('Xcoord'); ylabel('Ycoord');
    title(sprintf('%s: draw axis (crypt base to villus tip)', groupName), ...
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
            error('Axis selection aborted for group %s.', groupName);
        end
        if numel(xClick) < 2
            close(fScatter);
            error('Two points are required to define the axis for group %s.', groupName);
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
