function crypt_to_villus_gradients_from_IntTau
%CRYPT_TO_VILLUS_GRADIENTS_FROM_INTTAU Analyze crypt->villus gradients for Int/TauPhase/TauModulation.
%   Interactively selects CSV files (one row per cell) with required columns:
%   Int, TauPhase, TauModulation, Xcoord, Ycoord.
%   The user defines the crypt-to-villus axis by clicking two points on a
%   scatter plot (crypt first, villus second). The function projects all cells
%   onto this axis, bins the normalized positions, computes pooled mean and
%   SEM per bin for each metric, and saves a single Prism-like figure with
%   the three pooled mean ± SEM curves plus a CSV summary of bin statistics.
%
%   Output directory: ./output_IntTau/
%
%   This function relies only on base MATLAB (R2021b+).

%% --------------------- USER SETTINGS ---------------------
numBins   = 50;         % number of bins between crypt (0) and villus (1)
outFolder = './output_IntTau';

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

%% --------------------- LOAD & CONCATENATE ---------------------
requiredVars = {'Int','TauPhase','TauModulation','Xcoord','Ycoord'};
allRows = table();
for k = 1:numel(fileList)
    fpath = fileList{k};
    T = readtable(fpath);
    assert(all(ismember(requiredVars, T.Properties.VariableNames)), ...
        'File %s must contain columns: %s', fpath, strjoin(requiredVars, ', '));
    [~,fname,ext] = fileparts(fpath);
    T.Sample = repmat(string([fname ext]), height(T), 1);
    allRows = [allRows; T]; %#ok<AGROW>
end

if isempty(allRows)
    disp('No data loaded. Exiting.');
    return;
end

%% --------------------- AXIS DEFINITION (USER CLICKS) ---------------------
fScatter = figure('Color','w');
scatter(allRows.Xcoord, allRows.Ycoord, 8, 'k', 'filled');
axis equal;
xlabel('Xcoord'); ylabel('Ycoord');
title('Click crypt base (first) then villus tip (second)');
set(gca,'FontName','Arial','FontWeight','bold','LineWidth',2,'TickDir','out');
try
    [xClick,yClick] = ginput(2);
catch
    close(fScatter);
    error('Axis selection aborted.');
end
if numel(xClick) < 2
    close(fScatter);
    error('Two points are required to define the axis.');
end
baseXY = [xClick(1), yClick(1)]; % crypt
tipXY  = [xClick(2), yClick(2)]; % villus
close(fScatter);

v = tipXY - baseXY;
vlen2 = sum(v.^2);
assert(vlen2 > 0, 'Crypt and villus clicks are identical.');

pts = [allRows.Xcoord, allRows.Ycoord];
proj = ((pts - baseXY) * v.') / vlen2; % projection scalar
proj01 = min(max(proj,0),1);           % clip to [0,1]
allRows.Pos01 = proj01;

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
    m = stats.(metrics{i}).Mean;
    s = stats.(metrics{i}).SEM;
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

end
