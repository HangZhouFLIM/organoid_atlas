function crypt_to_villus_gradients_combined
%CRYPTO_VILLUS_GRADIENTS_COMBINED Interactive analysis for multiple metrics.
%   Combines previous standalone scripts for MeanInt, MeanTauPhase and
%   MeanTauModulation into one pipeline. Users select CSV files containing
%   per-cell measurements and define the crypt base and villus tip cells.
%   The function computes normalized positions along the crypt->villus axis
%   and generates summaries and plots for each metric.
%
%   Output structure:
%     ./output/<Metric>/ ... files and figures for each metric.
%
%   Expected input columns in each CSV:
%       Id, CentroidX, CentroidY, MeanInt, MeanTauPhase, MeanTauModulation
%
%   Example:
%       crypt_to_villus_gradients_combined
%
% Author: OpenAI Assistant

%% --------------------- USER SETTINGS ---------------------
numBins        = 50;
pctList        = [10 25 50 75 90];
outFolder      = './output';
lowessSpan     = 0.2;

if ~exist(outFolder,'dir'), mkdir(outFolder); end

%% Metrics to analyse
metrics = {
    struct('field','MeanInt',         'label','Mean intensity');
    struct('field','MeanTauPhase',    'label','MeanTauPhase');
    struct('field','MeanTauModulation','label','MeanTauModulation')
    };

%% --------------------- SELECT FILES + AXES ---------------------
[fileList, axisDef] = chooseFilesAndAxes();
if isempty(fileList), disp('No files selected. Exiting.'); return; end

%% --------------------- PROCESS SAMPLES (projection only) ---------------------
numFiles = numel(fileList);
tables = cell(numFiles,1);
sampleNames = cell(numFiles,1);

% Expand required metric fields
req = [{'Id','CentroidX','CentroidY'}];
for m = 1:numel(metrics)
    req{end+1} = metrics{m}.field; %#ok<AGROW>
end

for k = 1:numFiles
    fpath = fileList{k};
    [~,fname,ext] = fileparts(fpath);
    sampleName = [fname ext];
    T = readtable(fpath);

    assert(all(ismember(req, T.Properties.VariableNames)), ...
        'File %s must contain columns: %s', fpath, strjoin(req, ', '));

    cryptId  = axisDef{k}.cryptId;
    villusId = axisDef{k}.villusId;
    [baseXY, tipXY] = getAxisFromIds(T, cryptId, villusId);

    v = (tipXY - baseXY);
    vlen2 = sum(v.^2);
    assert(vlen2 > 0, 'Chosen crypt and villus points are identical for %s', sampleName);

    pts = [T.CentroidX, T.CentroidY];
    proj = ((pts - baseXY) * v.') / vlen2;
    proj01 = min(max(proj,0),1);

    edges = linspace(0,1,numBins+1);
    [~,~,binIdx] = histcounts(proj01, edges);

    T.Sample = repmat(string(sampleName), height(T), 1);
    T.Proj01 = proj01;
    T.Bin    = binIdx;

    tables{k} = T;
    sampleNames{k} = sampleName;
end

%% --------------------- ANALYSIS PER METRIC ---------------------
edges = linspace(0,1,numBins+1);
binCenters = edges(1:end-1) + diff(edges)/2;

for m = 1:numel(metrics)
    metric = metrics{m};
    subOut = fullfile(outFolder, metric.field);
    if ~exist(subOut,'dir'), mkdir(subOut); end

    allRows = table();
    perSampleSummary = table();
    sampleInfo = struct();

    for k = 1:numFiles
        T = tables{k};
        sampleName = sampleNames{k};
        allRows = [allRows; T(:, {'Sample','Id','Proj01','Bin', metric.field})]; %#ok<AGROW>

        S = summarizeBins(T, binCenters, pctList, metric.field);
        S.Sample = repmat(string(sampleName), height(S), 1);
        perSampleSummary = [perSampleSummary; S]; %#ok<AGROW>

        [~,ord] = sort(T.Proj01);
        sampleInfo(k).sampleName = sampleName; %#ok<SAGROW>
        sampleInfo(k).Proj01 = T.Proj01(ord);
        sampleInfo(k).values = T.(metric.field)(ord);

        sampleTableOrdered = table(T.Id(ord), T.Proj01(ord), T.(metric.field)(ord), ...
            'VariableNames', {'Id','Proj01', metric.field});
        writetable(sampleTableOrdered, fullfile(subOut, sprintf('heatmap_%s_ordered.csv', matlab.lang.makeValidName(sampleName))));
    end

    writetable(allRows, fullfile(subOut,'all_cells_with_projection.csv'));
    writetable(perSampleSummary, fullfile(subOut,'bin_summary_per_sample.csv'));

    overall = summarizeBins(allRows, binCenters, pctList, metric.field);
    writetable(overall, fullfile(subOut,'bin_summary_overall.csv'));

    boxLong = allRows(:, {'Bin','Sample','Id', metric.field});
    writetable(boxLong, fullfile(subOut,'boxplot_long_per_bin.csv'));

    [~, ordAll] = sort(allRows.Proj01);
    heat_all = table(allRows.Sample(ordAll), allRows.Id(ordAll), allRows.Proj01(ordAll), allRows.(metric.field)(ordAll), ...
        'VariableNames', {'Sample','Id','Proj01', metric.field});
    writetable(heat_all, fullfile(subOut,'heatmap_allcells_ordered.csv'));

    lowessX = overall.BinCenter;
    lowessY = overall.Median;
    smY = smooth(lowessX, lowessY, lowessSpan, 'lowess');
    lowessT = table(lowessX, smY, 'VariableNames', {'BinCenter','SmoothedMedian'});
    writetable(lowessT, fullfile(subOut,'lowess_smoothed.csv'));

    makePlotsAndHeatmaps(allRows, overall, sampleInfo, numBins, subOut, metric.field, metric.label);
end
end

%% --------------------- INTERACTIVE FILE + AXIS SELECTION ---------------------
function [files, axisDef] = chooseFilesAndAxes()
    [fn,pth] = uigetfile('*.csv','Select CSV files','MultiSelect','on');
    if isequal(fn,0), files={}; axisDef={}; return; end
    if ischar(fn), fn={fn}; end
    files = fullfile(pth, fn);
    axisDef = cell(numel(files),1);
    for k=1:numel(files)
        T = readtable(files{k});
        [~,basename,ext] = fileparts(files{k});
        sampName = [basename ext];
        fprintf('\n--- File %d: %s ---\n', k, sampName);
        previewTbl = T(1:min(height(T),30), {'Id','CentroidX','CentroidY'});
        disp(previewTbl);
        cryptId = input(sprintf('Enter crypt_base Id for %s: ', sampName));
        villusId = input(sprintf('Enter villus_tip Id for %s: ', sampName));
        if ~any(T.Id==cryptId), error('cryptId %d not found', cryptId); end
        if ~any(T.Id==villusId), error('villusId %d not found', villusId); end
        centroids=[T.CentroidX,T.CentroidY];
        f = figure('Name',sprintf('Centroids preview: %s',sampName),'Color','w');
        styleDefaults(f);
        scatter(centroids(:,1),centroids(:,2),12,'k','filled'); hold on;
        rC=T(T.Id==cryptId,:); rV=T(T.Id==villusId,:);
        scatter(rC.CentroidX,rC.CentroidY,80,'b','filled','MarkerEdgeColor','w');
        scatter(rV.CentroidX,rV.CentroidY,80,'r','filled','MarkerEdgeColor','w');
        legend({'cells','crypt','villus'},'Box','off','Location','best');
        axis equal; set(gca,'YDir','reverse');
        title(sampName,'FontWeight','bold');
        prismAxes(gca);
        input('Press Enter to continue...','s'); close(f);
        axisDef{k}=struct('cryptId',cryptId,'villusId',villusId);
    end
end

%% --------------------- AXIS LOOKUP ---------------------
function [baseXY, tipXY] = getAxisFromIds(T, cryptId, villusId)
    rowC = T(T.Id==cryptId,:); rowV = T(T.Id==villusId,:);
    baseXY=[rowC.CentroidX,rowC.CentroidY];
    tipXY=[rowV.CentroidX,rowV.CentroidY];
end

%% --------------------- BIN SUMMARY ---------------------
function S = summarizeBins(T, binCenters, pctList, fieldName)
    numBins = numel(binCenters);
    Mean=nan(numBins,1); SEM=nan(numBins,1); N=zeros(numBins,1);
    Pcts=nan(numBins,numel(pctList));
    for b=1:numBins
        mask=(T.Bin==b); x=T.(fieldName)(mask); x=x(~isnan(x));
        N(b)=numel(x);
        if N(b)>0
            Mean(b)=mean(x); SEM(b)=std(x)/sqrt(N(b));
            Pcts(b,:)=prctile(x,pctList);
        end
    end
    S=table((1:numBins).',binCenters(:),N(:),Mean(:),SEM(:), ...
        'VariableNames',{'Bin','BinCenter','N','Mean','SEM'});
    for i=1:numel(pctList)
        S.(sprintf('P%02d',pctList(i)))=Pcts(:,i);
    end
    if any(strcmp('P50',S.Properties.VariableNames))
        S.Median=S.P50; else S.Median=NaN(size(Mean)); end
end

%% --------------------- PLOTTING ---------------------
function makePlotsAndHeatmaps(allRows, overall, sampleInfo, numBins, outFolder, fieldName, fieldLabel)
    [~,ordAll]=sort(allRows.Proj01); valsAll=allRows.(fieldName)(ordAll);
    clim = prctile(valsAll,[5 95]);

    %% 1) Mean + SEM (pooled; Prism-like)
    f1 = figure('Color','w','Position',[200 200 1000 600]); styleDefaults(f1);
    ax1 = axes('Parent',f1); hold(ax1,'on'); box(ax1,'off');
    p = plot(ax1, overall.BinCenter, overall.Mean, 'k-', 'LineWidth', 2);
    fb = fillBetween(ax1, overall.BinCenter, overall.Mean-overall.SEM, overall.Mean+overall.SEM, [0 0 0], 0.15);
    uistack(fb,'bottom'); uistack(p,'top');
    xlabel(ax1,'Normalized position (0 = crypt, 1 = villus)','FontWeight','bold');
    ylabel(ax1,fieldLabel,'FontWeight','bold');
    title(ax1,'Pooled mean ± SEM','FontWeight','bold');
    prismAxes(ax1);
    export600(f1, fullfile(outFolder,'fig_mean_sem.png'));

    %% 2) Combined single-cell heatmap
    f2 = figure('Color','w','Position',[200 200 1200 320]); styleDefaults(f2);
    ax2 = axes('Parent',f2); hold(ax2,'on'); box(ax2,'off');
    H = repmat(valsAll.',24,1);
    imagesc(ax2, [0 1],[0 1], H);
    set(ax2,'YTick',[],'YColor','none');
    xlabel(ax2,'Normalized position (all samples)','FontWeight','bold');
    title(ax2,['Single-cell ' fieldLabel ' (combined)'],'FontWeight','bold');
    c = colorbar(ax2); ylabel(c,fieldLabel,'FontWeight','bold');
    caxis(ax2,clim); colormap(ax2,parula);
    prismAxes(ax2,true);
    export600(f2, fullfile(outFolder,'fig_heatmap_allcells.png'));

    %% 3) Per-sample panel heatmap
    ns = numel(sampleInfo);
    f3 = figure('Color','w','Position',[200 200 1200, max(300, 160+130*ns)]); styleDefaults(f3);
    for i=1:ns
        ax = subplot(ns,1,i,'Parent',f3);
        prismAxes(ax,true);
        vals_i = sampleInfo(i).values;
        if isempty(vals_i)
            imagesc(ax, [0 1],[0 1], zeros(10,1));
        else
            imagesc(ax, [0 1],[0 1], repmat(vals_i.',20,1));
        end
        set(ax,'YTick',[],'YColor','none');
        title(ax, sampleInfo(i).sampleName, 'Interpreter','none','FontWeight','bold');
        colormap(ax,parula); caxis(ax,clim);
        if i==ns
            xlabel(ax,'Normalized position (within sample)','FontWeight','bold');
        else
            set(ax,'XTickLabel',[]);
        end
    end
    try
        sgtitle(f3,'Per-sample single-cell heatmaps','FontWeight','bold');
    catch
    end
    export600(f3, fullfile(outFolder,'fig_heatmap_per_sample_panel.png'));

    %% 4) Bin-averaged intensity heatmap
    f4 = figure('Color','w','Position',[200 200 1200 320]); styleDefaults(f4);
    ax4 = axes('Parent',f4); prismAxes(ax4,true); hold(ax4,'on');
    binMeans = overall.Mean(:)';
    H2 = repmat(binMeans,20,1);
    imagesc(ax4, [0 1],[0 1], H2);
    xlabel(ax4,'Normalized position (bin centers)','FontWeight','bold');
    set(ax4,'YTick',[],'YColor','none');
    title(ax4,['Bin-averaged ' fieldLabel ' heatmap'],'FontWeight','bold');
    c = colorbar(ax4); ylabel(c,fieldLabel,'FontWeight','bold');
    caxis(ax4,clim); colormap(ax4,parula);
    export600(f4, fullfile(outFolder,'fig_heatmap_bin_means.png'));

    %% 5) Boxplot per bin
    f5 = figure('Color','w','Position',[200 200 1400 520]); styleDefaults(f5);
    ax5 = axes('Parent',f5); hold(ax5,'on'); box(ax5,'off');
    valid = ~isnan(allRows.Bin) & allRows.Bin > 0;
    boxplot(ax5, allRows.(fieldName)(valid), allRows.Bin(valid), ...
        'PlotStyle','compact', 'Labels',[]);
    set(findobj(ax5,'Tag','Box'),'LineWidth',1.2);
    set(findobj(ax5,'Tag','Median'),'LineWidth',1.2);
    set(findobj(ax5,'Tag','Whisker'),'LineWidth',1.2);
    set(findobj(ax5,'Tag','Outliers'),'MarkerSize',3);

    edges = linspace(0,1,numBins+1);
    xticksLoc = round(linspace(1,numBins,11));
    xticks(ax5, xticksLoc);
    xticklabels(ax5, arrayfun(@(x)sprintf('%.2f',edges(x)), xticksLoc,'UniformOutput',false));
    xlabel(ax5,'Normalized position (bin index/edges)','FontWeight','bold');
    ylabel(ax5,fieldLabel,'FontWeight','bold');
    title(ax5,['Per-bin ' fieldLabel ' distributions'],'FontWeight','bold');
    prismAxes(ax5);
    export600(f5, fullfile(outFolder,'fig_boxplots_per_bin.png'));
end

%% --------------------- STYLING HELPERS ---------------------
function styleDefaults(figHandle)
    set(figHandle, ...
        'DefaultAxesFontName','Arial', ...
        'DefaultAxesFontSize',12, ...
        'DefaultAxesFontWeight','bold', ...
        'DefaultAxesLineWidth',2.5, ...
        'DefaultAxesTickDir','out', ...
        'DefaultAxesTickLength',[0.015 0.015], ...
        'DefaultAxesBox','off', ...
        'DefaultLineLineWidth',2, ...
        'DefaultFigurePaperPositionMode','auto');
end

function prismAxes(ax, isHeatmap)
    if nargin<2, isHeatmap=false; end
    set(ax, 'Box','off', ...
            'TickDir','out', ...
            'LineWidth',2.5, ...
            'FontWeight','bold', ...
            'FontName','Arial', ...
            'FontSize',12, ...
            'Layer','top');
    if isHeatmap
        grid(ax,'off');
        ax.LineWidth = 1.0;
    end
end

function export600(figHandle, outPath)
    [~,~,ext] = fileparts(outPath);
    if isempty(ext), outPath = [outPath '.png']; end
    try
        exportgraphics(figHandle, outPath, 'Resolution', 600);
    catch
        print(figHandle, outPath, '-dpng', '-r600');
    end
end

function h=fillBetween(ax,x,y1,y2,colorRGB,alphaVal)
    x=x(:)'; y1=y1(:)'; y2=y2(:)';
    hold(ax,'on');
    h=fill(ax, [x fliplr(x)], [y1 fliplr(y2)], colorRGB, ...
        'EdgeColor','none', 'FaceAlpha', alphaVal);
end
