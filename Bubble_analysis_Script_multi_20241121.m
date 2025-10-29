%% Bubble analysis, run multiple images/stacks (BINARY MASKS: pores=black, background=white) — all in one cell

clc
clear all
close all

% Fixed folder instead of uigetfile
filepath = 'C:\Users\walsh\Downloads\AFM Accuracy INTERNAL\AFM 1% INTERNAL\AFM 1% [1]';
D = dir(fullfile(filepath,'*.tif'));
filename = {D.name};
N_files = numel(filename);
filepath = [filepath filesep];

clear D   % add this line so you can reuse D later as a cell array

nspacing = 1;

% ADDED: output subfolders
out_root         = filepath; % change later if you want a different parent
out_multi_graphs = fullfile(out_root, 'BA MULTI GRAPHS');
out_multi_stats  = fullfile(out_root, 'BA MULTI STATS');
if ~exist(out_multi_graphs,'dir'), mkdir(out_multi_graphs); end
if ~exist(out_multi_stats,'dir'),  mkdir(out_multi_stats);  end


g=1;
for i=1:N_files
    fname = [filepath filename{i}];

    [All_bubbles_stack, px_size, Porosity, Pore_coverage] = ...
        Bubble_analysis_from_binary_masks(fname, nspacing);

    FinalBubble_radii{g}=All_bubbles_stack;
    Porosity_all{g}=Porosity;
    Pore_Coverage_all{g}=Pore_coverage;

    FinalBubble_diameters{g}=All_bubbles_stack.*px_size.*2;
    px_size_all{g}=px_size;

    TotalNofBubbles(g)=numel(FinalBubble_diameters{g});

    Avg_Poresize(g)=mean(FinalBubble_diameters{g});
    Std_Poresize(g)=std(FinalBubble_diameters{g});
    std_avg_ratio(g)=Std_Poresize(g)/Avg_Poresize(g);
    Median_Poresize(g)=median(FinalBubble_diameters{g});

    Avg_Porosity(g)=mean(Porosity_all{g});
    Std_Porosity(g)=std(Porosity_all{g});

    Avg_Pore_coverage(g)=mean(Pore_Coverage_all{g});
    Std_Pore_coverage(g)=std(Pore_Coverage_all{g});

    % ADDED: per-image stats & diameters 
    [~, base, ~] = fileparts(filename{g});
    T_indiv = table( ...
        Avg_Poresize(g), Std_Poresize(g), Median_Poresize(g), std_avg_ratio(g), ...
        Avg_Porosity(g), Avg_Pore_coverage(g), TotalNofBubbles(g), ...
        'VariableNames', {'Mean','SD','Median','SD_Mean_ratio','Porosity','Pore_Coverage','Total_N_bubbles'} );
    writetable(T_indiv, fullfile(out_multi_stats, [base '_stats.xlsx']));
    writetable(T_indiv, fullfile(out_multi_stats, [base '_stats.csv']));
    writematrix(FinalBubble_diameters{g}(:), fullfile(out_multi_stats, [base '_diameters_um.csv']));
    

    g=g+1;
end

% print out
nspacing
filename

Avg_Poresize
Std_Poresize
std_avg_ratio
Median_Poresize

% Plot individual normalized histograms (and save each)
g=1;
for i=1:N_files
    figure
    hold on
    box on

    b{g}=hist(FinalBubble_radii{i},[1:1:ceil(max(FinalBubble_radii{i}))]);
    P{g}=b{g}./sum(b{g});
    xx{g}=0:numel(b{g})-1;
    D{g}=xx{g}.*2.*px_size_all{g};
    plot(D{g},P{g},'.k--','MarkerSize',10)
    title({filename{g}},'FontSize',10)
    xlabel({'D (\mum)'},'FontSize',14);
    ylabel({'P(D)'},'FontSize',14);
    grid on

    % ADDED: save per-image histogram figure
    [~, base, ~] = fileparts(filename{g});
    exportgraphics(gcf, fullfile(out_multi_graphs, [base '_indiv_hist_view.png']), 'Resolution', 300);
    saveas(gcf, fullfile(out_multi_graphs, [base '_indiv_hist_view.fig']));
  

    hold off
    g=g+1;
end

% Plot all results in normalized histograms (legend outside, 11 distinct colors)
figure
set(gcf,'Position',[100 100 950 500]); % wider to fit the outside legend
hold on
box on

% Define 11 distinct colors (visually distinct)
cmap = [ ...
    0.000 0.447 0.741;  % blue
    0.850 0.325 0.098;  % orange-red
    0.929 0.694 0.125;  % yellow
    0.494 0.184 0.556;  % purple
    0.466 0.674 0.188;  % green
    0.301 0.745 0.933;  % light blue
    0.635 0.078 0.184;  % dark red
    0.000 0.500 0.000;  % dark green
    0.750 0.500 0.000;  % brownish orange
    0.250 0.250 0.250;  % dark gray
    0.900 0.600 0.900]; % light magenta

g = 1;
for i = 1:N_files
    % Pick color cyclically from the 11-color palette
    color_i = cmap(mod(i-1, size(cmap,1)) + 1, :);

    b{g}  = hist(FinalBubble_radii{i}, [1:1:ceil(max(FinalBubble_radii{i}))]);
    P{g}  = b{g} ./ sum(b{g});
    xx{g} = 0:numel(b{g})-1;
    D{g}  = xx{g} .* 2 .* px_size_all{g};

    plot(D{g}, P{g}, 'MarkerSize', 10, 'Marker', '.', 'LineStyle', '--', 'Color', color_i);
    leg(g) = {filename{g}}; %#ok<AGROW>
    g = g + 1;
end

% Legend outside on the right
lgd = legend(leg, 'Location', 'eastoutside', 'Interpreter', 'none');
xlabel({'D (\mum)'}, 'FontSize', 14);
ylabel({'P(D)'}, 'FontSize', 14);
grid on
hold off

% ADDED: save summary line plot 
exportgraphics(gcf, fullfile(out_multi_graphs, 'ALL_normalized_histograms.png'), 'Resolution', 300);
saveas(gcf, fullfile(out_multi_graphs, 'ALL_normalized_histograms.fig'))

% Tidy layout
set(gca, 'LooseInset', get(gca, 'TightInset'));

figure
hold on
box on
if exist('boxplotGroup','file') == 2
    % Use boxplotGroup with filename labels
    boxplotGroup(FinalBubble_diameters, 'PrimaryLabels', filename);
else
    % Fallback: grouped boxplot with labels = filenames
    allD = []; grp = [];
    for ii = 1:N_files
        di = FinalBubble_diameters{ii}(:);
        allD = [allD; di];
        grp = [grp; repmat(ii, numel(di), 1)];
    end
    boxplot(allD, grp, 'Labels', filename, 'LabelOrientation', 'inline');
end
ylabel({'D (\mum)'},'FontSize',14);
set(gca,'TickLabelInterpreter','none'); % show raw filenames
xtickangle(45);
hold off

% ADDED: save grouped boxplot
exportgraphics(gcf, fullfile(out_multi_graphs, 'ALL_boxplots.png'), 'Resolution', 300);
saveas(gcf, fullfile(out_multi_graphs, 'ALL_boxplots.fig'));


%%%%%% Results array (Rows: avg, std, median, sd/avg, Porosity, Pore_Coverage, Total N)
ResultsArray=[Avg_Poresize; Std_Poresize; Median_Poresize; ...
              std_avg_ratio; Avg_Porosity; Avg_Pore_coverage; TotalNofBubbles]';
ResultsTable=array2table(ResultsArray);
ResultsTable.Properties.VariableNames = ...
    {'Mean','SD','Median','SD_Mean_ratio','Porosity','Pore_Coverage','Total_N_bubbles'};
ResultsTable.Properties.RowNames = filename;

% CHANGED: save to BA MULTI STATS 
writetable(ResultsTable, fullfile(out_multi_stats, 'ResultsTable.xlsx'), 'WriteRowNames', true);
writetable(ResultsTable, fullfile(out_multi_stats, 'ResultsTable.csv'),  'WriteRowNames', true);
 

% Helper functions 
function [radii_eff_px, px_size_eff, Porosity, Pore_coverage] = Bubble_analysis_from_binary_masks(fname, nspacing)
    info = imfinfo(fname);
    radii_eff_px = [];
    Porosity = [];
    Pore_coverage = [];

    [W_um, H_um, px_size_from_tags_ok] = fov_from_metadata_or_name(fname);

    for z = 1:nspacing:numel(info)
        I = imread(fname, z);
        [H_px, W_px, ~] = size(I);

        if px_size_from_tags_ok
            [px_x_um_tag, px_y_um_tag] = pixel_size_from_tags(info(1));
            if ~isnan(px_x_um_tag) && ~isnan(px_y_um_tag)
                W_um = px_x_um_tag * W_px;
                H_um = px_y_um_tag * H_px;
            end
        end

        px_x_um = W_um / W_px;
        px_y_um = H_um / H_px;
        px_size_eff = sqrt(px_x_um * px_y_um);

        BW_pores = (I == 0);      % pores = black
        BW_bg    = ~BW_pores;     % background = not black

        CC = bwconncomp(BW_pores, 8);
        if CC.NumObjects > 0
            S = regionprops(CC,'Area');
            Apx = [S.Area]';
            r_eff = sqrt(Apx .* (px_x_um * px_y_um) / pi) ./ px_size_eff;
            radii_eff_px = [radii_eff_px; r_eff];
            area_pore = sum(Apx);
        else
            area_pore = 0;
        end

        Porosity(end+1) = area_pore / numel(BW_pores);
        Pore_coverage(end+1) = Porosity(end);
    end
end

function [px_x_um, px_y_um] = pixel_size_from_tags(tag)
    px_x_um = NaN; px_y_um = NaN;
    try
        if isfield(tag,'XResolution') && isfield(tag,'YResolution') && isfield(tag,'ResolutionUnit')
            xr = tag.XResolution; yr = tag.YResolution; unit = tag.ResolutionUnit;
            if ischar(unit); unit = string(unit); end
            if isstruct(xr); xr = xr.Value; end
            if isstruct(yr); yr = yr.Value; end
            if isnumeric(xr) && numel(xr)==2, xr = xr(1)/xr(2); end
            if isnumeric(yr) && numel(yr)==2, yr = yr(1)/yr(2); end
            if strcmpi(unit,'Inch')
                um_per_unit = 25400;
            elseif strcmpi(unit,'Centimeter')
                um_per_unit = 10000;
            else
                um_per_unit = NaN;
            end
            if ~isnan(um_per_unit)
                px_x_um = um_per_unit / xr;
                px_y_um = um_per_unit / yr;
            end
        end
    catch
        px_x_um = NaN; px_y_um = NaN;
    end
end

function [W_um, H_um, had_tags] = fov_from_metadata_or_name(fname)
    had_tags = false;
    info = imfinfo(fname);
    if ~isempty(info)
        [px_x_um_tag, px_y_um_tag] = pixel_size_from_tags(info(1));
        if ~isnan(px_x_um_tag) && ~isnan(px_y_um_tag)
            had_tags = true;
            I0 = imread(fname, 1);
            [H_px0, W_px0, ~] = size(I0);
            W_um = px_x_um_tag * W_px0;
            H_um = px_y_um_tag * H_px0;
            return
        end
    end

    [~, base, ext] = fileparts(fname);
    key = [base ext];
    W_um = 1000; H_um = 1000;
    switch key
        case '1percent.tif'
            W_um = 10.00; H_um = 10.00;
        case 'Agarose_type_A_1percent.tif'
            W_um = 10.00; H_um = 10.00;
        case '1percent_region1.tif'
            W_um = 4.00;  H_um = 2.75;
        case '1_percent_region_2.tif'
            W_um = 4.00;  H_um = 2.00;
        case 'Agarose_type_B_1percent.tif'
            W_um = 10.00; H_um = 10.02;
    end
end