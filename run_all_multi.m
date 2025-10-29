%% MULTI analysis — CRYO-SEM X10000 (robust save version)

clc
clear all
close all

folderPath = 'C:\Users\walsh\Downloads\CRYO-SEM Accuracy INTERNAL\CRYO-SEM X10000\CRYO-SEM X10000 [1]';
fprintf('\n=== MULTI: %s ===\n', folderPath);

D = dir(fullfile(folderPath,'*.tif'));
filename = {D.name};
N_files  = numel(filename);
folderPath = [folderPath filesep];
clear D

githubOutRoot_multi = 'C:\Users\walsh\Documents\GitHub\AGAROSE-HYDROGEL-TRENDS-USING-AI-ML\results_multi';
if ~exist(githubOutRoot_multi,'dir'), mkdir(githubOutRoot_multi); end

out_multi_graphs = fullfile(githubOutRoot_multi, 'BA_MULTI_GRAPHs');
out_multi_stats  = fullfile(githubOutRoot_multi, 'BA_MULTI_STATs');
if ~exist(out_multi_graphs,'dir'), mkdir(out_multi_graphs); end
if ~exist(out_multi_stats,'dir'),  mkdir(out_multi_stats);  end


if N_files == 0
    warning('No .tif files found in: %s', folderPath);
end

nspacing = 1;

% Containers
FinalBubble_radii     = cell(1, N_files);
Porosity_all          = cell(1, N_files);
Pore_Coverage_all     = cell(1, N_files);
FinalBubble_diameters = cell(1, N_files);
px_size_all           = cell(1, N_files);

Avg_Poresize      = NaN(1, N_files);
Std_Poresize      = NaN(1, N_files);
std_avg_ratio     = NaN(1, N_files);
Median_Poresize   = NaN(1, N_files);
Avg_Porosity      = NaN(1, N_files);
Std_Porosity      = NaN(1, N_files);
Avg_Pore_coverage = NaN(1, N_files);
Std_Pore_coverage = NaN(1, N_files);
TotalNofBubbles   = zeros(1, N_files);

% Run per file
for g = 1:N_files
    fname = [folderPath filename{g}];
    [~, base, ~] = fileparts(filename{g});
    fprintf('   >>> MULTI processing: %s\n', filename{g});

    [All_bubbles_stack, px_size, Porosity, Pore_coverage] = ...
        Bubble_analysis_stack_v2_1mod(fname, nspacing);

    % Store raw outputs
    FinalBubble_radii{g}     = All_bubbles_stack;
    Porosity_all{g}          = Porosity;
    Pore_Coverage_all{g}     = Pore_coverage;
    FinalBubble_diameters{g} = All_bubbles_stack .* px_size .* 2;   % µm
    px_size_all{g}           = px_size;

    % Check if we actually got detections
    detected_any = ~isempty(FinalBubble_radii{g}) && any(isfinite(FinalBubble_radii{g}));

    if detected_any
        Avg_Poresize(g)      = mean(FinalBubble_diameters{g});
        Std_Poresize(g)      = std(FinalBubble_diameters{g});
        std_avg_ratio(g)     = Std_Poresize(g) / Avg_Poresize(g);
        Median_Poresize(g)   = median(FinalBubble_diameters{g});
        TotalNofBubbles(g)   = numel(FinalBubble_diameters{g});
    else
        Avg_Poresize(g)      = NaN;
        Std_Poresize(g)      = NaN;
        std_avg_ratio(g)     = NaN;
        Median_Poresize(g)   = NaN;
        TotalNofBubbles(g)   = 0;
    end

    if ~isempty(Porosity_all{g})
        Avg_Porosity(g)      = mean(Porosity_all{g});
        Std_Porosity(g)      = std(Porosity_all{g});
    else
        Avg_Porosity(g)      = NaN;
        Std_Porosity(g)      = NaN;
    end

    if ~isempty(Pore_Coverage_all{g})
        Avg_Pore_coverage(g) = mean(Pore_Coverage_all{g});
        Std_Pore_coverage(g) = std(Pore_Coverage_all{g});
    else
        Avg_Pore_coverage(g) = NaN;
        Std_Pore_coverage(g) = NaN;
    end

    % Per-image stats table
    T_indiv = table( ...
        Avg_Poresize(g), Std_Poresize(g), Median_Poresize(g), std_avg_ratio(g), ...
        Avg_Porosity(g), Avg_Pore_coverage(g), TotalNofBubbles(g), ...
        'VariableNames', {'Mean','SD','Median','SD_Mean_ratio','Porosity','Pore_Coverage','Total_N_bubbles'} );

    writetable(T_indiv, fullfile(out_multi_stats, [base '_stats.xlsx']));
    writetable(T_indiv, fullfile(out_multi_stats, [base '_stats.csv']));

    if detected_any
        writematrix(FinalBubble_diameters{g}(:), fullfile(out_multi_stats, [base '_diameters_um.csv']));
    else
        writematrix([], fullfile(out_multi_stats, [base '_diameters_um.csv']));
    end

    fprintf('      wrote stats for %s\n', base);
end

% === Summary normalized histograms across all images ===
fig_hist = figure('Position', [100 100 950 500]); hold on; box on
cmap = [ ...
    0.000 0.447 0.741;
    0.850 0.325 0.098;
    0.929 0.694 0.125;
    0.494 0.184 0.556;
    0.466 0.674 0.188;
    0.301 0.745 0.933;
    0.635 0.078 0.184;
    0.000 0.500 0.000;
    0.750 0.500 0.000;
    0.250 0.250 0.250;
    0.900 0.600 0.900];

leg = cell(1, N_files);

for i = 1:N_files
    radii_i = FinalBubble_radii{i};
    if isempty(radii_i) || all(~isfinite(radii_i))
        % draw placeholder so you still get a legend entry
        plot(nan, nan, '.', 'LineStyle', '--', 'MarkerSize', 10, ...
            'Color', cmap(mod(i-1, size(cmap,1)) + 1, :));
    else
        b = hist(radii_i, [1:1:ceil(max(radii_i))]);
        P = b ./ max(1,sum(b));
        xx = 0:numel(b)-1;
        D  = xx .* 2 .* px_size_all{i};

        color_i = cmap(mod(i-1, size(cmap,1)) + 1, :);
        plot(D, P, 'MarkerSize', 10, 'Marker', '.', 'LineStyle', '--', 'Color', color_i);
    end
    leg{i} = filename{i};
end

legend(leg, 'Location', 'eastoutside', 'Interpreter', 'none');
xlabel({'D (\mum)'}, 'FontSize', 14);
ylabel({'P(D)'}, 'FontSize', 14);
grid on; hold off
set(gca, 'LooseInset', get(gca, 'TightInset'));

exportgraphics(fig_hist, fullfile(out_multi_graphs, 'ALL_normalized_histograms.png'), 'Resolution', 300);
close(fig_hist);
fprintf('   saved ALL_normalized_histograms.png\n');

% === Grouped boxplot across all files ===
fig_box = figure; hold on; box on

% Check if boxplotGroup is available
if exist('boxplotGroup','file') == 2
    % Give empty cells something harmless so boxplotGroup doesn't choke
    FinalBubble_diameters_fixed = FinalBubble_diameters;
    for ii = 1:numel(FinalBubble_diameters_fixed)
        if isempty(FinalBubble_diameters_fixed{ii})
            FinalBubble_diameters_fixed{ii} = NaN; %#ok<AGROW>
        end
    end
    boxplotGroup(FinalBubble_diameters_fixed, 'PrimaryLabels', filename);
else
    % Fallback manual grouped boxplot
    allD = []; grp = [];
    for ii = 1:N_files
        di = FinalBubble_diameters{ii}(:);
        if isempty(di)
            di = NaN; % placeholder so the group appears
        end
        allD = [allD; di];
        grp  = [grp; repmat(ii, numel(di), 1)];
    end
    boxplot(allD, grp, 'Labels', filename, 'LabelOrientation', 'inline');
end

ylabel({'D (\mum)'},'FontSize',14);
set(gca,'TickLabelInterpreter','none');
xtickangle(35);
hold off

exportgraphics(fig_box, fullfile(out_multi_graphs, 'ALL_boxplots.png'), 'Resolution', 300);
close(fig_box);
fprintf('   saved ALL_boxplots.png\n');

% === Combined results table ===
ResultsArray=[Avg_Poresize; Std_Poresize; Median_Poresize; ...
              std_avg_ratio; Avg_Porosity; Avg_Pore_coverage; TotalNofBubbles]';

ResultsTable=array2table(ResultsArray);
ResultsTable.Properties.VariableNames = ...
    {'Mean','SD','Median','SD_Mean_ratio','Porosity','Pore_Coverage','Total_N_bubbles'};
ResultsTable.Properties.RowNames = filename;

writetable(ResultsTable, fullfile(out_multi_stats, 'ResultsTable.xlsx'), 'WriteRowNames', true);
writetable(ResultsTable, fullfile(out_multi_stats, 'ResultsTable.csv'),  'WriteRowNames', true);
fprintf('   wrote ResultsTable.xlsx / .csv\n');

fprintf('\nDone MULTI: %s\n', folderPath);
disp('All MULTI analyses finished.');
