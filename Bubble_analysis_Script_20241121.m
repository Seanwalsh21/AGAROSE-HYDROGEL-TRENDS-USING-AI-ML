%% Bubble analysis — run the single-image workflow over all .tif in a chosen folder
% Outputs:
%   - Figures per image  ->  BA INDIVIDUAL GRAPHS
%   - Stats/CSVs per image -> BA INDIVIDUAL STATS

clc
clear all
close all

% SELECT OR DEFINE FOLDER 
% Option 1: Pick folder interactively
filepath = uigetdir(pwd, 'Select folder containing .tif images');
if filepath == 0
    error('No folder selected. Script aborted.');
end

% Option 2 (manual): uncomment and replace with your path if preferred
% filepath = 'C:\path\to\your\images';

D = dir(fullfile(filepath, '*.tif'));
filename = {D.name};
N_files = numel(filename);
filepath = [filepath filesep];
clear D

% OUTPUT SUBFOLDERS (created next to input folder)
out_indiv_graphs = fullfile(filepath, 'BA INDIVIDUAL GRAPHS');
out_indiv_stats  = fullfile(filepath, 'BA INDIVIDUAL STATS');
if ~exist(out_indiv_graphs,'dir'), mkdir(out_indiv_graphs); end
if ~exist(out_indiv_stats,'dir'),  mkdir(out_indiv_stats);  end

nspacing = 10;   % same as your single-cell analysis

for i = 1:N_files
    close all  % keep figures tidy per file

    % RUN the original single-image workflow 
    fname = [filepath filename{i}];
    [All_bubbles_stack, px_size, Porosity, Pore_coverage] = Bubble_analysis_stack_v2_1mod(fname, nspacing);

    % Analyze pore size distribution
    FinalBubble_radii = double(All_bubbles_stack);
    FinalBubble_diameters = FinalBubble_radii .* px_size .* 2;  % micrometers

    % Histogram (same binning approach)
    b = hist(FinalBubble_radii, [1:1:ceil(max(FinalBubble_radii))]);
    hist_norm = b ./ sum(b);
    P = hist_norm;

    xx = 0:numel(b)-1;          % x values for histogram
    D_um = xx .* 2 .* px_size;  % diameter axis in micrometers

    % Summary numbers
    Avg_Poresize      = mean(FinalBubble_diameters);
    Std_Poresize      = std(FinalBubble_diameters);
    std_avg_ratio     = Std_Poresize / Avg_Poresize;
    Median_Poresize   = median(FinalBubble_diameters);

    % SAVE STATS/CSVs for this image 
    [~, base, ~] = fileparts(filename{i});

    % Stats table (per image)
    T_stats = table( ...
        Avg_Poresize, Std_Poresize, Median_Poresize, std_avg_ratio, ...
        mean(Porosity), mean(Pore_coverage), numel(FinalBubble_diameters), ...
        'VariableNames', {'Mean','SD','Median','SD_Mean_ratio','Porosity','Pore_Coverage','Total_N_bubbles'} ...
    );
    writetable(T_stats, fullfile(out_indiv_stats, [base '_stats.xlsx']));
    writetable(T_stats, fullfile(out_indiv_stats, [base '_stats.csv']));

    % Raw outputs
    writematrix(FinalBubble_diameters(:), fullfile(out_indiv_stats, [base '_diameters_um.csv']));
    writematrix(Porosity(:),            fullfile(out_indiv_stats, [base '_porosity_per_slice.csv']));
    writematrix(Pore_coverage(:),       fullfile(out_indiv_stats, [base '_porecoverage_per_slice.csv']));
    writematrix([D_um(:) P(:)],         fullfile(out_indiv_stats, [base '_hist_DvP.csv']));

    % FIGURE 1: Gaussian fit of normalized histogram
    try
        [xData, yData] = prepareCurveData(D_um, P);
        [fitresult, gof] = fit(xData, yData, 'gauss1'); %#ok<NASGU>
        fh1 = figure('Visible','off'); hold on; box on;
        h = plot(fitresult, xData, yData);
        h(1).MarkerSize = 10;  % data markers
        legend(h, 'data', 'Gaussian fit', 'Location', 'NorthEast');
    catch
        % If Curve Fitting Toolbox isn't available
        fh1 = figure('Visible','off'); hold on; box on;
        plot(D_um, P, 'k.--', 'MarkerSize', 10);
        legend('data', 'Location', 'NorthEast');
    end
    xlabel({'D (\mum)'}, 'FontSize', 14);
    ylabel({'P(D)'},    'FontSize', 14);
    grid on; hold off;

    exportgraphics(fh1, fullfile(out_indiv_graphs, [base '_hist_fit.png']), 'Resolution', 300);
    saveas(fh1,       fullfile(out_indiv_graphs, [base '_hist_fit.fig']));
    close(fh1);

    % FIGURE 2: Boxplot (labeled with filename) 
    fh2 = figure('Visible','off'); hold on; box on;
    boxplot(FinalBubble_diameters, 'Labels', {filename{i}});
    ylabel({'D (\mum)'},'FontSize',14);
    set(gca,'TickLabelInterpreter','none'); xtickangle(0);
    grid on; hold off;

    exportgraphics(fh2, fullfile(out_indiv_graphs, [base '_boxplot.png']), 'Resolution', 300);
    saveas(fh2,       fullfile(out_indiv_graphs, [base '_boxplot.fig']));
    close(fh2);

    % (Optional) simple data echo to console 
    fprintf('Done: %s | Mean=%.3f µm, SD=%.3f µm, Median=%.3f µm, SD/Mean=%.3f | Porosity=%.4f\n', ...
        filename{i}, Avg_Poresize, Std_Poresize, Median_Poresize, std_avg_ratio, mean(Porosity));
end

disp('All individual analyses complete.');