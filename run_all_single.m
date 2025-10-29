%% SINGLE analysis — CRYO-SEM X10000 (robust save version)

clc
clear all
close all

folderPath = 'C:\Users\walsh\Downloads\CRYO-SEM Accuracy INTERNAL\CRYO-SEM X10000\CRYO-SEM X10000 [1]';
fprintf('\n=== SINGLE: %s ===\n', folderPath);

D = dir(fullfile(folderPath,'*.tif'));
filename = {D.name};
N_files  = numel(filename);
folderPath = [folderPath filesep];
clear D

githubOutRoot = 'C:\Users\walsh\Documents\GitHub\AGAROSE-HYDROGEL-TRENDS-USING-AI-ML\results_single';
if ~exist(githubOutRoot,'dir'), mkdir(githubOutRoot); end

out_indiv_graphs = fullfile(githubOutRoot, 'BA_IND_GRAPHs');
out_indiv_stats  = fullfile(githubOutRoot, 'BA_IND_STATs');
if ~exist(out_indiv_graphs,'dir'), mkdir(out_indiv_graphs); end
if ~exist(out_indiv_stats,'dir'),  mkdir(out_indiv_stats);  end


if N_files == 0
    warning('No .tif files found in: %s', folderPath);
end

nspacing = 10;

for i = 1:N_files
    fname = [folderPath filename{i}];
    [~, base, ~] = fileparts(filename{i});
    fprintf('   >>> Processing: %s\n', filename{i});

    [All_bubbles_stack, px_size, Porosity, Pore_coverage] = ...
        Bubble_analysis_stack_v2_1mod(fname, nspacing);

    % Convert to radii/diameters (µm)
    FinalBubble_radii     = double(All_bubbles_stack);
    FinalBubble_diameters = FinalBubble_radii .* px_size .* 2;  % microns

    % Handle completely empty detections gracefully
    detected_any = ~isempty(FinalBubble_radii) && any(isfinite(FinalBubble_radii));

    % ================== Save stats even if empty ==================
    if detected_any
        Avg_Poresize      = mean(FinalBubble_diameters);
        Std_Poresize      = std(FinalBubble_diameters);
        std_avg_ratio     = Std_Poresize / Avg_Poresize;
        Median_Poresize   = median(FinalBubble_diameters);
        TotalNofBubbles   = numel(FinalBubble_diameters);
    else
        Avg_Poresize      = NaN;
        Std_Poresize      = NaN;
        std_avg_ratio     = NaN;
        Median_Poresize   = NaN;
        TotalNofBubbles   = 0;
    end

    if ~isempty(Porosity)
        Avg_Porosity      = mean(Porosity);
    else
        Avg_Porosity      = NaN;
    end

    if ~isempty(Pore_coverage)
        Avg_Pore_coverage = mean(Pore_coverage);
    else
        Avg_Pore_coverage = NaN;
    end

    T_stats = table(Avg_Poresize, Std_Poresize, Median_Poresize, std_avg_ratio, ...
                    Avg_Porosity, Avg_Pore_coverage, TotalNofBubbles, ...
                    'VariableNames', {'Mean','SD','Median','SD_Mean_ratio','Porosity','Pore_Coverage','Total_N_bubbles'});

    writetable(T_stats, fullfile(out_indiv_stats, [base '_stats.xlsx']));
    writetable(T_stats, fullfile(out_indiv_stats, [base '_stats.csv']));

    % Only write these detail CSVs if we actually had detections
    if detected_any
        writematrix(FinalBubble_diameters(:), fullfile(out_indiv_stats, [base '_diameters_um.csv']));
    else
        writematrix([], fullfile(out_indiv_stats, [base '_diameters_um.csv']));
    end

    if ~isempty(Porosity)
        writematrix(Porosity(:),              fullfile(out_indiv_stats, [base '_porosity_per_slice.csv']));
    else
        writematrix([], fullfile(out_indiv_stats, [base '_porosity_per_slice.csv']));
    end

    if ~isempty(Pore_coverage)
        writematrix(Pore_coverage(:),         fullfile(out_indiv_stats, [base '_porecoverage_per_slice.csv']));
    else
        writematrix([], fullfile(out_indiv_stats, [base '_porecoverage_per_slice.csv']));
    end

    % ================== P(D) histogram data ==================
    if detected_any
        b    = hist(FinalBubble_radii, [1:1:ceil(max(FinalBubble_radii))]);
        P    = b ./ max(1,sum(b));
        xx   = 0:numel(b)-1;
        D_um = xx .* 2 .* px_size;
        writematrix([D_um(:) P(:)], fullfile(out_indiv_stats, [base '_hist_DvP.csv']));
    else
        D_um = [];
        P    = [];
        writematrix([], fullfile(out_indiv_stats, [base '_hist_DvP.csv']));
    end

    % ================== FIGURE 1: Normalized histogram ==================
    fh_view = figure('Visible','off'); hold on; box on
    if detected_any && ~isempty(D_um)
        plot(D_um, P, '.k--', 'MarkerSize', 10)
    else
        text(0.5,0.5,'No bubbles found','HorizontalAlignment','center');
        xlim([0 1]); ylim([0 1]);
    end
    title({filename{i}}, 'FontSize', 10, 'Interpreter','none');
    xlabel({'D (\mum)'}, 'FontSize', 14);
    ylabel({'P(D)'},     'FontSize', 14);
    grid on; hold off
    exportgraphics(fh_view, fullfile(out_indiv_graphs, [base '_hist_view.png']), 'Resolution', 300);
    close(fh_view);
    fprintf('      saved: %s_hist_view.png\n', base);

    % ================== FIGURE 2: Gaussian-fit histogram ==================
    fh_fit = figure('Visible','off'); hold on; box on
    if detected_any && numel(D_um) > 1
        try
            [xData, yData] = prepareCurveData(D_um, P);
            [fitresult, gof] = fit(xData, yData, 'gauss1'); %#ok<NASGU>
            h = plot(fitresult, xData, yData);
            h(1).MarkerSize = 10;
            legend(h, 'data', 'Gaussian fit', 'Location', 'NorthEast');
        catch
            plot(D_um, P, 'k.--', 'MarkerSize', 10);
            legend('data', 'Location', 'NorthEast');
        end
    else
        text(0.5,0.5,'No fit (not enough data)','HorizontalAlignment','center');
        xlim([0 1]); ylim([0 1]);
    end
    xlabel({'D (\mum)'}, 'FontSize', 14);
    ylabel({'P(D)'},     'FontSize', 14);
    grid on; hold off
    exportgraphics(fh_fit, fullfile(out_indiv_graphs, [base '_hist_fit.png']), 'Resolution', 300);
    close(fh_fit);
    fprintf('      saved: %s_hist_fit.png\n', base);

    % ================== FIGURE 3: Boxplot ==================
    fh_box = figure('Visible','off'); hold on; box on
    if detected_any
        boxplot(FinalBubble_diameters, 'Labels', {filename{i}});
        ylabel({'D (\mum)'}, 'FontSize', 14);
        set(gca,'TickLabelInterpreter','none'); xtickangle(0);
        grid on;
    else
        text(0.5,0.5,'No diameters to boxplot','HorizontalAlignment','center');
        xlim([0 1]); ylim([0 1]);
    end
    hold off
    exportgraphics(fh_box, fullfile(out_indiv_graphs, [base '_boxplot.png']), 'Resolution', 300);
    close(fh_box);
    fprintf('      saved: %s_boxplot.png\n', base);

    % ================== FIGURE 4: whatever Bubble_analysis_stack_v2_1mod opened ==================
    figs = findall(0, 'Type', 'figure');
    k = 1;
    for f = figs'
        exportgraphics(f, fullfile(out_indiv_graphs, [base '_Fig' num2str(k) '.png']), 'Resolution', 300);
        close(f);
        k = k + 1;
    end
    fprintf('      saved: %s_Fig*.png\n', base);

    fprintf('   <<< Done: %s\n', filename{i});
end

fprintf('\nDone SINGLE: %s\n', folderPath);
disp('All SINGLE analyses finished.');
