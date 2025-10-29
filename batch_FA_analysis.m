%% FA analysis for ALL specified folders (exact MATLAB figs + overlays + full stats)
% Requires: createFit_Single_exp.m, createFit_double_exp.m on the MATLAB path.
% Optional: fitDwelltimehist_poremod2022.m for the Weitz model.

% ABSOLUTE OUTPUT TARGET (everything goes here)
githubRoot = 'C:\Users\walsh\Documents\GitHub\AGAROSE-HYDROGEL-TRENDS-USING-AI-ML';

% Source folders containing .tif images you want to analyze
targetDirs = { ...
    'C:\Users\walsh\Downloads\CRYO-SEM Accuracy INTERNAL\CRYO-SEM X10000\CRYO-SEM X10000 [1]' ...
};

% Safety: make sure githubRoot exists
if ~exist(githubRoot,'dir')
    error('Output folder does not exist: %s', githubRoot);
end

% Make figures render off-screen
oldVis = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','off');
cleanupObj = onCleanup(@() set(0,'DefaultFigureVisible',oldVis)); %#ok<NASGU>

for r = 1:numel(targetDirs)
    thisDir = targetDirs{r};
    if ~isfolder(thisDir)
        warning('Source folder not found: %s', thisDir);
        continue;
    end

    datasetLabel = getLastFolderName(thisDir); % e.g. "CRYO-SEM X10000 [1]"
    dataTag = sanitizeFigName(datasetLabel);   % safe for filenames

    tifList = dir(fullfile(thisDir,'*.tif'));
    if isempty(tifList)
        fprintf('[%s] No .tif files.\n', thisDir);
        continue;
    end

    pooled_void_um = [];
    pooled_pxsize  = [];
    pooled_names   = strings(0,1);
    allRows        = table();  % start as table, not []

    fprintf('\nAnalyzing dataset "%s"\n', datasetLabel);
    fprintf('Source images:   %s\n', thisDir);
    fprintf('Saving ALL outputs directly to: %s\n\n', githubRoot);

    for k = 1:numel(tifList)
        fpath = fullfile(thisDir, tifList(k).name);
        baseName  = stripExt(tifList(k).name);         % file without extension
        baseTag   = sanitizeFigName(baseName);         % safe for filenames
        prefixTag = [dataTag '_' baseTag];             % prefix for outputs

        try
            % Load image stack and compute on/off bins
            InfoImage = imfinfo(fpath);
            mImage = InfoImage(1).Width;
            nImage = InfoImage(1).Height;
            N = length(InfoImage); %#ok<NASGU>
            Image = zeros(nImage,mImage,N,'uint16');

            if isfield(InfoImage(1),'XResolution') && InfoImage(1).XResolution > 0
                px_size = 1/InfoImage(1).XResolution; % microns per pixel
            else
                error('Missing XResolution in %s', fpath);
            end

            TifLink = Tiff(fpath,'r');
            for ii = 1:N
                TifLink.setDirectory(ii);
                Image(:,:,ii) = TifLink.read();
            end
            TifLink.close();

            FinalImage = Image(:,:,1);
            FinalImage(FinalImage==255) = 1;

            RowToAnalyze = FinalImage'; % walk across x

            OnOffTh = struct('Ontimes',[],'Offtimes',[]);
            for w = 1:size(RowToAnalyze,2)
                On_diff  = diff([0; RowToAnalyze(:,w); 0] == 1);
                Off_diff = diff([1; RowToAnalyze(:,w); 1] == 0);

                OnbinsStart  = find(On_diff==1);
                OnbinsEnd    = find(On_diff==-1);
                OffbinsStart = find(Off_diff==1);
                OffbinsEnd   = find(Off_diff==-1);

                OnOffTh(w).Ontimes  = OnbinsEnd  - OnbinsStart;
                OnOffTh(w).Offtimes = OffbinsEnd - OffbinsStart;
            end

            All_Voidspace = cell2mat({OnOffTh.Offtimes}');  % pixels
            binsize = 5; % pixels per bin

            if isempty(All_Voidspace) || max(All_Voidspace) < binsize
                warning('  %s: not enough Off events to build a histogram (skipping).', tifList(k).name);
                continue;
            end

            OffHist = hist(All_Voidspace, binsize:binsize:max(All_Voidspace)); %#ok<HIST>
            OffHist = OffHist(:);
            OffHist_norm = OffHist ./ max(1,numel(All_Voidspace));

            xaxis_px        = (binsize:binsize:max(All_Voidspace))';
            x_axis_units    = xaxis_px .* px_size; % microns
            mean_Off_hist_units = mean(All_Voidspace) .* px_size;

            % Fit single and double exponential to histogram
            preFigs = findall(0,'Type','figure'); % track open figs before fitting

            [fitresult_single, gof_single] = createFit_Single_exp(x_axis_units, OffHist_norm);
            [fitresult_double, gof_double] = createFit_double_exp(x_axis_units, OffHist_norm);

            newFigs = setdiff(findall(0,'Type','figure'), preFigs);
            for hf = reshape(newFigs,1,[])
                try
                    figName = get(hf,'Name');
                    if isempty(figName), figName = 'Figure'; end
                    safeName = sprintf('%s_%s', prefixTag, sanitizeFigName(figName));
                    print(hf, fullfile(githubRoot, [safeName '.png']), '-dpng','-r300');
                    savefig(hf, fullfile(githubRoot, [safeName '.fig']));
                catch MEfig
                    warning('    (figure save) %s', MEfig.message);
                end
                closeSafely(hf);
            end

            % Extract fit parameters
            p_single = coeffvalues(fitresult_single);
            poresize_single = (-1/p_single(2)); % microns

            p_double = coeffvalues(fitresult_double);
            pore1_double = (-1/p_double(2));    % microns
            pore2_double = (-1/p_double(4));    % microns
            PoreCombined_double = (pore1_double*pore2_double)/(pore1_double+pore2_double);

            fracAmp1 = abs(p_double(1))/(abs(p_double(1))+abs(p_double(3)));
            fracAmp2 = abs(p_double(3))/(abs(p_double(1))+abs(p_double(3)));

            % Optional Weitz model fit
            Poresize_Weitz = NaN;
            cu_curve = [];
            if exist('fitDwelltimehist_poremod2022','file') == 2
                st = 1;
                st_end = numel(OffHist_norm);
                [p_weitz, ~, cu_curve, ~] = fitDwelltimehist_poremod2022( ...
                    x_axis_units(st:st_end), OffHist_norm(st:st_end), 5);
                Poresize_Weitz = p_weitz(2); % microns
            end

            % Overlay figure: data + fits
            fOverlay = figure('Visible','off','Color','w');
            ax = axes('Parent',fOverlay,'YScale','log','YMinorTick','on','FontSize',12);
            hold(ax,'on');
            semilogy(x_axis_units, OffHist_norm, '.', 'DisplayName','Data');
            hS = plot(fitresult_single, x_axis_units, OffHist_norm);
            set(hS,'DisplayName','Single exponential fit');
            hD = plot(fitresult_double, x_axis_units, OffHist_norm);
            set(hD,'DisplayName','Double exponential fit');
            if ~isempty(cu_curve)
                semilogy(x_axis_units(1:numel(cu_curve{1})), cu_curve{1}, '-', ...
                    'DisplayName','Weitz model fit');
            end
            xlabel('\xi [\mum]','FontSize',14);
            ylabel('P(\xi)','FontSize',14);
            title(sprintf('%s %s', datasetLabel, baseName), ...
                  'Interpreter','none','FontSize',12,'FontWeight','bold');
            legend('Location','northeastoutside');
            grid on;
            hold(ax,'off');
            print(fOverlay, fullfile(githubRoot, [prefixTag '_HIST_OVERLAY.png']), '-dpng','-r300');
            savefig(fOverlay, fullfile(githubRoot, [prefixTag '_HIST_OVERLAY.fig']));
            closeSafely(fOverlay);

            % Residuals figure for double exponential fit
            y_fit_double = feval(fitresult_double, x_axis_units);
            residuals = OffHist_norm - y_fit_double;

            fRes = figure('Visible','off','Color','w');
            ax2 = axes('Parent',fRes,'FontSize',12);
            hold(ax2,'on');
            plot(x_axis_units, residuals, '-', 'DisplayName','Residuals (Data - Double exp)');
            yline(0,'k:','DisplayName','Zero line');
            xlabel('\xi [\mum]','FontSize',14);
            ylabel('Residual','FontSize',14);
            title(sprintf('%s %s (Residuals)', datasetLabel, baseName), ...
                  'Interpreter','none','FontSize',12,'FontWeight','bold');
            legend('Location','northeastoutside');
            grid on;
            hold(ax2,'off');
            print(fRes, fullfile(githubRoot, [prefixTag '_RESIDUALS.png']), '-dpng','-r300');
            savefig(fRes, fullfile(githubRoot, [prefixTag '_RESIDUALS.fig']));
            closeSafely(fRes);

            % Build per-file stats row with consistent types
            stats_struct = struct( ...
                'dataset_label',      string(datasetLabel), ...
                'file_label',         string(baseName), ...
                'filename_full',      string(fpath), ...
                'px_size_um',         px_size, ...
                'mean_Off_hist_um',   mean_Off_hist_units, ...
                'poresize_single_um', poresize_single, ...
                'pore1_double_um',    pore1_double, ...
                'pore2_double_um',    pore2_double, ...
                'PoreCombined_double_um', PoreCombined_double, ...
                'fracAmp1',           fracAmp1, ...
                'fracAmp2',           fracAmp2, ...
                'Poresize_Weitz_um',  Poresize_Weitz, ...
                'gof_single_R2',      safeField(gof_single,'rsquare'), ...
                'gof_double_R2',      safeField(gof_double,'rsquare'), ...
                'n_voids',            numel(All_Voidspace) ...
            );

            T = struct2table(stats_struct);

            % save per-image stats (csv + mat)
            writetable(T, fullfile(githubRoot, [prefixTag '_stats.csv']));
            save(fullfile(githubRoot, [prefixTag '_stats.mat']), 'stats_struct');

            % save per-image histogram
            H = table(x_axis_units(:), OffHist_norm(:), ...
                'VariableNames',{'xi_um','P_xi'});
            writetable(H, fullfile(githubRoot, [prefixTag '_histogram.csv']));

            % Update pooled summary arrays for this dataset
            tmp_um = All_Voidspace .* px_size; % convert px to microns
            pooled_void_um        = [pooled_void_um; tmp_um(:)];
            pooled_pxsize(end+1,1) = px_size;
            pooled_names(end+1,1)  = string(tifList(k).name);

            % Append this row to allRows table (vertcat-safe now)
            allRows = [allRows; T]; %#ok<AGROW>

            fprintf('  %s -> exported to %s\n', tifList(k).name, githubRoot);

        catch ME
            warning('  %s failed: %s', tifList(k).name, ME.message);
        end
    end

    % Save combined per-image stats summary for this dataset
    if ~isempty(allRows)
        writetable(allRows, fullfile(githubRoot, [dataTag '_ALL_INDIVIDUAL_STATS.csv']));
    end

    % Dataset-level pooled analysis
    if ~isempty(pooled_void_um)
        try
            if numel(unique(round(pooled_pxsize,6))) > 1
                warning('[%s] Mixed pixel sizes; pooling in microns with median-based bins.', datasetLabel);
            end

            binsize_um = median(pooled_pxsize) * 5;  % 5 px equivalent, in microns
            max_um     = max(pooled_void_um);
            edges      = (binsize_um:binsize_um:max_um);
            if isempty(edges)
                edges = binsize_um;
            end

            [counts, centers] = hist(pooled_void_um, edges); %#ok<HIST>
            OffHist_norm_m = counts(:) ./ max(1,numel(pooled_void_um));
            x_axis_units_m = centers(:); % microns

            [fit_single_m, gof_single_m] = createFit_Single_exp(x_axis_units_m, OffHist_norm_m); %#ok<NASGU>
            [fit_double_m, gof_double_m] = createFit_double_exp(x_axis_units_m, OffHist_norm_m); %#ok<NASGU>

            p_single_m = coeffvalues(fit_single_m);
            p_double_m = coeffvalues(fit_double_m);

            pore_single_um   = (-1/p_single_m(2));
            pore1_double_um  = (-1/p_double_m(2));
            pore2_double_um  = (-1/p_double_m(4));
            pore_combined_um = (pore1_double_um*pore2_double_um)/(pore1_double_um+pore2_double_um);

            fracAmp1_m = abs(p_double_m(1))/(abs(p_double_m(1))+abs(p_double_m(3)));
            fracAmp2_m = abs(p_double_m(3))/(abs(p_double_m(1))+abs(p_double_m(3)));

            pore_weitz_um = NaN;
            cu_multi = [];
            if exist('fitDwelltimehist_poremod2022','file') == 2
                [p_weitz,~,cu_multi,~] = fitDwelltimehist_poremod2022(x_axis_units_m, OffHist_norm_m, 5); %#ok<ASGLU>
                pore_weitz_um = p_weitz(2);
            end

            % pooled figure (dataset-level)
            fpool = figure('Visible','off','Color','w');
            axm = axes('Parent',fpool,'YScale','log','YMinorTick','on','FontSize',12);
            hold(axm,'on');
            semilogy(x_axis_units_m, OffHist_norm_m, '.', 'DisplayName','Data');
            pmS = plot(fit_single_m, x_axis_units_m, OffHist_norm_m);
            set(pmS,'DisplayName','Single exponential fit');
            pmD = plot(fit_double_m, x_axis_units_m, OffHist_norm_m);
            set(pmD,'DisplayName','Double exponential fit');
            if ~isempty(cu_multi)
                semilogy(x_axis_units_m(1:numel(cu_multi{1})), cu_multi{1}, '-', ...
                    'DisplayName','Weitz model fit');
            end
            xlabel('\xi [\mum]','FontSize',14);
            ylabel('P(\xi)','FontSize',14);
            title(datasetLabel,'Interpreter','none','FontSize',12,'FontWeight','bold');
            legend('Location','northeastoutside');
            grid on;
            hold(axm,'off');

            print(fpool, fullfile(githubRoot, [dataTag '_MULTI_POOLED.png']), '-dpng','-r300');
            savefig(fpool, fullfile(githubRoot, [dataTag '_MULTI_POOLED.fig']));
            closeSafely(fpool);
            closeSafely(fit_single_m);
            closeSafely(fit_double_m);

            % pooled stats table
            multiStats = table( ...
                pore_single_um, pore1_double_um, pore2_double_um, pore_combined_um, ...
                fracAmp1_m, fracAmp2_m, pore_weitz_um, ...
                sum(counts)', mean(pooled_void_um), max(pooled_void_um), ...
                'VariableNames', { ...
                    'Pore_single_um','Pore1_double_um','Pore2_double_um','Pore_combined_um', ...
                    'FracAmp1','FracAmp2','Pore_Weitz_um','TotalCounts','MeanVoid_um','MaxVoid_um' ...
                } ...
            );
            writetable(multiStats, fullfile(githubRoot, [dataTag '_MULTI_STATS.csv']));
            save(fullfile(githubRoot, [dataTag '_MULTI_STATS.mat']), 'multiStats', 'pooled_names');

            % pooled histogram export
            Hmulti = table(x_axis_units_m, OffHist_norm_m, ...
                'VariableNames',{'xi_um','P_xi'});
            writetable(Hmulti, fullfile(githubRoot, [dataTag '_MULTI_histogram.csv']));

            fprintf('  pooled summary exported to %s\n', githubRoot);

        catch ME
            warning('[%s] pooled analysis failed: %s', datasetLabel, ME.message);
        end
    end
end

fprintf('\nALL OUTPUTS SAVED TO: %s\n', githubRoot);
disp('Done. Figures and stats saved to GitHub root.');

%% helper functions

function nm = stripExt(fn)
    [~, nm, ~] = fileparts(fn);
end

function name = getLastFolderName(p)
    [~, name] = fileparts(p);
end

function v = safeField(s, f)
    if isstruct(s) && isfield(s, f)
        v = s.(f);
    else
        v = NaN;
    end
end

function closeSafely(h)
    if isempty(h), return; end
    if isgraphics(h)
        try
            close(h);
        catch
        end
    elseif isvector(h)
        for ii = 1:numel(h)
            if isgraphics(h(ii))
                try
                    close(h(ii));
                catch
                end
            end
        end
    end
end

function s = sanitizeFigName(nm)
    if isempty(nm), nm = 'Figure'; end
    s = regexprep(nm,'[^A-Za-z0-9_\- ]',''); % keep alnum, space, _ , -
    s = strtrim(strrep(s,' ','_'));          % spaces -> _
end
