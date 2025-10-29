function [All_bubbles_stack, px_size, Porosity, Pore_coverage] = Bubble_analysis_stack_v2_1mod(fname, nspacing)
% Bubble analysis on binary masks
% Goal for this version:
%   - Show image with colors INVERTED (prev white -> black, prev black -> white)
%   - Detect and draw circles in the NEW BLACK regions (i.e., previously white)

    if nargin < 1
        [filename, filepath] = uigetfile('*.tif');
        fname = [filepath filename];
    end

    % --- Read stack (robust to RGB / different classes) ---
    info    = imfinfo(fname);
    mImage  = info(1).Width;      % columns (x)
    nImage  = info(1).Height;     % rows (y)
    nSlices = numel(info);

    % Preallocate as grayscale uint16
    Image = zeros(nImage, mImage, nSlices, 'uint16');

    % Pixel size [um]  (fallback-safe if XResolution missing/zero)
    if isfield(info(1),'XResolution') && ~isempty(info(1).XResolution) && info(1).XResolution > 0
        px_size = 1 / info(1).XResolution;
    else
        px_size = 1;  % fallback (avoid NaN/Inf)
    end

    tif = Tiff(fname, 'r');
    for k = 1:nSlices
        tif.setDirectory(k);
        imgk = tif.read();  % could be gray, RGB, logical, uint8, uint16, etc.

        % If RGB, convert to grayscale in native class
        if ndims(imgk) == 3 && size(imgk,3) == 3
            imgk = rgb2gray(imgk);
        end

        % Convert to uint16 without relying on Image Processing Toolbox
        switch class(imgk)
            case 'uint16'
                % ok
            case 'uint8'
                imgk = uint16(imgk) * 257;          % 0..255 -> 0..65535
            case 'logical'
                imgk = uint16(imgk) * 65535;
            case {'double','single'}
                % Scale to 0..65535 safely
                imgk = uint16(65535 * mat2gray(imgk));
            otherwise
                % Fallback for other integer types
                imgk = uint16(65535 * mat2gray(double(imgk)));
        end

        % Size sanity (some odd pages can differ); center-crop/pad if needed
        [hk, wk] = size(imgk);
        if hk ~= nImage || wk ~= mImage
            imgk = imresize(imgk, [nImage mImage], 'nearest');
        end

        Image(:,:,k) = imgk;
    end
    tif.close();

    slice_counter = 1;

    for z = 1:nspacing:nSlices

        % -------------------------------------------------
        % 1) Define masks: target pores = PREVIOUSLY WHITE
        % -------------------------------------------------
        Iraw      = Image(:,:,z);       % original binary slice
        BW_pores  = (Iraw > 0);         % previously white -> NOW the "pores" we target
        BW_solid  = ~BW_pores;          % previously black

        % Porosity over this slice (counting new "pores" = previously white)
        V_cube                  = numel(BW_pores);
        V_fluid                 = nnz(BW_pores);
        Porosity(slice_counter) = V_fluid / V_cube;

        % -------------------------------------------------
        % 2) Distance transform: peaks inside NEW BLACK (prev white)
        % -------------------------------------------------
        EDM          = bwdist(BW_solid);                % distance to solids (prev black)
        smoothed_EDM = imfilter(EDM, fspecial('gaussian',5,1));
        local_maxima = imregionalmax(smoothed_EDM);

        % Maxima to coordinates (x = col, y = row)
        [rows_max, cols_max] = find(local_maxima);
        bubble_coord = [cols_max, rows_max];            % [x, y]
        bubble_radii = EDM(local_maxima);

        n_bubble = numel(bubble_radii);

        % Parametric circles for plotting
        th = 0:pi/50:2*pi;
        xunit = zeros(n_bubble, numel(th));
        yunit = zeros(n_bubble, numel(th));
        for iB = 1:n_bubble
            xunit(iB,:) = bubble_radii(iB) * cos(th) + bubble_coord(iB,1);
            yunit(iB,:) = bubble_radii(iB) * sin(th) + bubble_coord(iB,2);
        end

        % -------------------------------------------------
        % 3) Filters (edge, tangent, overlap)
        % -------------------------------------------------
        EdgeX = round(0.05 * mImage);
        EdgeY = round(0.05 * nImage);

        EdgecoordsX = (bubble_coord(:,1) < EdgeX) | (bubble_coord(:,1) > (mImage - EdgeX));
        EdgecoordsY = (bubble_coord(:,2) < EdgeY) | (bubble_coord(:,2) > (nImage - EdgeY));
        Edge_keep   = ~(EdgecoordsX | EdgecoordsY);

        idx_keep             = find(Edge_keep);
        bubble_coord_Filter1 = bubble_coord(idx_keep,:);
        bubble_radii_Filter1 = bubble_radii(idx_keep);
        xunit_Filter1        = xunit(idx_keep,:);
        yunit_Filter1        = yunit(idx_keep,:);
        n_bubble_Filter1     = numel(bubble_radii_Filter1);

        % Tangent filtering (Molteni 2013)
        xunit_r = round(xunit_Filter1);
        yunit_r = round(yunit_Filter1);

        g = 1;
        bubble_radii_Filter1_zeros = zeros(n_bubble_Filter1,1);
        for iB = 1:n_bubble_Filter1
            % Define R1 up front so it's always available
            R1 = bubble_radii_Filter1(iB);

            % Sample the smoothed EDM along the circle
            k = 1;
            BubbleXY = [];
            for j = 1:numel(th)
                xj = xunit_r(iB,j);
                yj = yunit_r(iB,j);
                if xj >= 1 && xj <= mImage && yj >= 1 && yj <= nImage
                    BubbleXY(k) = smoothed_EDM(yj, xj);  % (row=y, col=x)
                    k = k + 1;
                end
            end
            BubbleXY = double(BubbleXY);

            % Find minima in EDM (as peaks in -EDM); guard against empty
            if ~isempty(BubbleXY)
                [PKS, ~]   = findpeaks(-BubbleXY);
                PKS_sorted = sort(-PKS, 'ascend');  % minima values
            else
                PKS_sorted = [];
            end

            % Keep decision
            if numel(PKS_sorted) >= 3
                Ravg = mean([R1, R1 + PKS_sorted(2), R1 + PKS_sorted(3)]);
                dR   = 0.05 * Ravg;
                keep = (R1 + PKS_sorted(2) <= Ravg + dR) && (R1 + PKS_sorted(3) <= Ravg + dR);
            else
                keep = true;
            end

            % Store kept radius (or 0)
            bubble_radii_Filter1_zeros(g) = R1 * double(keep);
            g = g + 1;
        end

        % Proceed to Filter2
        idx2                 = find(bubble_radii_Filter1_zeros);
        bubble_radii_Filter2 = bubble_radii_Filter1(idx2);
        bubble_coord_Filter2 = bubble_coord_Filter1(idx2,:);
        xunit_Filter2        = xunit_Filter1(idx2,:);
        yunit_Filter2        = yunit_Filter1(idx2,:);
        n_bubble_Filter2     = numel(bubble_radii_Filter2);

        % Overlap filtering (keep larger first)
        [bubble_radii_Filter2_sorted, ind_sorted] = sort(bubble_radii_Filter2, 'descend');
        bubble_coord_Filter2_sorted = bubble_coord_Filter2(ind_sorted,:);
        xunit_Filter2_sorted        = xunit_Filter2(ind_sorted,:);
        yunit_Filter2_sorted        = yunit_Filter2(ind_sorted,:);

        bubble_radii_Filter3 = bubble_radii_Filter2_sorted(1);
        bubble_coord_Filter3 = bubble_coord_Filter2_sorted(1,:);
        xunit_Filter3        = xunit_Filter2_sorted(1,:);
        yunit_Filter3        = yunit_Filter2_sorted(1,:);

        n = 1; 
        g = 2;
        for iB = 2:n_bubble_Filter2
            eD = zeros(n,1);
            for j = 1:n
                eD(j) = CalcDistance( ...
                    bubble_coord_Filter2_sorted(iB,1), bubble_coord_Filter2_sorted(iB,2), ...
                    bubble_coord_Filter3(j,1),         bubble_coord_Filter3(j,2));
            end
            if min(eD - bubble_radii_Filter3(1:n)) > 0
                bubble_radii_Filter3(g,1) = bubble_radii_Filter2_sorted(iB);
                bubble_coord_Filter3(g,:) = bubble_coord_Filter2_sorted(iB,:);
                xunit_Filter3(g,:)        = xunit_Filter2_sorted(iB,:);
                yunit_Filter3(g,:)        = yunit_Filter2_sorted(iB,:);
                g = g + 1; 
                n = n + 1;
            end
        end

        n_bubble_Filter3 = numel(bubble_radii_Filter3);

        % -------------------------------------------------
        % 4) DISPLAY: invert colors & keep circles in NEW BLACK
        % -------------------------------------------------
        figure; 
        hold on;
        imagesc(BW_solid);                % background white, NEW pores black
        colormap(gray); axis image; axis ij;
        xlim([0 mImage]); ylim([0 nImage]);
        scatter(bubble_coord_Filter3(:,1), bubble_coord_Filter3(:,2), 'x', 'MarkerEdgeColor',[0 0.5 1]);
        for iB = 1:n_bubble_Filter3
            plot(xunit_Filter3(iB,:), yunit_Filter3(iB,:), 'r', 'LineWidth', 1);
        end
        title('Final bubbles'); 
        hold off;

        % -------------------------------------------------
        % 5) Coverage using NEW pore mask (previously white)
        % -------------------------------------------------
        bubble_mask       = createCirclesMask([mImage, nImage], bubble_coord_Filter3, bubble_radii_Filter3);
        bubble_mask_units = bubble_mask * 255;

        BW_pores_crop = BW_pores(EdgeX:end-EdgeX, EdgeY:end-EdgeY);
        Area_Fluid    = nnz(BW_pores_crop);  % pore pixels (new pores)
        Area_Bubbles  = nnz(uint16(bubble_mask_units(EdgeX:end-EdgeX, EdgeY:end-EdgeY)));
        Pore_coverage(slice_counter) = Area_Bubbles / max(1, Area_Fluid);

        % Store radii
        FinalBubble_radii_STACK(slice_counter).bubble_radii_Filter3 = bubble_radii_Filter3;

        slice_counter = slice_counter + 1;

        % housekeeping
        clear rows_max cols_max xunit yunit xunit_Filter* yunit_Filter* idx_keep idx2
    end

    % Combine all slices
    c = {FinalBubble_radii_STACK.bubble_radii_Filter3};
    All_bubbles_stack = cell2mat(c');

end
