# Auto-generated snapshot of analysis pipeline from Jupyter
# Commit this to GitHub alongside the results it produced.


OUTPUT_BASE = 'C:\\Users\\walsh\\Documents\\GitHub\\AGAROSE-HYDROGEL-TRENDS-USING-AI-ML\\cryo_sem_results'

DEFAULT_MIN_BRANCH_UM = 0.0
DEFAULT_HEAD_KEEP_PERCENT = 100
DEFAULT_MIN_RADIUS_PX = 6
MAX_CANVAS_SIDE = 1024
LARGE_N_THRESHOLD = 1500
FILE_PARAMS = {}
ROOTS = ['C:\\Users\\walsh\\Downloads\\CRYO-SEM Accuracy INTERNAL\\CRYO-SEM X10000\\CRYO-SEM X10000 [1]\\CRYO-SEM X10000 [1] STATS']



def _root_label_from_path(root_path):
    """
    Turn a long ROOT path into a short folder name for saving results in GitHub.
    Example:
    '...\\CRYO-SEM X3000 [2] STATS' -> 'CRYO-SEM_X3000_[2]_STATS'
    """
    tail = os.path.basename(root_path)  # last folder name
    label = tail.replace(" ", "_")
    return label



def exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c



def _cap_canvas(width_px, height_px, x, y):
    """
    This function keeps the working image from getting too large.

    Sometimes the input area is very big (many thousands of pixels wide).
    Big images are slow and use a lot of memory.

    What this does:
    - If the image would be bigger than MAX_CANVAS_SIDE in either direction,
      it scales the whole thing down.
    - It also scales the coordinate arrays x and y so they still line up.

    It returns:
    new_width, new_height, new_x, new_y, scale_factor
    """
    W = int(width_px)
    H = int(height_px)
    scale = 1.0

    # If either dimension is above the limit, scale everything down
    if max(W, H) > MAX_CANVAS_SIDE:
        scale = float(MAX_CANVAS_SIDE) / float(max(W, H))
        W = max(32, int(np.round(W * scale)))
        H = max(32, int(np.round(H * scale)))
        x = x * scale
        y = y * scale

    return W, H, x, y, scale



def _safe_grid_from_coords(x, y, width_px, height_px):
    """
    We have pore locations in continuous coordinates (floats).
    We want to place them onto an image grid (integer pixel locations).

    Steps:
    - Look at the span of x and y
    - Scale them so they fit inside an image of size width_px x height_px
    - Clip anything that goes slightly out of range
    - Convert to integers, because pixels are integer positions

    Returns:
    xs, ys  = pixel coordinates for each pore
    eff_x, eff_y = effective ranges used for scaling (useful for later)
    """
    if x.size:
        rng_x = float(x.max())
    else:
        rng_x = 0.0

    if y.size:
        rng_y = float(y.max())
    else:
        rng_y = 0.0

    tiny = 1e-6  # protects against divide-by-zero on very small data
    eff_x = max(rng_x, tiny)
    eff_y = max(rng_y, tiny)

    xs = (x / eff_x) * (max(1, width_px) - 1)
    ys = (y / eff_y) * (max(1, height_px) - 1)

    xs = np.clip(xs, 0, max(1, width_px) - 1).astype(int)
    ys = np.clip(ys, 0, max(1, height_px) - 1).astype(int)

    return xs, ys, eff_x, eff_y



def _knn_graph_sparse(pts, k):
    """
    Build a simple map to determine which pores and coordinates are near each other.

    For each point, we find its k nearest neighbours.
    Then record edges between those points.

    Return:
    rows, cols, data, n

    rows[i], cols[i] tell you which two points are linked.
    data[i] is the distance between those two points.
    n is how many points total.
    """
    tree = cKDTree(pts)

    # For each point, query the nearest neighbours.
    # Check for k+1 because the closest "neighbour" is always itself at distance 0.
    dists, idxs = tree.query(pts, k=min(k + 1, len(pts)))

    rows = []
    cols = []
    data = []
    n = len(pts)

    for i in range(n):
        # skip index 0 because that's the point itself
        for dist_ij, j in zip(dists[i][1:], idxs[i][1:]):
            # add edge i -> j
            rows.append(i)
            cols.append(j)
            data.append(dist_ij)

            # and add j -> i as well, to make it symmetric
            rows.append(j)
            cols.append(i)
            data.append(dist_ij)

    return rows, cols, data, n



def _mst_edges_sparse_knn(pts, k=6):
    """
    Create a minimal set of edges that keep all points connected.

    How it works:
    1. Build the k-nearest-neighbour graph (so we know which points are close).
    2. Compute a minimum spanning tree (MST) on that graph.
       The MST is basically "the cheapest set of links that still keeps
       everything in one connected network".

    We return a list of pairs (i, j) meaning "connect point i to point j".
    Each pair is sorted so (2, 5) and not (5, 2).
    """
    if len(pts) < 2:
        return []

    rows, cols, data, n = _knn_graph_sparse(pts, k=k)
    if not data:
        return []

    # Build a sparse graph where edge weights are distances
    G = coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()

    # Compute MST on that sparse graph
    T = minimum_spanning_tree(G).tocoo()

    edges = set()
    for i, j in zip(T.row, T.col):
        # Store edges in (small_index, large_index) form
        if i < j:
            a = int(i)
            b = int(j)
        else:
            a = int(j)
            b = int(i)
        edges.add((a, b))

    return list(edges)



def _rasterize_edges(pts_xy, edges, shape, thickness=1):
    """
    Turn a list of edges into a binary image.

    pts_xy is an array of point locations: [[x0, y0], [x1, y1], ...]
    edges is a list like [(0, 5), (2, 7), ...] telling us which points to connect.
    shape is (H, W), the image size.

    Draw each edge as a straight line of pixels.
    Return a 2D image with 1s where the lines are, and 0 elsewhere.

    If thickness > 1, we "fatten" the lines slightly.
    """
    H, W = shape
    canvas = np.zeros((H, W), dtype=np.uint8)

    for (i, j) in edges:
        x0, y0 = pts_xy[i]
        x1, y1 = pts_xy[j]

        rr, cc = draw_line(int(y0), int(x0), int(y1), int(x1))

        inside = (rr >= 0) & (rr < H) & (cc >= 0) & (cc < W)
        canvas[rr[inside], cc[inside]] = 1

    if thickness > 1:
        canvas = binary_dilation(canvas, disk(thickness // 2)).astype(np.uint8)

    return canvas



def build_network_fallback(xs,
                           ys,
                           base_shape,
                           upscale=2,
                           knn_k=3,
                           use_mst=True,
                           extra_dilate_px=2,
                           large_fast=False):
    """
    Build an approximate skeleton (network) from scattered pore points.
    We use this as a fallback method in two situations:
    1. When normal skeletonisation gives almost no usable branches.
    2. When there are so many points that we need a faster method.

    What this does, step by step:

    1. Take the pore coordinates (xs, ys).
    2. Put them onto an image canvas. We can optionally upscale
       to get a bit more detail.
    3. Work out which points are near each other. Connect nearby points
       with straight lines, so we get a rough network.
       Optionally, also add MST links to make sure the network is not broken
       into separate islands.
    4. Draw all those lines into an empty image.
    5. Slightly thicken the lines so there are no gaps.
    6. Run skeletonize() to thin those fat lines down to a 1-pixel-wide "wire".
       That is our final skeleton.

    Arguments:
    base_shape        : (width_pixels, height_pixels) for the starting canvas
    upscale           : how much to scale the canvas up
    knn_k             : how many neighbours each point tries to link to
    use_mst           : if True, also add MST edges for connectivity
    extra_dilate_px   : how much to thicken the lines before thinning
    large_fast        : if True, use a lighter, faster mode for huge datasets

    Returns:
    skel              : a 2D boolean array (True where skeleton is present)
    """
    Wb, Hb = base_shape

    # If in "large_fast" mode, try to keep things cheap:
    # - do not upscale
    # - use fewer neighbour links
    # - skip the MST step
    if large_fast:
        upscale = 1
        knn_k = min(2, knn_k)
        use_mst = False

    # Work out the actual canvas size we will draw on
    W = max(32, int(Wb * upscale))
    H = max(32, int(Hb * upscale))

    pts = np.column_stack([xs, ys]).astype(float)

    # If there are no points at all, just return an empty mask
    if pts.size == 0:
        return np.zeros((H, W), dtype=bool)

    # Move coordinates so the smallest x and y become 0
    px = pts[:, 0] - pts[:, 0].min()
    py = pts[:, 1] - pts[:, 1].min()

    # Scale them so they fill the canvas
    span_x = max(px.max(), 1e-6)
    span_y = max(py.max(), 1e-6)

    px = (px / span_x) * (W - 1)
    py = (py / span_y) * (H - 1)

    pts_on_canvas = np.column_stack([px, py])

    # Work out which points should be linked
    edges = []
    if len(pts_on_canvas) >= 2:
        tree = cKDTree(pts_on_canvas)
        dummy_dists, idxs = tree.query(
            pts_on_canvas,
            k=min(knn_k + 1, len(pts_on_canvas))
        )

        # Build a set of unique edges
        edge_set = set()
        for i in range(len(pts_on_canvas)):
            # idxs[i][0] is the point itself
            for j in idxs[i][1:]:
                if i < j:
                    a = i
                    b = j
                else:
                    a = j
                    b = i
                edge_set.add((a, b))

        edges = list(edge_set)

        # Optionally add MST edges for better global connectivity
        if use_mst:
            mst_edges = _mst_edges_sparse_knn(pts_on_canvas, k=max(knn_k, 6))
            edges = list(set(edges) | set(mst_edges))

    # Draw these edges into an image
    canvas = _rasterize_edges(pts_on_canvas, edges, (H, W), thickness=1)

    # Slightly thicken the lines, to make sure there are no holes
    if extra_dilate_px > 0:
        canvas = binary_dilation(canvas, disk(extra_dilate_px)).astype(np.uint8)

    # Finally, thin the lines back down to a single-pixel skeleton
    skel = skeletonize(canvas > 0).astype(bool)

    return skel



def process_csv(csv_path, out_skel_dir, out_skel_csv_dir, out_decay_dir, params):
    """
    Take one CSV file (csv_path), extract pore coordinates (X, Y), build and measure
    the pore network skeleton for that sample, and write out:
    - a skeleton preview image,
    - CSVs of pore segment lengths and summary stats,
    - and an exponential decay fit plot of those segment lengths.

    The params argument can include extra info like physical size (µm per pixel)
    and overrides for behaviour.
    """

    # Read the CSV into a table
    df = pd.read_csv(csv_path)

    # We only accept files that contain columns called "X" and "Y"
    if not {"X", "Y"}.issubset(df.columns):
        print("Skipping (no X/Y): " + csv_path)
        return

    # Pull optional metadata from the params dictionary, if present
    width_um  = float(params.get("width_um", 0))
    height_um = float(params.get("height_um", 0))
    width_px  = int(params.get("width_px", 0))
    height_px = int(params.get("height_px", 0))
    radius_override   = params.get("radius_override", None)
    MIN_BRANCH_UM     = float(params.get("min_branch_um", DEFAULT_MIN_BRANCH_UM))
    HEAD_KEEP_PERCENT = int(params.get("head_keep_percent", DEFAULT_HEAD_KEEP_PERCENT))

    # Get coordinates as an array of [X, Y] pairs (floats)
    coords = df[["X", "Y"]].astype(float).to_numpy()

    # Remove any rows that contain NaN
    coords = coords[~np.isnan(coords).any(axis=1)]

    if coords.size == 0:
        print("No valid coords in " + csv_path)
        return

    # Shift coordinates so they start at (0,0)
    x = coords[:, 0] - coords[:, 0].min()
    y = coords[:, 1] - coords[:, 1].min()

    # If the pixel size of the field of view was not provided,
    # estimate a canvas size from the data
    if width_px <= 0 or height_px <= 0:
        width_px  = max(64, int(np.ceil(max(x.max(), 1.0)))) + 1
        height_px = max(64, int(np.ceil(max(y.max(), 1.0)))) + 1

    # Reduce canvas size if it is huge (for speed)
    width_px, height_px, x, y, _ = _cap_canvas(width_px, height_px, x, y)

    # Convert from real coordinates to integer pixel locations in the canvas
    xs, ys, eff_x, eff_y = _safe_grid_from_coords(x, y, width_px, height_px)

    # Work out micrometers per pixel (um_per_px)
    # Case 1: we know the physical field of view in micrometers
    if width_um > 0 and height_um > 0:
        um_per_px_x = width_um  / float(width_px)
        um_per_px_y = height_um / float(height_px)
        um_per_px = float((um_per_px_x + um_per_px_y) / 2.0)
    else:
        # Case 2: we approximate using the spread of the coordinates themselves
        um_per_px = max(
            (eff_x / max(1, (width_px  - 1))),
            (eff_y / max(1, (height_px - 1)))
        )

    # Build a binary map (image) where each pore point sets a pixel to 1
    pore_map = np.zeros((height_px, width_px), dtype=np.uint8)
    pore_map[ys, xs] = 1

    # Decide how much to dilate the pore map before skeletonisation.
    # This roughly controls "thickness" before thinning.
    pts = np.column_stack([xs, ys])
    if radius_override is not None:
        radius = int(max(1, int(radius_override)))
    else:
        if len(pts) >= 2:
            # Look at the typical distance to nearest neighbour.
            # Use that to pick a sensible dilation radius.
            dists, _ = cKDTree(pts).query(pts, k=2)
            nn_med = float(np.median(dists[:, 1]))
            radius = max(DEFAULT_MIN_RADIUS_PX, int(np.ceil(nn_med / 2.0)))
        else:
            # If we only have one point, just fall back to a default radius
            radius = DEFAULT_MIN_RADIUS_PX

    # First skeleton attempt:
    # 1. Dilate pore_map so regions become blobs
    # 2. Skeletonize those blobs down to 1-pixel-wide traces
    pore_blobs = binary_dilation(pore_map, disk(radius))
    skeleton = skeletonize(pore_blobs).astype(bool)

    # Helper to measure skeleton branch lengths in micrometers
    def summarize_lengths(skel_bool):
        if not np.any(skel_bool):
            return [], 0

        sk = Skeleton(skel_bool)
        branch_df = summarize(sk, separator="_")

        # "summarize" sometimes names this column slightly differently
        if "branch_distance" in branch_df.columns:
            col_name = "branch_distance"
        else:
            col_name = "branch-distance"

        lengths_px = branch_df[col_name].to_numpy(dtype=float)

        # Convert from pixels to micrometers
        path_um = (lengths_px * um_per_px)

        # Filter out any: NaN, zeros, and branches shorter than MIN_BRANCH_UM
        mask_valid = np.isfinite(path_um) & (path_um > MIN_BRANCH_UM)
        path_um = path_um[mask_valid]

        return path_um.tolist(), len(path_um)

    # Get branch lengths from the first skeleton
    path_lengths_um, total_paths = summarize_lengths(skeleton)

    # If we barely got any branches, try a "fallback" method that
    # force-builds a network by linking nearby points.
    if total_paths <= 3:
        large_fast = (len(pts) >= LARGE_N_THRESHOLD)

        skel_fb = build_network_fallback(
            xs,
            ys,
            base_shape=(width_px, height_px),
            upscale=2 if not large_fast else 1,
            knn_k=3 if not large_fast else 2,
            use_mst=(not large_fast),
            extra_dilate_px=max(1, radius // 2),
            large_fast=large_fast
        )

        fb_paths_um, fb_total = summarize_lengths(skel_fb)

        # Only replace the original skeleton if the fallback found more usable paths
        if fb_total > total_paths:
            skeleton = skel_fb
            path_lengths_um = fb_paths_um
            total_paths = fb_total

    # Count junctions.
    # A "junction" is roughly where multiple branches meet.
    # The small kernel below highlights pixels that have lots of neighbours.
    branch_kernel = np.array([[1, 1, 1],
                              [1,10, 1],
                              [1, 1, 1]])
    convolved = convolve(skeleton.astype(int), branch_kernel, mode="constant")
    junctions = int(np.sum((convolved >= 13) & (skeleton == 1)))

    # Basic stats on branch lengths (in micrometers)
    if path_lengths_um:
        mean_um = float(np.mean(path_lengths_um))
        max_um  = float(np.max(path_lengths_um))
        min_um  = float(np.min(path_lengths_um))
    else:
        mean_um = 0.0
        max_um  = 0.0
        min_um  = 0.0

    # Use the filename (without .csv) for naming outputs
    prefix = os.path.splitext(os.path.basename(csv_path))[0]

    # Save an image of the skeleton network for visual QC
    plt.imshow(skeleton, cmap="gray")
    plt.title(
        "Skeleton Network\n"
        + "Junctions: " + str(junctions)
        + " | Paths: " + str(total_paths)
        + " | Mean Length: " + str(round(mean_um, 2)) + " µm"
    )
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(os.path.join(out_skel_dir, prefix + "_skeleton_network_detailed.png"), dpi=300)
    plt.close()

    # Save a CSV of all individual path lengths (in micrometers)
    pd.DataFrame({"Path_Length_um": path_lengths_um}).to_csv(
        os.path.join(out_skel_csv_dir, prefix + "_skeleton_path_lengths.csv"),
        index=False
    )

    # Save a CSV of summary metrics for this image
    pd.DataFrame({
        "Metric": [
            "Total Junctions",
            "Total Path Segments",
            "Mean Path Length (µm)",
            "Max Path Length (µm)",
            "Min Path Length (µm)",
            "µm per pixel (used)",
            "Dilation radius (px)",
            "Head kept for fit (%)",
            "Min branch used (µm)",
        ],
        "Value": [
            junctions,
            total_paths,
            round(mean_um, 3),
            round(max_um, 3),
            round(min_um, 3),
            round(float(um_per_px), 6),
            int(radius),
            int(DEFAULT_HEAD_KEEP_PERCENT),
            float(DEFAULT_MIN_BRANCH_UM),
        ],
    }).to_csv(
        os.path.join(out_skel_csv_dir, prefix + "_skeleton_summary_metrics.csv"),
        index=False
    )

    # Now try to fit an exponential decay curve to the branch lengths.
    # Idea: sort branch lengths from longest to shortest,
    # and see if they follow a "long tail" decay pattern.
    lengths = np.array(path_lengths_um, dtype=float)
    lengths = lengths[np.isfinite(lengths) & (lengths > 0)]

    if lengths.size >= 4:
        # Sort from longest to shortest
        lengths_sorted = np.sort(lengths)[::-1]

        # "ranks" is basically 0, 1, 2, ... (index in the sorted list)
        ranks = np.arange(len(lengths_sorted))

        x_fit = ranks
        y_fit = lengths_sorted

        # Rough initial guesses for the curve fit
        if len(y_fit):
            a0_guess = max(y_fit) - min(y_fit)
            c0_guess = min(y_fit)
        else:
            a0_guess = lengths_sorted[0]
            c0_guess = lengths_sorted[-1]

        # Guess a decay rate based on "half-life" around 10% of the series length
        half_life = max(1, int(0.10 * len(x_fit)))
        b0_guess = np.log(2.0) / float(half_life)

        p0 = (a0_guess, b0_guess, c0_guess)

        try:
            popt, _ = curve_fit(
                exp_decay,
                x_fit,
                y_fit,
                p0=p0,
                bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
                maxfev=20000
            )

            # Also build a simple "reference" decay curve using our guesses
            b_ref = np.log(2.0) / float(half_life)
            a_ref = lengths_sorted[0]
            c_ref = lengths_sorted[-1]
            reference_curve = exp_decay(ranks, a_ref, b_ref, c_ref)

            # Plot:
            #  - Blue points: actual sorted branch lengths
            #  - Red line: fitted exponential decay
            #  - Green dashed: simple reference curve
            fig, ax = plt.subplots(figsize=(9, 6))
            ax.scatter(ranks, lengths_sorted, label="Raw Data", s=18, color="blue")
            ax.plot(ranks, exp_decay(ranks, *popt), label="Exponential Fit", linewidth=2, color="red")
            ax.plot(ranks, reference_curve, linestyle="--", linewidth=2, label="Reference", color="green")

            ax.set_title(prefix + " — Exponential Decay of Skeleton Segment Lengths", pad=16)
            ax.set_xlabel("Segment Rank", labelpad=14)
            ax.set_ylabel("Segment Length (µm)", labelpad=14)

            ax.grid(True, alpha=0.4)
            ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0)

            plt.tight_layout(rect=[0, 0, 0.80, 1])
            plt.savefig(os.path.join(out_decay_dir, prefix + "_decay_fit.png"), dpi=300, bbox_inches="tight")
            plt.close()

        except Exception as e:
            print("Fit failed for " + csv_path + ": " + str(e))

    else:
        print(csv_path + ": not enough segments for fit (n=" + str(lengths.size) + ").")

    print("Done: " + csv_path)



def process_root(ROOT):
    if not os.path.isdir(ROOT):
        print("Missing folder: " + ROOT)
        return

    label = _root_label_from_path(ROOT)
    dest_base = os.path.join(OUTPUT_BASE, label)

    # NEW: save the current code snapshot (works in Jupyter)
    save_notebook_code(dest_base)

    out_skel_dir = os.path.join(dest_base, "SKELETON")
    out_skel_csv_dir = os.path.join(out_skel_dir, "CSV RESULTS")
    out_decay_dir = os.path.join(dest_base, "EXPONENTIAL DECAY FITTING")

    os.makedirs(out_skel_dir, exist_ok=True)
    os.makedirs(out_skel_csv_dir, exist_ok=True)
    os.makedirs(out_decay_dir, exist_ok=True)
