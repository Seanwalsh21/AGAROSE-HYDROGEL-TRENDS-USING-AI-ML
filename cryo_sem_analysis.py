
import os
import csv
import json
import numpy as np
from PIL import Image
from scipy.ndimage import binary_erosion

# This script:
# 1. Reads GOLD STANDARD.tif and each method .tif from BASE_DIR
# 2. Converts them to binary pore masks (1=pore, 0=background)
# 3. Calculates Dice, IoU, MCC, etc.
# 4. Measures pore fraction and pore area in µm²
# 5. Saves:
#       - per-method overlay TIFFs
#       - per-method CSVs
#       - one summary CSV for all methods
#       - a copy of this code (.py)
#       - a minimal notebook (.ipynb)
#
# All outputs go into your GitHub repo folder.

BASE_DIR = r"C:\Users\walsh\Downloads\CRYO-SEM Accuracy INTERNAL\CRYO-SEM X10000\CRYO-SEM X10000 [1]"
GT_FILENAME = "GOLD STANDARD.tif"

GITHUB_DIR = r"C:\Users\walsh\Documents\GitHub\AGAROSE-HYDROGEL-TRENDS-USING-AI-ML"
os.makedirs(GITHUB_DIR, exist_ok=True)

METRICS_CSV = os.path.join(GITHUB_DIR, "metrics_summary.csv")

# pixel calibration: 640 x 480 image covers 5.98 µm x 4.49 µm
PIX_SIZE_X_UM = 5.98 / 640.0
PIX_SIZE_Y_UM = 4.49 / 480.0
PIX_AREA_UM2  = PIX_SIZE_X_UM * PIX_SIZE_Y_UM


def _read_tiff_any(path):
    """Read a TIFF from disk using Pillow. Return as numpy array."""
    if not os.path.exists(path):
        return None
    with Image.open(path) as im:
        im.load()
        # convert to grayscale 16-bit or 8-bit
        if "I;16" in im.mode:
            im = im.convert("I;16")
        else:
            im = im.convert("L")
        arr = np.array(im)
    return arr


def _to_gray(arr):
    """Ensure we have a single 2D grayscale array."""
    if arr is None:
        return None
    if arr.ndim == 2:
        return arr
    if arr.ndim == 3:
        # if RGB slipped in, average channels
        if arr.shape[2] >= 3:
            return np.mean(arr[:, :, :3], axis=2).astype(arr.dtype)
        else:
            return arr[:, :, 0]
    return arr


def _ensure_uint(arr):
    """Force array into uint8 or uint16, for consistent thresholding."""
    if arr.dtype == np.uint8 or arr.dtype == np.uint16:
        return arr

    if np.issubdtype(arr.dtype, np.floating):
        a_min = float(arr.min())
        a_max = float(arr.max())
        rng = (a_max - a_min) + 1e-12
        scaled = (arr - a_min) / rng
        scaled = (scaled * 255.0 + 0.5).astype(np.uint8)
        return scaled

    maxv = float(arr.max())
    if maxv > 255.0:
        return arr.astype(np.uint16)
    return arr.astype(np.uint8)


def otsu_thresh_uint8(gray_u8):
    """Manual Otsu threshold. Returns mask of dark pixels as 1."""
    hist = np.bincount(gray_u8.flatten(), minlength=256).astype(float)
    total = gray_u8.size
    prob = hist / float(total)

    cum_prob = np.cumsum(prob)
    cum_mean = np.cumsum(prob * np.arange(256))
    global_mean = cum_mean[-1]

    best_t = 0
    best_score = -1.0

    for t in range(256):
        w0 = cum_prob[t]
        w1 = 1.0 - w0
        if w0 == 0.0 or w1 == 0.0:
            continue
        mu0 = cum_mean[t] / w0
        mu1 = (global_mean - cum_mean[t]) / w1
        diff = mu0 - mu1
        score = w0 * w1 * diff * diff
        if score > best_score:
            best_score = score
            best_t = t

    # pores are darker, so <= threshold is pore
    mask_dark = (gray_u8 <= best_t).astype(np.uint8)
    return mask_dark


def binarize_pores_black(img_gray):
    """Make a binary mask where pores (dark) are 1 and background is 0."""
    g = _ensure_uint(img_gray)
    uvals = np.unique(g)

    # common case where the mask is already binary grayscale
    if uvals.size == 2:
        darker = int(uvals[0])
        pores = (g == darker)
        return pores.astype(np.uint8)

    # not binary already: do Otsu on 8-bit version
    if g.dtype == np.uint16:
        g8 = (g / 257).astype(np.uint8)
    else:
        g8 = g.astype(np.uint8)

    pores = otsu_thresh_uint8(g8)
    return pores


def read_mask_as_binary(path):
    """Load tiff, turn into a binary pore mask (0/1)."""
    raw = _read_tiff_any(path)
    if raw is None:
        raise FileNotFoundError("Cannot read TIFF: " + path)
    gray = _to_gray(raw)
    if gray is None or gray.ndim != 2:
        raise ValueError("Not single-channel grayscale: " + path)
    return binarize_pores_black(gray)


def _safe_div(n, d):
    if d == 0:
        return 0.0
    return float(n) / float(d)


def calc_confusion(gt, pr):
    """Return tp, fp, tn, fn for two binary masks."""
    gt = gt.astype(np.uint8)
    pr = pr.astype(np.uint8)

    tp = int(np.sum((gt == 1) & (pr == 1)))
    tn = int(np.sum((gt == 0) & (pr == 0)))
    fp = int(np.sum((gt == 0) & (pr == 1)))
    fn = int(np.sum((gt == 1) & (pr == 0)))

    return tp, fp, tn, fn


def compute_metrics(gt, pr):
    """Return a dict of Dice, IoU, MCC, etc."""
    tp, fp, tn, fn = calc_confusion(gt, pr)

    acc  = _safe_div(tp + tn, tp + tn + fp + fn)
    prec = _safe_div(tp, tp + fp)
    rec  = _safe_div(tp, tp + fn)
    spec = _safe_div(tn, tn + fp)
    ba   = 0.5 * (rec + spec)
    dice = _safe_div(2 * tp, 2 * tp + fp + fn)
    iou  = _safe_div(tp, tp + fp + fn)

    # MCC denominator
    prod_val = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if prod_val > 0.0:
        den = prod_val ** 0.5
    else:
        den = 0.0
    if den > 0.0:
        mcc_val = (tp * tn - fp * fn) / den
    else:
        mcc_val = 0.0

    out = {
        "accuracy": acc,
        "precision": prec,
        "recall": rec,
        "specificity": spec,
        "balanced_accuracy": ba,
        "f1_dice": dice,
        "iou_jaccard": iou,
        "mcc": mcc_val,
        "TP": tp,
        "FP": fp,
        "TN": tn,
        "FN": fn
    }
    return out


def pore_fraction(mask):
    """Fraction of pixels that are pore=1."""
    return float(np.mean(mask))


def pore_area_um2(mask):
    """Total pore area in square microns."""
    pore_px = int(np.sum(mask == 1))
    return pore_px * PIX_AREA_UM2


def make_pil_overlay(gt_mask, pr_mask, save_path):
    """
    Make an overlay using Pillow .paste(mask=...),
    similar to the GeeksforGeeks example.

    Colors:
      red   = false positive (pred says pore, GT says no pore)
      blue  = false negative (pred missed pore GT has)
      green = outline of true positive
    """
    base_gray = (1 - gt_mask) * 255.0
    base_gray = base_gray.astype(np.uint8)
    base_img = Image.fromarray(base_gray, mode="L").convert("RGBA")

    h, w = gt_mask.shape

    # false positives
    fp_mask = ((gt_mask == 0) & (pr_mask == 1)).astype(np.uint8) * 255
    red_img = Image.new("RGBA", (w, h), (255, 0, 0, 180))
    red_mask = Image.fromarray(fp_mask.astype(np.uint8), mode="L")
    base_img.paste(red_img, (0, 0), mask=red_mask)

    # false negatives
    fn_mask = ((gt_mask == 1) & (pr_mask == 0)).astype(np.uint8) * 255
    blue_img = Image.new("RGBA", (w, h), (0, 0, 255, 180))
    blue_mask = Image.fromarray(fn_mask.astype(np.uint8), mode="L")
    base_img.paste(blue_img, (0, 0), mask=blue_mask)

    # true positives, outline only
    tp_region = ((gt_mask == 1) & (pr_mask == 1)).astype(np.uint8)
    tp_eroded = binary_erosion(tp_region, border_value=0)
    tp_edge = tp_region.astype(np.uint8) - tp_eroded.astype(np.uint8)
    tp_edge_mask = (tp_edge > 0).astype(np.uint8) * 255

    green_img = Image.new("RGBA", (w, h), (0, 255, 0, 255))
    green_mask = Image.fromarray(tp_edge_mask.astype(np.uint8), mode="L")
    base_img.paste(green_img, (0, 0), mask=green_mask)

    # save RGB TIFF
    final_rgb = base_img.convert("RGB")
    final_rgb.save(save_path, format="TIFF")
    return save_path


def write_notebook_copy(code_text, ipynb_out_path):
    """Save a tiny 1-cell .ipynb that just contains this script."""
    nb_obj = {
        "cells": [
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": code_text.splitlines(True)
            }
        ],
        "metadata": {
            "language_info": {"name": "python"},
            "kernelspec": {
                "display_name": "Python",
                "language": "python",
                "name": "python"
            }
        },
        "nbformat": 4,
        "nbformat_minor": 5
    }

    with open(ipynb_out_path, "w", encoding="utf-8") as f_out:
        json.dump(nb_obj, f_out, indent=2)


def write_single_method_csv(row_dict, repo_dir):
    """
    Save metrics for a single method as <method>_metrics.csv
    so each method is traceable in Git.
    """
    method_name = row_dict["method"]
    out_path = os.path.join(repo_dir, method_name + "_metrics.csv")

    cols = [
        "method",
        "f1_dice",
        "iou_jaccard",
        "mcc",
        "precision",
        "recall",
        "specificity",
        "balanced_accuracy",
        "accuracy",
        "TP",
        "FP",
        "TN",
        "FN",
        "gt_pore_fraction",
        "pred_pore_fraction",
        "pore_fraction_bias",
        "gt_pore_area_um2",
        "pred_pore_area_um2",
        "pore_area_bias_um2"
    ]

    with open(out_path, "w", newline="", encoding="utf-8") as f_one:
        w = csv.DictWriter(f_one, fieldnames=cols)
        w.writeheader()
        w.writerow(row_dict)

    return out_path


def analyze_folder(folder, repo_dir, code_text):
    """
    Do the full run:
    - read GT
    - loop over each method .tif
    - compute metrics and pore stats
    - save overlay.tif
    - save <method>_metrics.csv
    - save metrics_summary.csv
    - save cryo_sem_analysis.py and cryo_sem_analysis.ipynb
    """

    gt_path = os.path.join(folder, GT_FILENAME)
    if not os.path.exists(gt_path):
        raise FileNotFoundError("Ground truth not found: " + gt_path)

    gt_mask = read_mask_as_binary(gt_path)

    all_rows = []

    for fname in os.listdir(folder):
        low = fname.lower()
        if not (low.endswith(".tif") or low.endswith(".tiff")):
            continue
        if fname == GT_FILENAME:
            continue

        pr_path = os.path.join(folder, fname)
        pr_mask = read_mask_as_binary(pr_path)

        if gt_mask.shape != pr_mask.shape:
            raise ValueError("Shape mismatch: GT " + str(gt_mask.shape) +
                             " vs " + fname + " " + str(pr_mask.shape))

        m = compute_metrics(gt_mask, pr_mask)

        gt_frac = pore_fraction(gt_mask)
        pr_frac = pore_fraction(pr_mask)
        frac_bias = pr_frac - gt_frac

        gt_area = pore_area_um2(gt_mask)
        pr_area = pore_area_um2(pr_mask)
        area_bias = pr_area - gt_area

        method_name = os.path.splitext(fname)[0]

        overlay_file = os.path.join(repo_dir, method_name + "_overlay.tif")
        make_pil_overlay(gt_mask, pr_mask, overlay_file)

        row = {
            "method": method_name,
            "f1_dice": m["f1_dice"],
            "iou_jaccard": m["iou_jaccard"],
            "mcc": m["mcc"],
            "precision": m["precision"],
            "recall": m["recall"],
            "specificity": m["specificity"],
            "balanced_accuracy": m["balanced_accuracy"],
            "accuracy": m["accuracy"],
            "TP": m["TP"],
            "FP": m["FP"],
            "TN": m["TN"],
            "FN": m["FN"],
            "gt_pore_fraction": gt_frac,
            "pred_pore_fraction": pr_frac,
            "pore_fraction_bias": frac_bias,
            "gt_pore_area_um2": gt_area,
            "pred_pore_area_um2": pr_area,
            "pore_area_bias_um2": area_bias
        }

        all_rows.append(row)

        single_csv_path = write_single_method_csv(row, repo_dir)

        print(method_name + ": Dice=%.3f IoU=%.3f Bias=%.2f%% -> %s / %s"
              % (m["f1_dice"], m["iou_jaccard"],
                 frac_bias * 100.0, overlay_file, single_csv_path))

    # now write combined CSV for all methods
    if len(all_rows) > 0:
        fieldnames = [
            "method",
            "f1_dice",
            "iou_jaccard",
            "mcc",
            "precision",
            "recall",
            "specificity",
            "balanced_accuracy",
            "accuracy",
            "TP",
            "FP",
            "TN",
            "FN",
            "gt_pore_fraction",
            "pred_pore_fraction",
            "pore_fraction_bias",
            "gt_pore_area_um2",
            "pred_pore_area_um2",
            "pore_area_bias_um2"
        ]
    else:
        fieldnames = []

    with open(METRICS_CSV, "w", newline="", encoding="utf-8") as f_all:
        w_all = csv.DictWriter(f_all, fieldnames=fieldnames)
        w_all.writeheader()
        for r in all_rows:
            w_all.writerow(r)

    print("metrics_summary.csv -> " + METRICS_CSV)

    # save a copy of the code and a notebook version into the repo
    script_path_out = os.path.join(repo_dir, "cryo_sem_analysis.py")
    nb_path_out     = os.path.join(repo_dir, "cryo_sem_analysis.ipynb")

    with open(script_path_out, "w", encoding="utf-8") as f_py:
        f_py.write(code_text)

    write_notebook_copy(code_text, nb_path_out)

    print("saved: " + script_path_out)
    print("saved: " + nb_path_out)
    print("overlays and per-method metrics saved in: " + repo_dir)


# run it immediately when this script is executed
analyze_folder(BASE_DIR, GITHUB_DIR, SELF_SOURCE_CODE)
