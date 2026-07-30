import multiprocessing
import time

import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

# import my_py_generic_utils as ut


def roc_curve(truth: np.ndarray, scores: np.ndarray, *,
              points: int = 1000):
    # Sort scores and corresponding truth values
    # Sort ascending where signal favours lower values
    asc_score_indices = np.argsort(scores)
    sorted_scores = scores[asc_score_indices]
    sorted_truth = truth[asc_score_indices]

    # Calculate TPR and FPR and AUC
    tpr = np.cumsum(sorted_truth) / np.sum(sorted_truth)
    fpr = np.cumsum(1 - sorted_truth) / np.sum(1 - sorted_truth)

    # Remove non-unique FPR
    unique_fpr, inv = np.unique(fpr, return_inverse=True)
    indices = np.zeros(len(unique_fpr), dtype=int)
    indices[inv] = np.arange(len(fpr))
    tpr = tpr[indices]
    fpr = fpr[indices]
    sorted_scores = sorted_scores[indices]

    # Choose the ROC region of interest
    indices = tpr > 0.5
    tpr = tpr[indices]
    fpr = fpr[indices]
    sorted_scores = sorted_scores[indices]

    # Sample randomly proportional to inverse frequency of weights
    res = 0.01
    binned_tpr = np.round(tpr / res) * res
    unique_tpr_bins, counts_tpr = np.unique(binned_tpr, return_counts=True)
    tpr_count_map = dict(zip(unique_tpr_bins, counts_tpr))
    weights_tpr = np.array([1 / tpr_count_map[x] for x in binned_tpr])
    weights_tpr /= weights_tpr.sum()
    n_samples = points
    indices = np.random.choice(len(tpr), size=n_samples, replace=True, p=weights_tpr)

    tpr = tpr[indices]
    fpr = fpr[indices]
    auc = np.trapezoid(tpr, fpr)
    sorted_scores = sorted_scores[indices]

    return fpr, tpr, auc, sorted_scores


def make_roc_png(roc_res, *,
                 filename: str = "roc_curve.png",
                 scale: str = "default",
                 xlim: tuple[float] = (0.1, 1.1), ylim: tuple[float] = (0.1, 1.1), **kwargs):


    markers = ['o', '*', 'v', '^']
    cmaps = ['viridis', 'plasma', 'inferno', 'coolwarm']

    # Plot ROC curves
    print(f"Saving ROC curve to {filename}...")
    plt.figure(figsize=(8, 6))
    for i, (fpr, tpr, thr, sample) in enumerate(roc_res):
        print(fpr[1])
        scatter = plt.scatter(fpr, tpr,
                                c=thr, cmap = cmaps[i], norm='log', # Colour maps can slow larger arrays
                                marker=markers[i],
                            #   label=f'{vals[i][2]} (AUC = {auc:.4f})', **kwargs)
                                label=sample, **kwargs)
        plt.colorbar(scatter)

    if scale == "piecewise_linear":
        # Adjust plotting to show 0.1 to 0.9 and 0.9 to 0.99 and 0.99 to 1.0 and 1.0 to 1.1 regions clearly
        # Optional: vertical guides at the two boundaries for clarity
        plt.xscale('piecewise_0p1_0p9_0p999')
        plt.yscale('piecewise_0p1_0p9_0p999')
        plt.xlim(xlim[0], xlim[1])
        plt.ylim(ylim[0], ylim[1])
        for v in [0.9, 0.99, 1.0]:
            plt.axvline(v, color='0.7', lw=1, ls='--')
    else:
        plt.xlim(left=xlim[0], right=xlim[1])
        plt.ylim(bottom=ylim[0], top=ylim[1])

    plt.grid(True, which='both', ls='--', alpha=0.4)

    plt.xlabel('False Positive Rate (log scale)')
    plt.ylabel('True Positive Rate (log scale)')
    plt.title('ROC Curve')

    plt.legend(loc='lower right')
    plt.savefig(filename, dpi=300)
    print(f"ROC curve saved to {filename}")
    plt.close()
    

# @ut.time_eval
def make_roc(vals: list[list[np.ndarray]], *,
                         filename: str = "roc_curve.png",
                         scale: str = "default",
                         points: int = 1000, 
                         xlim: tuple[float] = (0.1, 1.1), ylim: tuple[float] = (0.1, 1.1), **kwargs):

    roc_res = []

    for (si, bi, sample) in vals:
        si_flat = ak.to_numpy(ak.flatten(si, axis=None))
        bi_flat = ak.to_numpy(ak.flatten(bi, axis=None))

        true_val = np.concatenate([np.ones_like(si_flat), np.zeros_like(bi_flat)])
        pred_scores = np.concatenate([si_flat, bi_flat]) 

        fpr, tpr, auc, thr = roc_curve(true_val, pred_scores, points = points)
        print(f"AUC: {auc:.4f}")
        roc_res.append([fpr, tpr, thr, sample])

    make_roc_png(roc_res, filename=filename, scale=scale,
                 xlim=xlim, ylim=ylim, **kwargs)


def make_roc_per_event(vals, *,
                       thrvs: np.ndarray = np.arange(0, 10.1, 0.1),
                       si_pt: np.ndarray | None = None,
                       bi_pt: np.ndarray | None = None,
                       ptcuts: list[float] | None = None):

    thrvs_sorted = np.sort(thrvs)  # searchsorted requires sorted array

    def _compute_roc(si, bi):
        # Filter to non-empty events only
        si_nonempty = si[ak.num(si) > 0]
        bi_nonempty = bi[ak.num(bi) > 0]

        nsi = len(si_nonempty)
        nbi = len(bi_nonempty)

        # Precompute per-event minimum — O(n_events)
        si_mins_np = np.sort(ak.to_numpy(ak.min(si_nonempty, axis=1)))
        bi_mins_np = np.sort(ak.to_numpy(ak.min(bi_nonempty, axis=1)))

        # Count events passing each threshold — O(n_thresholds log n_events)
        tpr_counts = np.searchsorted(si_mins_np, thrvs_sorted, side='right')
        fpr_counts = np.searchsorted(bi_mins_np, thrvs_sorted, side='right')

        tprvs = tpr_counts / nsi
        fprvs = fpr_counts / nbi

        return fprvs.tolist(), tprvs.tolist()

    if ptcuts is None or si_pt is None or bi_pt is None:
        # Original behaviour — no PT filtering
        roc_res = []
        for (si, bi, sample) in vals:
            print(sample)
            nsi = len(si)
            nbi = len(bi)
            print(f"Total signal events: {nsi}, background: {nbi}")
            fprvs, tprvs = _compute_roc(si, bi)
            roc_res.append([fprvs, tprvs, thrvs_sorted, sample])
        return roc_res

    # PT filtering active
    ptcuts = sorted(ptcuts)
    all_roc_res = []

    for ptcut in ptcuts:
        roc_res = []
        for (si, bi, sample) in vals:
            print(f"PT cut = {ptcut} GeV, sample = {sample}")

            si_filtered = si[si_pt > ptcut]
            bi_filtered = bi[bi_pt > ptcut]

            nsi = ak.sum(ak.num(si_filtered) > 0)
            nbi = ak.sum(ak.num(bi_filtered) > 0)
            print(f"  Signal events: {nsi} / {len(si)}")
            print(f"  Background events: {nbi} / {len(bi)}")

            fprvs, tprvs = _compute_roc(si_filtered, bi_filtered)
            label = sample
            roc_res.append([fprvs, tprvs, thrvs_sorted, label])
        all_roc_res.append(roc_res)
    return all_roc_res