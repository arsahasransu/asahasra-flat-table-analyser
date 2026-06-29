import multiprocessing
import time

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

        true_val = np.concatenate([np.ones_like(si), np.zeros_like(bi)])
        pred_scores = np.concatenate([si, bi]) 

        fpr, tpr, auc, thr = roc_curve(true_val, pred_scores, points = points)
        print(f"AUC: {auc:.4f}")
        roc_res.append([fpr, tpr, thr, sample])

    make_roc_png(roc_res, filename=filename, scale=scale,
                 xlim=xlim, ylim=ylim, **kwargs)


def make_roc_per_event(vals, *,
                       thrvs: np.ndarray = np.arange(0, 10.1, 0.1)):

    roc_res = []
    thrvs_sorted = np.sort(thrvs)  # searchsorted requires sorted array

    for (si, bi, sample) in vals:
        print(sample)

        # Filter to non-empty events only
        si_nonempty = [ev for ev in si if len(ev) > 0]
        bi_nonempty = [ev for ev in bi if len(ev) > 0]

        nsi = len(si_nonempty)
        nbi = len(bi_nonempty)

        print("Total event count: ",si.shape[0], ". Non-empty: ", nsi)
        print("Total event count: ",bi.shape[0], ". Non-empty: ", nbi)

        # Precompute per-event minimum — O(n_events)
        si_mins = np.array([ev.min() for ev in si_nonempty])
        bi_mins = np.array([ev.min() for ev in bi_nonempty])

        # Sort mins once — O(n_events log n_events)
        si_mins_sorted = np.sort(si_mins)
        bi_mins_sorted = np.sort(bi_mins)

        # Count events passing each threshold — O(n_thresholds log n_events)
        # searchsorted gives index of first element >= thrv,
        # so events below thrv = that index count
        tpr_counts = np.searchsorted(si_mins_sorted, thrvs_sorted, side='right')
        fpr_counts = np.searchsorted(bi_mins_sorted, thrvs_sorted, side='right')

        tprvs = tpr_counts / nsi
        fprvs = fpr_counts / nbi

        roc_res.append([fprvs.tolist(), tprvs.tolist(), thrvs_sorted, sample])

    return roc_res