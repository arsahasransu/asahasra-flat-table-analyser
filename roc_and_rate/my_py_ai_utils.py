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

        true_val = np.concatenate([np.ones_like(si), np.zeros_like(bi)])
        pred_scores = np.concatenate([si, bi]) 

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
        si_nonempty = [ev for ev in si if len(ev) > 0]
        bi_nonempty = [ev for ev in bi if len(ev) > 0]

        nsi = len(si_nonempty)
        nbi = len(bi_nonempty)

        # Precompute per-event minimum — O(n_events)
        si_mins = np.array([ev.min() for ev in si_nonempty])
        bi_mins = np.array([ev.min() for ev in bi_nonempty])

        # Sort mins once — O(n_events log n_events)
        si_mins_sorted = np.sort(si_mins)
        bi_mins_sorted = np.sort(bi_mins)

        # Count events passing each threshold — O(n_thresholds log n_events)
        tpr_counts = np.searchsorted(si_mins_sorted, thrvs_sorted, side='right')
        fpr_counts = np.searchsorted(bi_mins_sorted, thrvs_sorted, side='right')

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

    def generate_mask(arr, ptlim):
        return np.vectorize(lambda x: x > ptlim, otypes=[object])(arr)

    def apply_mask(arr, mask):
        masked_arr = np.vectorize(lambda x, m: x[m], otypes=[object])(arr, mask)
        return np.array([x for x in masked_arr if x.size > 0], dtype=object)

    for ptcut in ptcuts:
        roc_res = []
        for (si, bi, sample) in vals:
            print(f"PT cut = {ptcut} GeV, sample = {sample}")

            # start = time.time()
            s_mask = generate_mask(si_pt, ptcut)
            si_filtered = apply_mask(si, s_mask)

            b_mask = generate_mask(bi_pt, ptcut)
            bi_filtered = apply_mask(bi, b_mask)
            # end = time.time()
            # print(f"Time taken to filter events: {end - start} seconds")

            nsi = si_filtered.shape[0]
            nbi = bi_filtered.shape[0]
            print(f"  Signal events: {nsi} / {si.shape[0]}")
            print(f"  Background events: {nbi} / {bi.shape[0]}")

            fprvs, tprvs = _compute_roc(si_filtered, bi_filtered)
            label = sample
            roc_res.append([fprvs, tprvs, thrvs_sorted, label])
        all_roc_res.append(roc_res)
    return all_roc_res