import itertools
import multiprocessing
import time

import numpy as np
import matplotlib.pyplot as plt

import my_py_generic_utils as ut


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


@ut.time_eval
def make_roc(vals: list[list[np.ndarray]], *,
                         filename: str = "roc_curve.png",
                         scale: str = "default",
                         points: int = 1000, 
                         xlim: tuple[float] = (0.1, 1.1), ylim: tuple[float] = (0.1, 1.1), **kwargs):

    roc_res = []
    markers = ['o', '*', 'v', '^']
    cmaps = ['viridis', 'plasma', 'inferno', 'coolwarm']

    for i, (si, bi, _) in enumerate(vals):

        true_val = np.concatenate([np.ones_like(si), np.zeros_like(bi)])
        pred_scores = np.concatenate([si, bi]) 

        fpr, tpr, auc, thr = roc_curve(true_val, pred_scores, points = points)
        print(f"AUC: {auc:.4f}")
        roc_res.append([fpr, tpr, auc, thr])

    # Plot ROC curves
    print(f"Saving ROC curve to {filename}...")
    plt.figure(figsize=(8, 6))
    for i, (fpr, tpr, auc, thr) in enumerate(roc_res):
        scatter = plt.scatter(fpr, tpr,
                              c=thr, cmap = cmaps[i], norm='log', # Colour maps can slow larger arrays
                              marker=markers[i],
                            #   label=f'{vals[i][2]} (AUC = {auc:.4f})', **kwargs)
                              label=f'{vals[i][2]}', **kwargs)
        cbar = plt.colorbar(scatter)
    
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


def _init_worker(valarr):
    global _VALARRAY
    _VALARRAY = valarr


@ut.time_eval
def count_events_passing_threshold(valarr, thrv):
    return sum(any(vals < thrv for vals in ev_valarr) for ev_valarr in valarr)


def count_events_passing_threshold_multiprocessing(thrv):
    return sum(any(vals < thrv for vals in ev_valarr) for ev_valarr in _VALARRAY)


@ut.time_eval
def make_roc_per_event(vals: list[list[np.ndarray]], *,
                       points: int = 1000,
                       thrvs: np.array = np.arange(0, 10, 0.1)):
    
    roc_res = []
    for i, (si, bi, sample) in enumerate(vals):
        print(sample)

        starttime = time.time()
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()-2,
                                  initializer=_init_worker,
                                  initargs=(si,)) as mp_pool:
            tprvs = mp_pool.map(count_events_passing_threshold_multiprocessing, thrvs)
        endtime = time.time()
        print(f"Execution time: {endtime - starttime} seconds")

        nsi = si.shape[0]
        tprvs = [tprv/nsi for tprv in tprvs]

        starttime = time.time()
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()-2,
                                  initializer=_init_worker,
                                  initargs=(bi,)) as mp_pool:
            fprvs = mp_pool.map(count_events_passing_threshold_multiprocessing, thrvs)
        endtime = time.time()
        print(f"Execution time: {endtime - starttime} seconds")

        nbi = bi.shape[0]
        fprvs = [fprv/nbi for fprv in fprvs]

        roc_res.append([fprvs, tprvs, thrvs, sample])

    return roc_res