from my_py_ai_utils import *

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

import sys
import os
from os.path import join

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import rdf_generic as rdf_g


def main():
    base = join(os.path.dirname(__file__), '..')
    fsig = rdf_g.load_rdf_snapshot_from_root(join(base, 'DY_PU200_EB_snapshot.root'))
    fbkg = rdf_g.load_rdf_snapshot_from_root(join(base, 'MinBias_EB_snapshot.root'))
    fsigee = rdf_g.load_rdf_snapshot_from_root(join(base, 'DY_PU200_EE_snapshot.root'))
    fbkgee = rdf_g.load_rdf_snapshot_from_root(join(base, 'MinBias_EE_snapshot.root'))

    # Each row is an event entry in the loaded data below

    ebspt = fsig["TkEleL2_EB_MCH_pt"]
    ebsiso = fsig["TkEleL2_EB_MCH_tkIso"]
    ebspiso = fsig["TkEleL2_EB_MCH_puppiIso"]
    ebspreiso = fsig["TkEleL2_EB_MCH_reisotot_dRmin0_03_puppiIso"]
    ebbpt = fbkg["TkEleL2_Pt5_EB_pt"]
    ebbiso = fbkg["TkEleL2_Pt5_EB_tkIso"]
    ebbpiso = fbkg["TkEleL2_Pt5_EB_puppiIso"]
    ebbpreiso = fbkg["TkEleL2_Pt5_EB_reisotot_dRmin0_03_puppiIso"]
    print(ebspt.shape, ebsiso.shape, ebspiso.shape, ebspreiso.shape)
    print(ebbpt.shape, ebbiso.shape, ebbpiso.shape, ebbpreiso.shape)


    eespt = fsigee["TkEleL2_EE_MCH_pt"]
    eesiso = fsigee["TkEleL2_EE_MCH_tkIso"]
    eespiso = fsigee["TkEleL2_EE_MCH_puppiIso"]
    eespreiso = fsigee["TkEleL2_EE_MCH_reisotot_dRmin0_03_puppiIso"]
    eebpt = fbkgee["TkEleL2_Pt5_EE_pt"]
    eebiso = fbkgee["TkEleL2_Pt5_EE_tkIso"]
    eebpiso = fbkgee["TkEleL2_Pt5_EE_puppiIso"]
    eebpreiso = fbkgee["TkEleL2_Pt5_EE_reisotot_dRmin0_03_puppiIso"]
    print(eespt.shape, eesiso.shape, eespiso.shape, eespreiso.shape)
    print(eebpt.shape, eebiso.shape, eebpiso.shape, eebpreiso.shape)

    # With the data loaded as above, you will generate performance of iso per event
    # In order to generate per electron, change to the following code below
    # eespt = np.concatenate(ebspt)

    roc_res_eb = make_roc_per_event([
        [ebsiso, ebbiso, 'tkIso'], 
        [ebspiso, ebbpiso, 'puppiIso'], 
        [ebspreiso, ebbpreiso, 'RePuppiIso']
        ], thrvs = np.arange(0, 10.0125, 0.0125))

    make_roc_per_event_png(roc_res_eb, 
                           filename="TkEleL2_EB_tkIso_ROCperevent_piecewise_linear.png", 
                           scale = "piecewise_linear",
                           xlim=(0.98, 1.001), 
                           ylim=(0.1, 1.01), 
                           s=5)

    make_roc_per_event_png(roc_res_eb, 
                           filename="TkEleL2_EB_tkIso_ROCperevent_linear.png", 
                           xlim=(0.992, 1.001), 
                           ylim=(0.8, 1.01), 
                           s=5)

    roc_res_ee = make_roc_per_event([
        [eesiso, ebbiso, 'tkIso'], 
        [eespiso, eebpiso, 'puppiIso'], 
        [eespreiso, eebpreiso, 'RePuppiIso']
        ], thrvs = np.arange(0, 10, 0.0125))

    make_roc_per_event_png(roc_res_ee, 
                           filename="TkEleL2_EE_tkIso_ROCperevent_piecewise_linear.png", 
                           scale = "piecewise_linear",
                           xlim=(0.1, 1.001), 
                           ylim=(0.1, 1.01), 
                           s=5)

    make_roc_per_event_png(roc_res_ee, 
                           filename="TkEleL2_EE_tkIso_ROCperevent_linear.png", 
                           xlim=(0.992, 1.001), 
                           ylim=(0.8, 1.01), 
                           s=5)


# ── Global style ─────────────────────────────────────────────────────────────
plt.rcParams.update({
    # Font
    "font.family":        "serif",
    "font.serif":         ["Arial"],
    "font.size":          13,
    "axes.titlesize":     15,
    "axes.titleweight":   "bold",
    "axes.labelsize":     13,
    "xtick.labelsize":    11,
    "ytick.labelsize":    11,
    "legend.fontsize":    11,

    # Axes
    "axes.linewidth":     1.0,
    "axes.spines.top":    True,
    "axes.spines.right":  True,

    # Ticks – inward, on all four sides
    "xtick.direction":    "in",
    "ytick.direction":    "in",
    "xtick.top":          True,
    "ytick.right":        True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.width":  1.0,
    "ytick.major.width":  1.0,
    "xtick.minor.width":  0.75,
    "ytick.minor.width":  0.75,
    "xtick.major.size":   5,
    "ytick.major.size":   5,
    "xtick.minor.size":   2.5,
    "ytick.minor.size":   2.5,

    # Figure
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
})

# Perceptually uniform, print-safe, colour-blind-friendly colour maps
# (each is distinguishable in greyscale too)
_CMAPS   = ["viridis", "cividis", "plasma", "inferno"]
_MARKERS = ["o", "s", "^", "D"]
_MSIZES  = [18, 20, 22, 20]   # pt²; marker area passed to scatter


def make_roc_per_event_png(
    roc_res,
    *,
    filename: str = "roc_curve.png",
    scale: str = "default",
    xlim: tuple = (0.9980, 1.001),
    ylim: tuple = (0.90,   1.01),
    auc_labels: list | None = None,     # optional list of AUC strings
    xlabel: str = "Background Acceptance",
    ylabel: str = "Signal Efficiency",
    return_fig: bool = False,
    **kwargs,
):
    n = len(roc_res)

    # ── Figure & gridspec ────────────────────────────────────────────────────
    # Main axes takes 70 % of width; one narrow colorbar per curve in the rest
    fig = plt.figure(figsize=(7 + 0.65 * n, 7))
    gs  = fig.add_gridspec(
        1, 2,
        left=0.10, right=0.97,
        wspace=0.06,
        width_ratios=[7, 0.65 * n],
    )
    ax = fig.add_subplot(gs[0, 0])
        
    cbar_gs   = gs[0, 1].subgridspec(1, n, wspace=1.05)
    cbar_axes = [fig.add_subplot(cbar_gs[i]) for i in range(n)]

    # ── Scatter + colorbars ──────────────────────────────────────────────────
    leg_labels = []
    for i, (fpr, tpr, thr, label) in enumerate(roc_res):
        thr_arr = np.asarray(thr)
        thr_pos = np.where(thr_arr > 0, thr_arr, np.nan)   # hide ≤0 from LogNorm

        # Build per-curve norm so each colorbar spans its own range
        vmin, vmax = np.nanmin(thr_pos), np.nanmax(thr_pos)
        norm = LogNorm(vmin=vmin, vmax=vmax)

        leg_label = label if auc_labels is None else f"{label}  (AUC = {auc_labels[i]})"
        leg_labels.append(leg_label)

        # 's' (marker size) may arrive via **kwargs; honour the caller's value
        # if supplied, otherwise use our per-curve default.
        marker_size = kwargs.pop("s", _MSIZES[i % len(_MSIZES)])

        sc = ax.scatter(
            fpr, tpr,
            c=thr_pos,
            cmap=_CMAPS[i % len(_CMAPS)],
            norm=norm,
            marker=_MARKERS[i % len(_MARKERS)],
            s=marker_size,
            linewidths=0.0,         # no edge → cleaner at small sizes
            alpha=0.85,
            label=leg_label,
            zorder=3,
            **kwargs,
        )

        cb = fig.colorbar(sc, cax=cbar_axes[i], orientation="vertical")
        cb.ax.tick_params(labelsize=11, width=0.9, length=5, direction="in")

        # Label only the rightmost bar to avoid clutter
        if i == n - 1:
            cb.set_label("Threshold", fontsize=11, labelpad=5)
        else:
            cb.set_label(label, fontsize=11, labelpad=5, rotation=90)

        # Minor ticks on the colorbar
        cb.ax.yaxis.set_minor_locator(ticker.LogLocator(subs="all", numticks=10))
        cb.ax.yaxis.set_minor_formatter(ticker.NullFormatter())
        cb.ax.yaxis.set_major_locator(ticker.LogLocator())
        cb.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:g}'))

    # ── Axes scale & limits ──────────────────────────────────────────────────
    if scale == "piecewise_linear":
        ax.set_xscale("piecewise_0p1_0p9_0p999")
        ax.set_yscale("piecewise_0p1_0p9_0p999")
        for v in [0.9, 0.99, 1.0]:
            ax.axvline(v, color="0.75", lw=0.8, ls="--", zorder=1)
    else:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))  # tick every 0.1
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.005))  # tick every 0.1
    ax.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.3f'))

    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))  # tick every 0.1
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

    # ── Grid ─────────────────────────────────────────────────────────────────
    ax.grid(True, which="major", ls="--", lw=0.6, color="0.75", alpha=0.7, zorder=0)
    ax.grid(True, which="minor", ls=":",  lw=0.4, color="0.82", alpha=0.5, zorder=0)

    # ── Labels & title ───────────────────────────────────────────────────────
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # ── Legend ───────────────────────────────────────────────────────────────
    # Build custom Line2D handles so markers are guaranteed visible in the legend
    handles = []
    for i, (_, _, _, _) in enumerate(roc_res):
        handles.append(
            Line2D(
                [0], [0],
                marker=_MARKERS[i % len(_MARKERS)],
                linestyle='',
                color=mpl.colormaps[_CMAPS[i % len(_CMAPS)]](0.5),
                markersize=8,
                markeredgecolor='k',
                markeredgewidth=0.5,
                label=leg_labels[i],
            )
        )
    ax.legend(
        handles=handles,
        loc="lower right",
        frameon=False,
        borderpad=0.6,
        labelspacing=0.4,
        handletextpad=0.5,
    )

    # ── Output ───────────────────────────────────────────────────────────────
    if return_fig:
        return fig

    print(f"Saving ROC curve → {filename}")
    fig.savefig(filename)
    plt.close(fig)
    print("Done.")


if __name__ == '__main__':
    main()