from my_py_ai_utils import *

import numpy as np
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

    roc_res = make_roc_per_event([
        [ebsiso, ebbiso, 'tkIso'], 
        [ebspiso, ebbpiso, 'puppiIso'], 
        [ebspreiso, ebbpreiso, 'RePuppiIso']
        ], thrvs = np.arange(0, 10, 0.0125))

    make_roc_per_event_png(roc_res, filename="TkEleL2_EB_tkIso_ROCperevent_linear.png", xlim=(0.98, 1.001), ylim=(0.9, 1.01), s=5)

    roc_resee = make_roc_per_event([
        [eesiso, ebbiso, 'tkIso'], 
        [eespiso, eebpiso, 'puppiIso'], 
        [eespreiso, eebpreiso, 'RePuppiIso']
        ], thrvs = np.arange(0, 10, 0.0125))

    make_roc_per_event_png(roc_resee, filename="TkEleL2_EE_tkIso_ROCperevent_linear.png", xlim=(0.98, 1.001), ylim=(0.9, 1.01), s=5)


def make_roc_per_event_png(roc_res, *,
                            filename: str = "roc_curve.png",
                            scale: str = "default",
                            xlim: tuple[float] = (0.1, 1.1),
                            ylim: tuple[float] = (0.1, 1.1),
                            return_fig: bool = False,
                            **kwargs):


    markers = ['o', '*', 'v', '^']
    cmaps = ['viridis', 'plasma', 'inferno', 'coolwarm']

    # Plot ROC curves
    fig = plt.figure(figsize=(8, 6))
    for i, (fpr, tpr, thr, sample) in enumerate(roc_res):
        scatter = plt.scatter(fpr, tpr,
                                c=thr, cmap = cmaps[i], norm='log', # Colour maps can slow larger arrays
                                marker=markers[i],
                                label=sample, **kwargs)
        cbar = plt.colorbar(scatter)

    if scale == "piecewise_linear":
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
    plt.tight_layout()

    if return_fig:
        return fig
    else:
        print(f"Saving ROC curve to {filename}...")
        plt.savefig(filename, dpi=300)
        print(f"ROC curve saved to {filename}")
        plt.close()


if __name__ == '__main__':
    main()