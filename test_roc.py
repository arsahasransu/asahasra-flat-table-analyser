import awkward as ak
import numpy as np

si_pt = ak.Array([[10, 20], [5], [30]])
si = ak.Array([[0.1, 0.5], [0.9], [0.2]])

ptcut = 15
si_filtered = si[si_pt > ptcut]
print("si_filtered:", si_filtered)

si_nonempty = si_filtered[ak.num(si_filtered) > 0]
print("si_nonempty:", si_nonempty)

si_mins = ak.min(si_nonempty, axis=1)
print("si_mins:", si_mins)
print("si_mins numpy:", ak.to_numpy(si_mins))
