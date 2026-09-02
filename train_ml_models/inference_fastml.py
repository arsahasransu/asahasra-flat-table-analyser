"""
Inference for the fastml IsoSumNN soft-isolation models.

Analogous to the legacy ``train_ml_models/inference.py``, but driven by the
fastml package (``/mercury/data3/asahasra/RunLGNNDir/fastmlops/fastml``)
instead of the local ``SoftIsoSumNetwork``:

* The per-track-electron feature records are built with the same awkward
  pipeline as the training data (``ml_helper_functions``: cartesian
  tkel x puppi pairing with a dR cut, split by |pdgId|, keep the top-5
  pT candidates per species, pad with zeros, flatten).
* Feature selection (which columns the model consumes, and in what order)
  is delegated to the conventions of
  ``fastml.dataloaders.isosum_dataloader`` (``part_prefix`` /
  ``part_features``, ``bkg_prefix`` / ``bkg_features``, ``nbkg`` and its
  ``_match_cols`` / candidate-major layout checks), so the selection
  matches what ``IsoSumNN`` was trained on by construction.
* The model is ``fastml.models.isosum_nn.IsoSumNN``, checkpointed by
  ``fastml.training.train.train_model`` (``trained_models/best_model*.pt``)
  with the z-score statistics saved by
  ``fastml.dataloaders.isosum_dataloader.dataloader``
  (``norm_stats/*.pt``).

As in the legacy script, the per-track-electron weighted isolation sum is
written back into the snapshot ROOT file as the branch
``TkEleL2_{tkeltag}weighted_iso_score``.
"""

import os
import sys

import awkward as ak
import numpy as np
import torch
import uproot

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from ml_helper_functions import get_tkel_branches, get_puppi_branches
from ml_helper_functions import split_puppi_by_pdgid, flatten_puppi_collections
from ml_helper_functions import flatten_tkel_collection, pair_tkel_puppi_and_filter

# Root of the fastmlops repository containing the fastml package.
DEFAULT_FASTML_ROOT = "/mercury/data3/asahasra/RunLGNNDir/fastmlops"


def _load_fastml(fastml_root: str = DEFAULT_FASTML_ROOT):
    """Make the fastml package importable and return the inference parts.

    Returns:
        (isosum_dataloader, IsoSumNN, init_device)
    """
    root = os.path.abspath(fastml_root)
    if root not in sys.path:
        sys.path.insert(0, root)
    from fastml import isosum_dataloader
    from fastml.device import init_device
    from fastml.models.isosum_nn import IsoSumNN

    return isosum_dataloader, IsoSumNN, init_device


def _select_feature_indices(feature_names, dl):
    """Select model-input columns exactly as ``isosum_dataloader.dataloader`` does.

    Args:
        feature_names: Ordered list of column names (as in the training
            data 'features' list, or the awkward record fields here).
        dl: The ``fastml.dataloaders.isosum_dataloader`` module.

    Returns:
        (part_idx, bkg_idx): indices of the signal columns and of the
        background columns (candidate-major, ``len(bkg_features)`` per
        candidate) within ``feature_names``.

    Raises:
        ValueError: If the layout does not match the dataloader
            expectations (wrong column counts, pT not first in every
            candidate group).
    """
    part_idx = dl._match_cols(feature_names, dl.part_prefix, dl.part_features)
    bkg_idx = dl._match_cols(feature_names, dl.bkg_prefix, dl.bkg_features)

    if len(part_idx) != len(dl.part_features):
        raise ValueError(f"Expected {len(dl.part_features)} part columns, "
                         f"found {len(part_idx)}")
    expected_bkg = dl.nbkg * len(dl.bkg_features)
    if len(bkg_idx) != expected_bkg:
        raise ValueError(f"Expected {expected_bkg} bkg columns, "
                         f"found {len(bkg_idx)}")

    part_cols = [feature_names[i] for i in part_idx]
    if not part_cols[0].endswith("_pt"):
        raise ValueError(f"Expected pT as the first part feature, found {part_cols[0]}")
    n_bkg_feats = len(dl.bkg_features)
    bkg_cols = [feature_names[i] for i in bkg_idx]
    if not all(name.endswith("_pt") for name in bkg_cols[::n_bkg_feats]):
        raise ValueError("Expected pT as the first feature of each "
                         "background candidate group")

    return part_idx, bkg_idx


def add_soft_iso_to_root_file(file_path, tkeltag, puppitag, *,
                              tree_name='snapshot',
                              model_path, norm_stats_path,
                              fastml_root=DEFAULT_FASTML_ROOT,
                              dr_cut=0.5, max_items=5, batch_size=1024,
                              train_feature_names=None):
    """
    Reads a ROOT snapshot file with uproot, computes the Soft IsoSum with
    the fastml IsoSumNN model, and rewrites the file with a new
    'TkEleL2_{tkeltag}weighted_iso_score' branch added.

    Parameters:
      file_path: path to the *_snapshot.root file.
      tkeltag: tag for TkEle branches (e.g., 'Pt5_EB_MCH_').
      puppitag: tag for Puppi branches (e.g., 'Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_').
      tree_name: name of the tree in the ROOT file (default: 'snapshot').
      model_path: path to the IsoSumNN state dict saved by train_model
          (e.g., '<fastmlops>/trained_models/best_model_isosum_eb_nn.pt').
      norm_stats_path: path to the normalization stats saved by
          isosum_dataloader.dataloader / train_isosum
          (e.g., '<fastmlops>/norm_stats/norm_stats_eb_isosum_nn.pt').
      fastml_root: directory containing the fastml package.
      dr_cut: dR window for the tkel x puppi pairing (default: 0.5).
      max_items: max candidates kept per |pdgId| species (default: 5).
      batch_size: inference batch size.
      train_feature_names: optional 'features' list from a training .pt
          file; if given, the inference-time column selection is checked
          against it and a ValueError is raised on mismatch.
    """
    dl, IsoSumNN, init_device = _load_fastml(fastml_root)

    # 1. Load the snapshot tree and build the tkel / puppi records
    tkel_branches = get_tkel_branches(tkeltag)
    puppi_branches = get_puppi_branches(puppitag)

    with uproot.open(file_path) as f:
        data = f[tree_name].arrays()

    tkel = ak.zip({k: data[v] for k, v in tkel_branches.items()})
    puppi = ak.zip({k: data[v] for k, v in puppi_branches.items()})

    # 2. Same per-tkel pipeline as the training data
    tkel_with_puppi = pair_tkel_puppi_and_filter(tkel, puppi, dr_cut=dr_cut)
    tkel_with_puppi = split_puppi_by_pdgid(tkel_with_puppi)
    tkel_with_puppi = flatten_puppi_collections(tkel_with_puppi, max_items=max_items)
    flat_recs = flatten_tkel_collection(tkel_with_puppi)

    # Keep track of the event structure
    num_tkels = ak.num(flat_recs, axis=1)
    per_tkel_data = ak.flatten(flat_recs, axis=1)
    feature_names = list(per_tkel_data.fields)

    # 3. Select the model-input columns with the training dataloader's
    #    conventions (same order the model was trained on)
    part_idx, bkg_idx = _select_feature_indices(feature_names, dl)
    if train_feature_names is not None:
        t_part, t_bkg = _select_feature_indices(train_feature_names, dl)
        inferred = ([feature_names[i] for i in part_idx] +
                    [feature_names[i] for i in bkg_idx])
        trained = ([train_feature_names[i] for i in t_part] +
                   [train_feature_names[i] for i in t_bkg])
        if inferred != trained:
            raise ValueError("Inference feature selection does not match "
                             "the training data feature order")

    # 4. Stack features and split into signal / background
    features = [ak.to_numpy(per_tkel_data[f]) for f in feature_names]
    X = np.stack(features, axis=-1)
    x = torch.tensor(X, dtype=torch.float32)

    part = x[:, torch.tensor(part_idx)]
    bkg = x[:, torch.tensor(bkg_idx)].reshape(-1, dl.nbkg, len(dl.bkg_features))
    # Physical pT to sum over; cloned before the in-place log transform
    bkg_pt = bkg[:, :, 0].clone()

    # log(pT) transform, identical to the training-time dataloader
    part[:, 0] = torch.log(torch.clamp(part[:, 0], min=1e-6))
    bkg[:, :, 0] = torch.log(torch.clamp(bkg[:, :, 0], min=1e-6))

    # 5. z-score normalization with the saved training statistics
    norm_stats = torch.load(norm_stats_path, weights_only=True)
    tkel_mean, tkel_std = norm_stats['tkel_mean'], norm_stats['tkel_std']
    puppi_mean, puppi_std = norm_stats['puppi_mean'], norm_stats['puppi_std']
    if tuple(tkel_mean.shape) != (len(part_idx),):
        raise ValueError(f"tkel norm stats shape {tuple(tkel_mean.shape)} "
                         f"does not match {len(part_idx)} part features")
    if tuple(puppi_mean.shape) != (len(dl.bkg_features),):
        raise ValueError(f"puppi norm stats shape {tuple(puppi_mean.shape)} "
                         f"does not match {len(dl.bkg_features)} bkg features")

    part = (part - tkel_mean) / tkel_std
    bkg = (bkg - puppi_mean) / puppi_std

    # 6. Run inference with the fastml model
    device = init_device()
    model = IsoSumNN(part_dim=len(part_idx), bkg_dim=len(dl.bkg_features)).to(device)
    model.load_state_dict(torch.load(model_path, map_location=device, weights_only=True))
    model.eval()
    print(f"Loaded model from {model_path} on {device}")

    iso_sums = []
    dataset = torch.utils.data.TensorDataset(part, bkg, bkg_pt)
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)

    with torch.no_grad():
        for t_n, p_n, p_pt in loader:
            t_n, p_n, p_pt = t_n.to(device), p_n.to(device), p_pt.to(device)
            _, batch_iso_sum, _ = model(t_n, p_n, p_pt)
            iso_sums.append(batch_iso_sum.cpu().numpy())

    iso_sum_flat = (np.concatenate(iso_sums) if iso_sums
                    else np.array([], dtype=np.float32))

    # 7. Unflatten back to the event structure (N_events x N_tkels)
    iso_sum_unflat = ak.unflatten(iso_sum_flat, num_tkels)

    # 8. Add back to data and recreate the file
    new_branch_name = f"TkEleL2_{tkeltag}weighted_iso_score"
    data[new_branch_name] = iso_sum_unflat

    with uproot.recreate(file_path) as f:
        f[tree_name] = data

    print(f"Successfully added {new_branch_name} to {file_path}")


if __name__ == '__main__':
    # fastmlops_dir = DEFAULT_FASTML_ROOT
    fastmlops_dir = "/workspace/RunLGNNDir/fastmlops"
    model_dir = os.path.join(fastmlops_dir, 'trained_models')
    norm_dir = os.path.join(fastmlops_dir, 'norm_stats')
    data_dir = os.path.join(fastmlops_dir, 'DataForMaeen')

    snapshot_files = [
        # (snapshot file, tkel tag, puppi tag, region)
        ('../DoubleElectronGun_PU200_EB_snapshot.root', 'Pt5_EB_MCH_', 'Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_', 'eb'),
        ('../DY_PU200_EB_snapshot.root', 'Pt5_EB_MCH_', 'Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_', 'eb'),
        ('../MinBias_EB_snapshot.root', 'Pt5_EB_', 'Pt1_TkEleL2Pt5EB_0p0dR0p5_', 'eb'),
        ('../DoubleElectronGun_PU200_EE_snapshot.root', 'Pt5_EE_MCH_', 'Pt1_TkEleL2Pt5EEMCH_0p0dR0p5_', 'ee'),
        ('../DY_PU200_EE_snapshot.root', 'Pt5_EE_MCH_', 'Pt1_TkEleL2Pt5EEMCH_0p0dR0p5_', 'ee'),
        ('../MinBias_EE_snapshot.root', 'Pt5_EE_', 'Pt1_TkEleL2Pt5EE_0p0dR0p5_', 'ee'),
    ]

    # Load each region's training feature list once, so the inference-time
    # column selection can be checked against what the model actually saw.
    train_features = {}
    for region in ('eb', 'ee'):
        train_data_path = os.path.join(data_dir, f'{region}_train_data.pt')
        if os.path.exists(train_data_path):
            train_features[region] = torch.load(train_data_path,
                                                weights_only=True)['features']
        else:
            print(f"Training data {train_data_path} not found; "
                  f"skipping feature-order check for {region}.")

    for fname, tkeltag, puppitag, region in snapshot_files:
        if not os.path.exists(fname):
            print(f"File {fname} not found, skipping.")
            continue
        model_path = os.path.join(model_dir, f'best_model_isosum_{region}_nn.pt')
        norm_path = os.path.join(norm_dir, f'norm_stats_{region}_isosum_nn.pt')
        print(f"Processing {fname} ...")
        add_soft_iso_to_root_file(fname, tkeltag, puppitag,
                                  model_path=model_path,
                                  norm_stats_path=norm_path,
                                  train_feature_names=train_features.get(region),
                                  fastml_root=fastmlops_dir)
