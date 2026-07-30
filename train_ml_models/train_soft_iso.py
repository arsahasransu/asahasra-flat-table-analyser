"""
NN training script with the following constraints:
 1. The Shared MLP (Selector): A neural network 
    (Linear(num_features, _) -> ReLU -> ... -> Linear(_, 1) -> Sigmoid) 
    analyzes each of the PUPPI candidates by concatenating its features with the tkel 
    features. It outputs a weight w_i between 0 and 1 for each candidate.
 2. The Soft IsoSum: The network computes iso_sum = sum(w_i * puppi_pt_i) 
    over all puppi candidates, calculating the unnormalized, physical sum in GeV.
 3. The Option A Constraint: To map this physical isolation sum to a 
    binary classification probability, the network uses a learnable threshold and scale:
    Logits = (Threshold - IsoSum) * Scale
    Because IsoSum is subtracted, the network is strictly mathematically 
    forced to minimize IsoSum for Signal (to yield positive logits/probability > 0.5) 
    and maximize IsoSum for Background (to yield negative logits/probability < 0.5).
"""

import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader
import time
import math
import matplotlib.pyplot as plt
import datetime

import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from ml_helper_functions import SoftIsoSumNetwork, tkel_bs, puppi_bs, puppicands

def prepare_dataloaders(train_path, test_path, batch_size=1024):
    print("Loading data...")
    train_data = torch.load(train_path, weights_only=True)
    test_data = torch.load(test_path, weights_only=True)

    x_train, y_train = train_data['x'], train_data['y'].float()
    x_test, y_test = test_data['x'], test_data['y'].float()

    # Count the number of tkel features in data
    data_features = train_data['features']
    tkel_index = sum([feature.startswith('tkel') for feature in data_features])
    if tkel_index != tkel_bs:
        raise(f'Required exact match of tkel features: {tkel_bs}, found {tkel_index}')
    puppi_features = sum([feature.startswith('mpuppi') for feature in data_features])
    if puppi_features // puppicands != puppi_bs:
        raise(f'Required exact match of puppi features: {puppi_bs}, '\
              'found {puppi_features}')
    
    # Split into tkel and puppi
    tkel_train = x_train[:, 0:tkel_bs]
    puppi_train_flat = x_train[:, tkel_bs:]
    if puppi_features%puppicands == 0:
        puppi_train = puppi_train_flat.reshape(-1, puppicands, puppi_bs)
    else:
        raise UnboundLocalError(f"Unresolved splitting of puppi cands, found "\
                f"features {puppi_features} for {puppicands} puppi canddiates")
    
    tkel_test = x_test[:, 0:tkel_bs]
    puppi_test_flat = x_test[:, tkel_bs:]
    if puppi_features%puppicands == 0:
        puppi_test = puppi_test_flat.reshape(-1, puppicands, puppi_bs)
    else:
        raise UnboundLocalError(f"Unresolved splitting of puppi cands, found "\
                f"features {puppi_features} for {puppicands} puppi canddiates")
    
    # Extract unnormalized pT (which is the 0-th feature of the 14 puppi features)
    puppi_pt_train = puppi_train[:, :, 0].clone()
    puppi_pt_test = puppi_test[:, :, 0].clone()
    
    # Apply log(pT) transformation in-place on the feature arrays before calculating stats
    tkel_train[:, 0] = torch.log(torch.clamp(tkel_train[:, 0], min=1e-6))
    tkel_test[:, 0] = torch.log(torch.clamp(tkel_test[:, 0], min=1e-6))
    
    puppi_train[:, :, 0] = torch.log(torch.clamp(puppi_train[:, :, 0], min=1e-6))
    puppi_test[:, :, 0] = torch.log(torch.clamp(puppi_test[:, :, 0], min=1e-6))
    
    # Compute Normalization statistics on Train only
    tkel_mean = tkel_train.mean(dim=0)
    tkel_std = tkel_train.std(dim=0)
    tkel_std[tkel_std < 1e-6] = 1.0 # Prevent division by zero
    
    # Normalize all 25 candidates uniformly
    puppi_train_reshaped = puppi_train.reshape(-1, puppi_bs)
    puppi_mean = puppi_train_reshaped.mean(dim=0)
    puppi_std = puppi_train_reshaped.std(dim=0)
    puppi_std[puppi_std < 1e-6] = 1.0
    
    # Apply normalization
    tkel_train_norm = (tkel_train - tkel_mean) / tkel_std
    tkel_test_norm = (tkel_test - tkel_mean) / tkel_std
    
    puppi_train_norm = (puppi_train - puppi_mean) / puppi_std
    puppi_test_norm = (puppi_test - puppi_mean) / puppi_std
    
    train_dataset = TensorDataset(tkel_train_norm, puppi_train_norm, puppi_pt_train, y_train)
    test_dataset = TensorDataset(tkel_test_norm, puppi_test_norm, puppi_pt_test, y_test)
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, drop_last=True)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    
    # Save normalization stats for later inference
    norm_stats = {
        'tkel_mean': tkel_mean, 'tkel_std': tkel_std,
        'puppi_mean': puppi_mean, 'puppi_std': puppi_std
    }
    
    return train_loader, test_loader, norm_stats


def train_epoch(model, loader, optimizer, criterion, device):
    model.train()
    total_loss = 0.0
    correct = 0
    total = 0
    
    for tkel_norm, puppi_norm, puppi_pt, y in loader:
        tkel_norm = tkel_norm.to(device)
        puppi_norm = puppi_norm.to(device)
        puppi_pt = puppi_pt.to(device)
        y = y.to(device)
        
        optimizer.zero_grad()
        logits, iso_sum, weights = model(tkel_norm, puppi_norm, puppi_pt)
        
        loss = criterion(logits, y)
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item() * y.size(0)
        
        # Accuracy
        preds = (logits > 0).float()
        correct += (preds == y).sum().item()
        total += y.size(0)
        
    return total_loss / total, correct / total


@torch.no_grad()
def eval_epoch(model, loader, criterion, device):
    model.eval()
    total_loss = 0.0
    correct = 0
    total = 0
    
    iso_sig_sum = 0.0
    iso_bkg_sum = 0.0
    sig_count = 0
    bkg_count = 0
    
    for tkel_norm, puppi_norm, puppi_pt, y in loader:
        tkel_norm = tkel_norm.to(device)
        puppi_norm = puppi_norm.to(device)
        puppi_pt = puppi_pt.to(device)
        y = y.to(device)
        
        logits, iso_sum, weights = model(tkel_norm, puppi_norm, puppi_pt)
        
        loss = criterion(logits, y)
        total_loss += loss.item() * y.size(0)
        
        preds = (logits > 0).float()
        correct += (preds == y).sum().item()
        total += y.size(0)
        
        # Track avg iso_sum for monitoring
        sig_mask = (y == 1)
        bkg_mask = (y == 0)
        iso_sig_sum += iso_sum[sig_mask].sum().item()
        sig_count += sig_mask.sum().item()
        iso_bkg_sum += iso_sum[bkg_mask].sum().item()
        bkg_count += bkg_mask.sum().item()
        
    avg_iso_sig = iso_sig_sum / sig_count if sig_count > 0 else 0
    avg_iso_bkg = iso_bkg_sum / bkg_count if bkg_count > 0 else 0
        
    return total_loss / total, correct / total, avg_iso_sig, avg_iso_bkg


def plot_epoch(history, epochs, *,
               savefile_ext = ""):
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    plt.figure(figsize=(10, 5))
    plt.plot(range(1, epochs + 1), history['train_loss'], label='Train Loss')
    plt.plot(range(1, epochs + 1), history['val_loss'], label='Test Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title('Training and Testing Loss')
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'loss_plot{savefile_ext}_{timestamp}.png')
    plt.close()

    train_oneminusacc = [1-acc for acc in history['train_acc']]
    val_oneminusacc = [1-acc for acc in history['val_acc']]
    plt.figure(figsize=(10, 5))
    plt.plot(range(1, epochs + 1), train_oneminusacc, label='Train Accuracy')
    plt.plot(range(1, epochs + 1), val_oneminusacc, label='Test Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('1 - Accuracy')
    plt.title('Training and Testing Accuracy')
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'accuracy_plot{savefile_ext}_{timestamp}.png')
    plt.close()
    print(f"Saved plots with timestamp: {timestamp}")


def train_model(device, train_path, test_path, *,
                savefile_ext = "",
                epochs = 60,
                batch_size = 1024):
    # Increased batch size for better multi-core throughput
    train_loader, test_loader, norm_stats = prepare_dataloaders(train_path, test_path, batch_size=batch_size)
    torch.save(norm_stats, f'norm_stats{savefile_ext}.pt')
    
    model = SoftIsoSumNetwork(tkel_dim=tkel_bs, puppi_dim=puppi_bs, hidden_dims=[32, 32]).to(device)
    
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    criterion = nn.BCEWithLogitsLoss()
    
    epochs = epochs
    best_acc = 0.0
    
    history = {
        'train_loss': [],
        'train_acc': [],
        'val_loss': [],
        'val_acc': []
    }
    
    print("Starting training...")
    for epoch in range(epochs):
        t0 = time.time()
        
        train_loss, train_acc = train_epoch(model, train_loader, optimizer, criterion, device)
        val_loss, val_acc, avg_iso_sig, avg_iso_bkg = eval_epoch(model, test_loader, criterion, device)
        
        history['train_loss'].append(train_loss)
        history['train_acc'].append(train_acc)
        history['val_loss'].append(val_loss)
        history['val_acc'].append(val_acc)
        
        t1 = time.time()
        print(f"Epoch {epoch+1:02d}/{epochs} | Time: {t1-t0:.1f}s | "
              f"Train Loss: {train_loss:.4f} Acc: {train_acc:.4f} | "
              f"Val Loss: {val_loss:.4f} Acc: {val_acc:.4f} | "
              f"Iso (Sig/Bkg): {avg_iso_sig:.2f} / {avg_iso_bkg:.2f}")
              
        if val_acc > best_acc:
            best_acc = val_acc
            torch.save(model.state_dict(), f'best_soft_iso_model{savefile_ext}.pt')
            
    print(f"Training complete. Best Validation Accuracy: {best_acc:.4f}")
    
    # Print final learned parameters
    model.load_state_dict(torch.load(f'best_soft_iso_model{savefile_ext}.pt', weights_only=True))
    print(f"Learned Threshold (Iso Cut): {model.threshold.item():.3f}")
    print(f"Learned Scale: {F.softplus(model.scale_raw).item():.3f}")

    plot_epoch(history, epochs, savefile_ext = savefile_ext)


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    if device.type == 'cpu':
        # Optimize CPU threading for multi-core scaling
        cores = os.cpu_count()
        # Fallback to os.cpu_count(), but respect batch scheduler core limits if present
        num_threads = int(os.environ.get('SLURM_CPUS_PER_TASK', cores if cores else 4))
        torch.set_num_threads(num_threads)
        print(f"Using device: {device} with {num_threads} intra-op threads")
    else:
        print(f"Using device: {device}")
    
    train_path_eb = 'convert_root_to_torch/eb_train_data.pt'
    test_path_eb = 'convert_root_to_torch/eb_test_data.pt'
    train_model(device, train_path_eb, test_path_eb, savefile_ext = "_eb")
    
    train_path_ee = 'convert_root_to_torch/ee_train_data.pt'
    test_path_ee = 'convert_root_to_torch/ee_test_data.pt'
    train_model(device, train_path_ee, test_path_ee, savefile_ext = "_ee")

if __name__ == '__main__':
    main()
