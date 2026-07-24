import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader
import time
import math

class SoftIsoSumNetwork(nn.Module):
    def __init__(self, tkel_dim=11, puppi_dim=14, hidden_dims=[64, 32]):
        super().__init__()
        
        # Build shared MLP for candidates
        layers = []
        input_dim = tkel_dim + puppi_dim
        for h_dim in hidden_dims:
            layers.append(nn.Linear(input_dim, h_dim))
            layers.append(nn.ReLU())
            input_dim = h_dim
        layers.append(nn.Linear(input_dim, 1))
        
        self.mlp = nn.Sequential(*layers)
        
        # Learnable threshold and scale for BCE conversion
        # We want: logit = (threshold - iso_sum) * scale
        self.threshold = nn.Parameter(torch.tensor([5.0]))
        self.scale_raw = nn.Parameter(torch.tensor([0.0])) # F.softplus(0) ~ 0.69
        
    def forward(self, tkel_norm, puppi_norm, puppi_pt_unnorm):
        """
        tkel_norm: (batch, 11)
        puppi_norm: (batch, 25, 14)
        puppi_pt_unnorm: (batch, 25) - physical pT to sum over
        """
        batch_size = tkel_norm.shape[0]
        
        # Expand tkel: (batch, 25, 11)
        tkel_expanded = tkel_norm.unsqueeze(1).expand(-1, 25, -1)
        
        # Concat: (batch, 25, 11 + 14)
        combined = torch.cat([tkel_expanded, puppi_norm], dim=2)
        
        # MLP -> Weights: (batch, 25)
        weights_logit = self.mlp(combined).squeeze(-1)
        weights = torch.sigmoid(weights_logit)
        
        # Weighted sum of physical pT
        iso_sum = torch.sum(weights * puppi_pt_unnorm, dim=1) # (batch,)
        
        # Convert to BCE logit
        scale = F.softplus(self.scale_raw)
        logits = (self.threshold - iso_sum) * scale
        
        return logits, iso_sum, weights

def prepare_dataloaders(train_path, test_path, batch_size=1024):
    print("Loading data...")
    train_data = torch.load(train_path, weights_only=True)
    test_data = torch.load(test_path, weights_only=True)
    
    x_train, y_train = train_data['x'], train_data['y'].float()
    x_test, y_test = test_data['x'], test_data['y'].float()
    
    # Split into tkel and puppi
    tkel_train = x_train[:, 0:11]
    puppi_train_flat = x_train[:, 11:]
    puppi_train = puppi_train_flat.reshape(-1, 25, 14)
    
    tkel_test = x_test[:, 0:11]
    puppi_test_flat = x_test[:, 11:]
    puppi_test = puppi_test_flat.reshape(-1, 25, 14)
    
    # Extract unnormalized pT (which is the 0-th feature of the 14 puppi features)
    puppi_pt_train = puppi_train[:, :, 0].clone()
    puppi_pt_test = puppi_test[:, :, 0].clone()
    
    # Compute Normalization statistics on Train only
    tkel_mean = tkel_train.mean(dim=0)
    tkel_std = tkel_train.std(dim=0)
    tkel_std[tkel_std < 1e-6] = 1.0 # Prevent division by zero
    
    # Normalize all 25 candidates uniformly
    puppi_train_reshaped = puppi_train.reshape(-1, 14)
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
    
    train_path = 'convert_root_to_torch/train_data.pt'
    test_path = 'convert_root_to_torch/test_data.pt'
    
    # Increased batch size for better multi-core throughput
    train_loader, test_loader, norm_stats = prepare_dataloaders(train_path, test_path, batch_size=8192)
    torch.save(norm_stats, 'norm_stats.pt')
    
    model = SoftIsoSumNetwork(tkel_dim=11, puppi_dim=14, hidden_dims=[64, 32]).to(device)
    
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    criterion = nn.BCEWithLogitsLoss()
    
    epochs = 20
    best_acc = 0.0
    
    print("Starting training...")
    for epoch in range(epochs):
        t0 = time.time()
        
        train_loss, train_acc = train_epoch(model, train_loader, optimizer, criterion, device)
        val_loss, val_acc, avg_iso_sig, avg_iso_bkg = eval_epoch(model, test_loader, criterion, device)
        
        t1 = time.time()
        print(f"Epoch {epoch+1:02d}/{epochs} | Time: {t1-t0:.1f}s | "
              f"Train Loss: {train_loss:.4f} Acc: {train_acc:.4f} | "
              f"Val Loss: {val_loss:.4f} Acc: {val_acc:.4f} | "
              f"Iso (Sig/Bkg): {avg_iso_sig:.2f} / {avg_iso_bkg:.2f}")
              
        if val_acc > best_acc:
            best_acc = val_acc
            torch.save(model.state_dict(), 'best_soft_iso_model.pt')
            
    print(f"Training complete. Best Validation Accuracy: {best_acc:.4f}")
    
    # Print final learned parameters
    model.load_state_dict(torch.load('best_soft_iso_model.pt', weights_only=True))
    print(f"Learned Threshold (Iso Cut): {model.threshold.item():.3f}")
    print(f"Learned Scale: {F.softplus(model.scale_raw).item():.3f}")

if __name__ == '__main__':
    main()
