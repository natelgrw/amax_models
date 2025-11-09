#!/usr/bin/env python3
"""
train_AMAX_MLP1.py

Author: natelgrw
Date: 11/07/2025

Trains large PyTorch MLP with GPU acceleration on compound and solvent descriptors
(excluding Morgan fingerprints) to predict lambda_max.
"""

import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from scipy import stats
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.preprocessing import StandardScaler
from pathlib import Path
import copy

warnings.filterwarnings('ignore')


# ===== Configuration ===== #


RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)

DATA_DIR = Path("../data")
MODELS_DIR = Path("../models")
PERFORMANCE_DIR = Path("../performance/AMAX_MLP1")

HYPERPARAMETERS = {
    "hidden_layers": [1024, 512, 256, 128, 64],
    "dropout_rate": 0.3,
    "use_batch_norm": True,
    "use_residual": True,
    "learning_rate": 0.001,
    "batch_size": 128,
    "max_epochs": 500,
    "early_stopping_patience": 50,
    "early_stopping_min_delta": 1e-4,
    "lr_scheduler_patience": 20,
    "lr_scheduler_factor": 0.5,
    "lr_warmup_epochs": 5,
    "weight_decay": 1e-5,
    "grad_clip_norm": 1.0,
    "optimizer": "AdamW",
    "random_state": RANDOM_SEED
}

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {DEVICE}")

SPLIT_DIRS = {
    "scaffold": DATA_DIR / "scaffold_split",
    "solvent": DATA_DIR / "solvent_split", 
    "cluster": DATA_DIR / "cluster_split"
}

SPLIT_FILES = {
    "scaffold": ["fold_1.csv", "fold_2.csv", "fold_3.csv", "fold_4.csv", "fold_5.csv"],
    "solvent": ["solvents_1.csv", "solvents_2.csv", "solvents_3.csv", "solvents_4.csv", "solvents_5.csv"],
    "cluster": ["cluster_1.csv", "cluster_2.csv", "cluster_3.csv", "cluster_4.csv", "cluster_5.csv"]
}

# ===== Model Architecture ===== #


class ResidualBlock(nn.Module):
    """
    Residual block with optional batch normalization.
    """
    
    def __init__(self, in_dim, out_dim, dropout_rate=0.3, use_batch_norm=True):
        super(ResidualBlock, self).__init__()
        
        self.linear = nn.Linear(in_dim, out_dim)
        self.batch_norm = nn.BatchNorm1d(out_dim) if use_batch_norm else nn.Identity()
        self.activation = nn.ReLU()
        self.dropout = nn.Dropout(dropout_rate)
        
        self.use_residual = (in_dim == out_dim)
        if not self.use_residual and in_dim != out_dim:
            self.projection = nn.Linear(in_dim, out_dim)
        else:
            self.projection = None
    
    def forward(self, x):
        """
        Forward pass through residual block.
        """
        identity = x
        
        out = self.linear(x)
        out = self.batch_norm(out)
        out = self.activation(out)
        out = self.dropout(out)
        
        if self.projection is not None:
            identity = self.projection(identity)
        
        if self.use_residual or self.projection is not None:
            out = out + identity
        
        return out


class LargeMLP(nn.Module):
    """
    Large multi-layer perceptron with batch normalization, dropout, and optional residual connections.
    """
    
    def __init__(self, input_dim, hidden_layers, dropout_rate=0.3, use_batch_norm=True, use_residual=False):
        super(LargeMLP, self).__init__()
        
        self.use_residual = use_residual
        
        if use_residual:
            blocks = []
            prev_dim = input_dim
            
            for hidden_dim in hidden_layers:
                blocks.append(ResidualBlock(prev_dim, hidden_dim, dropout_rate, use_batch_norm))
                prev_dim = hidden_dim
            
            self.network = nn.Sequential(*blocks)
            self.output_layer = nn.Linear(prev_dim, 1)
        else:
            layers = []
            prev_dim = input_dim
            
            for hidden_dim in hidden_layers:
                layers.append(nn.Linear(prev_dim, hidden_dim))
                if use_batch_norm:
                    layers.append(nn.BatchNorm1d(hidden_dim))
                layers.append(nn.ReLU())
                layers.append(nn.Dropout(dropout_rate))
                prev_dim = hidden_dim
            
            layers.append(nn.Linear(prev_dim, 1))
            
            self.network = nn.Sequential(*layers)
            self.output_layer = None
        
        self.apply(self._init_weights)
    
    def _init_weights(self, module):
        """
        Initialize weights for linear layers.
        """
        if isinstance(module, nn.Linear):
            nn.init.kaiming_normal_(module.weight, mode='fan_in', nonlinearity='relu')
            if module.bias is not None:
                nn.init.constant_(module.bias, 0)
    
    def forward(self, x):
        """
        Forward pass through MLP.
        """
        out = self.network(x)
        if self.output_layer is not None:
            out = self.output_layer(out)
        return out.squeeze()


# ===== Helper Functions ===== #


def load_descriptors():
    """
    Load compound and solvent descriptors (excluding fingerprints).
    """
    print("Loading descriptors...")
    
    comp_desc_path = DATA_DIR / "compounds" / "comp_descriptors.csv"
    comp_desc = pd.read_csv(comp_desc_path)
    
    comp_desc = comp_desc.rename(columns={'smiles': 'compound'})
    
    comp_fp_cols = [col for col in comp_desc.columns if 'Morgan' in col or col.startswith('Morgan_FP_')]
    comp_desc = comp_desc.drop(columns=comp_fp_cols)
    print(f"  Loaded compound descriptors: {comp_desc.shape}")
    print(f"  Excluded {len(comp_fp_cols)} Morgan fingerprint columns")
    
    solv_desc_path = DATA_DIR / "solvents" / "solv_descriptors.csv"
    solv_desc = pd.read_csv(solv_desc_path)
    
    solv_desc = solv_desc.rename(columns={'smiles': 'solvent'})
    
    solv_fp_cols = [col for col in solv_desc.columns if 'Morgan' in col or col.startswith('Morgan_FP_')]
    solv_desc = solv_desc.drop(columns=solv_fp_cols)
    print(f"  Loaded solvent descriptors: {solv_desc.shape}")
    print(f"  Excluded {len(solv_fp_cols)} Morgan fingerprint columns")
    
    return comp_desc, solv_desc


def merge_descriptors(df, comp_desc, solv_desc):
    """
    Merge compound and solvent descriptors with dataset.
    """
    df = df.merge(comp_desc, on='compound', how='left', suffixes=('', '_comp'))
    
    df = df.merge(solv_desc, on='solvent', how='left', suffixes=('', '_solv'))
    
    return df


def prepare_data(df, comp_desc, solv_desc, scaler=None, fit_scaler=False):
    """
    Prepare features and target for training.
    """
    df_merged = merge_descriptors(df, comp_desc, solv_desc)
    
    exclude_cols = ['compound', 'solvent', 'lambda_max', 'scaffold', 'source']
    feature_cols = [col for col in df_merged.columns if col not in exclude_cols]
    
    X = df_merged[feature_cols]
    y = df_merged['lambda_max']
    
    X = X.fillna(X.mean())
    
    if fit_scaler:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
    else:
        X_scaled = scaler.transform(X)
    
    return X_scaled, y.values, feature_cols, df_merged, scaler


def calculate_metrics(y_true, y_pred):
    """
    Calculate regression metrics.
    """
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    pearson_r, _ = stats.pearsonr(y_true, y_pred)
    n_test = len(y_true)
    
    return {
        'RMSE': rmse,
        'MAE': mae,
        'R2': r2,
        'Pearson_r': pearson_r,
        'N_test': n_test
    }


def create_data_loader(X, y, batch_size, shuffle=True):
    """
    Create PyTorch DataLoader.
    """
    X_tensor = torch.FloatTensor(X).to(DEVICE)
    y_tensor = torch.FloatTensor(y).to(DEVICE)
    dataset = TensorDataset(X_tensor, y_tensor)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=shuffle)
    return loader


def train_epoch(model, train_loader, criterion, optimizer, grad_clip_norm=None):
    """
    Train for one epoch with optional gradient clipping.
    """
    model.train()
    total_loss = 0.0
    
    for X_batch, y_batch in train_loader:
        optimizer.zero_grad()
        outputs = model(X_batch)
        loss = criterion(outputs, y_batch)
        loss.backward()
        
        if grad_clip_norm is not None:
            torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip_norm)
        
        optimizer.step()
        total_loss += loss.item()
    
    return total_loss / len(train_loader)


def get_lr_warmup_multiplier(epoch, warmup_epochs):
    """
    Calculate learning rate multiplier for warmup phase.
    """
    if warmup_epochs == 0 or epoch >= warmup_epochs:
        return 1.0
    return (epoch + 1) / warmup_epochs


def evaluate(model, data_loader):
    """
    Evaluate model and return predictions.
    """
    model.eval()
    all_preds = []
    all_targets = []
    
    with torch.no_grad():
        for X_batch, y_batch in data_loader:
            outputs = model(X_batch)
            all_preds.append(outputs.cpu().numpy())
            all_targets.append(y_batch.cpu().numpy())
    
    y_pred = np.concatenate(all_preds)
    y_true = np.concatenate(all_targets)
    
    return y_true, y_pred


def train_single_split(split_type, holdout_id, comp_desc, solv_desc):
    """
    Train model on a single train/test split.
    """
    print(f"\n{'='*70}")
    print(f"Training: {split_type.capitalize()} Split - Holdout {holdout_id}")
    print(f"{'='*70}")
    
    split_dir = SPLIT_DIRS[split_type]
    split_files = SPLIT_FILES[split_type]
    
    folds = []
    for i, filename in enumerate(split_files):
        fold_path = split_dir / filename
        fold_df = pd.read_csv(fold_path)
        fold_df['fold_id'] = i + 1
        folds.append(fold_df)
    
    test_df = folds[holdout_id - 1]
    train_dfs = [folds[i] for i in range(5) if i != holdout_id - 1]
    train_df = pd.concat(train_dfs, ignore_index=True)
    
    print(f"  Train size: {len(train_df):,}")
    print(f"  Test size: {len(test_df):,}")
    
    X_train, y_train, feature_cols, train_merged, scaler = prepare_data(
        train_df, comp_desc, solv_desc, scaler=None, fit_scaler=True
    )
    X_test, y_test, _, test_merged, _ = prepare_data(
        test_df, comp_desc, solv_desc, scaler=scaler, fit_scaler=False
    )
    
    train_loader = create_data_loader(X_train, y_train, HYPERPARAMETERS['batch_size'], shuffle=True)
    test_loader = create_data_loader(X_test, y_test, HYPERPARAMETERS['batch_size'], shuffle=False)
    
    input_dim = X_train.shape[1]
    model = LargeMLP(
        input_dim=input_dim,
        hidden_layers=HYPERPARAMETERS['hidden_layers'],
        dropout_rate=HYPERPARAMETERS['dropout_rate'],
        use_batch_norm=HYPERPARAMETERS['use_batch_norm'],
        use_residual=HYPERPARAMETERS['use_residual']
    ).to(DEVICE)
    
    criterion = nn.MSELoss()
    
    if HYPERPARAMETERS['optimizer'] == 'AdamW':
        optimizer = optim.AdamW(
            model.parameters(),
            lr=HYPERPARAMETERS['learning_rate'],
            weight_decay=HYPERPARAMETERS['weight_decay']
        )
    else:
        optimizer = optim.Adam(
            model.parameters(),
            lr=HYPERPARAMETERS['learning_rate'],
            weight_decay=HYPERPARAMETERS['weight_decay']
        )
    
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode='min',
        factor=HYPERPARAMETERS['lr_scheduler_factor'],
        patience=HYPERPARAMETERS['lr_scheduler_patience'],
    )
    
    print(f"  Training MLP...")
    best_val_loss = float('inf')
    best_model_state = None
    patience_counter = 0
    warmup_epochs = HYPERPARAMETERS['lr_warmup_epochs']
    min_delta = HYPERPARAMETERS['early_stopping_min_delta']
    
    for epoch in range(HYPERPARAMETERS['max_epochs']):
        if epoch < warmup_epochs:
            lr_mult = get_lr_warmup_multiplier(epoch, warmup_epochs)
            for param_group in optimizer.param_groups:
                param_group['lr'] = HYPERPARAMETERS['learning_rate'] * lr_mult
        
        train_loss = train_epoch(model, train_loader, criterion, optimizer, 
                                grad_clip_norm=HYPERPARAMETERS['grad_clip_norm'])
        
        y_true_test, y_pred_test = evaluate(model, test_loader)
        val_loss = mean_squared_error(y_true_test, y_pred_test)
        
        if epoch >= warmup_epochs:
            scheduler.step(val_loss)
        
        if val_loss < (best_val_loss - min_delta):
            best_val_loss = val_loss
            best_model_state = copy.deepcopy(model.state_dict())
            patience_counter = 0
        else:
            patience_counter += 1
        
        if (epoch + 1) % 50 == 0:
            current_lr = optimizer.param_groups[0]['lr']
            print(f"    Epoch {epoch+1}/{HYPERPARAMETERS['max_epochs']}: "
                  f"Train Loss = {train_loss:.4f}, Val RMSE = {np.sqrt(val_loss):.4f}, LR = {current_lr:.6f}")
        
        if patience_counter >= HYPERPARAMETERS['early_stopping_patience']:
            print(f"    Early stopping at epoch {epoch+1}")
            break
    
    model.load_state_dict(best_model_state)
    
    y_true, y_pred = evaluate(model, test_loader)
    metrics = calculate_metrics(y_true, y_pred)
    
    print(f"  Results:")
    print(f"    RMSE: {metrics['RMSE']:.4f}")
    print(f"    MAE: {metrics['MAE']:.4f}")
    print(f"    R²: {metrics['R2']:.4f}")
    print(f"    Pearson r: {metrics['Pearson_r']:.4f}")
    
    predictions_df = test_merged[['compound', 'solvent']].copy()
    predictions_df['split_type'] = split_type
    predictions_df['holdout_id'] = holdout_id
    predictions_df['lmax_true'] = y_true
    predictions_df['lmax_pred'] = y_pred
    
    return {
        'split_type': split_type,
        'holdout_id': holdout_id,
        'metrics': metrics,
        'predictions': predictions_df,
        'model': model,
        'scaler': scaler,
        'feature_cols': feature_cols
    }


def train_all_splits():
    """
    Train models on all splits and collect results.
    """
    print("\n" + "="*70)
    print("STARTING CROSS-VALIDATION TRAINING")
    print("="*70)
    
    comp_desc, solv_desc = load_descriptors()
    
    all_results = []
    all_predictions = []
    
    for split_type in ['scaffold', 'solvent', 'cluster']:
        for holdout_id in range(1, 6):
            result = train_single_split(split_type, holdout_id, comp_desc, solv_desc)
            
            all_results.append({
                'Split_Type': split_type,
                'Holdout': holdout_id,
                **result['metrics']
            })
            
            all_predictions.append(result['predictions'])
    
    return all_results, all_predictions, comp_desc, solv_desc


def create_summary_table(all_results):
    """
    Create summary table with averages.
    """
    results_df = pd.DataFrame(all_results)
    
    summary_rows = []
    
    for split_type in ['scaffold', 'solvent', 'cluster']:
        split_results = results_df[results_df['Split_Type'] == split_type]
        avg_row = {
            'Split_Type': f"{split_type.capitalize()} (avg)",
            'Holdout': '—',
            'RMSE': split_results['RMSE'].mean(),
            'MAE': split_results['MAE'].mean(),
            'R2': split_results['R2'].mean(),
            'Pearson_r': split_results['Pearson_r'].mean(),
            'N_test': split_results['N_test'].sum()
        }
        summary_rows.append(avg_row)
    
    overall_avg = {
        'Split_Type': 'Overall (avg)',
        'Holdout': '—',
        'RMSE': results_df['RMSE'].mean(),
        'MAE': results_df['MAE'].mean(),
        'R2': results_df['R2'].mean(),
        'Pearson_r': results_df['Pearson_r'].mean(),
        'N_test': results_df['N_test'].sum()
    }
    summary_rows.append(overall_avg)
    
    summary_df = pd.concat([results_df, pd.DataFrame(summary_rows)], ignore_index=True)
    
    return summary_df


def train_final_model(comp_desc, solv_desc):
    """
    Train final model on entire dataset.
    """
    print("\n" + "="*70)
    print("TRAINING FINAL MODEL ON ENTIRE DATASET")
    print("="*70)
    
    full_data_path = DATA_DIR / "amax_dataset.csv"
    df = pd.read_csv(full_data_path)
    print(f"  Full dataset size: {len(df):,}")
    
    X, y, feature_cols, _, scaler = prepare_data(df, comp_desc, solv_desc, scaler=None, fit_scaler=True)
    
    train_loader = create_data_loader(X, y, HYPERPARAMETERS['batch_size'], shuffle=True)
    
    input_dim = X.shape[1]
    model = LargeMLP(
        input_dim=input_dim,
        hidden_layers=HYPERPARAMETERS['hidden_layers'],
        dropout_rate=HYPERPARAMETERS['dropout_rate'],
        use_batch_norm=HYPERPARAMETERS['use_batch_norm'],
        use_residual=HYPERPARAMETERS['use_residual']
    ).to(DEVICE)
    
    criterion = nn.MSELoss()
    
    if HYPERPARAMETERS['optimizer'] == 'AdamW':
        optimizer = optim.AdamW(
            model.parameters(),
            lr=HYPERPARAMETERS['learning_rate'],
            weight_decay=HYPERPARAMETERS['weight_decay']
        )
    else:
        optimizer = optim.Adam(
            model.parameters(),
            lr=HYPERPARAMETERS['learning_rate'],
            weight_decay=HYPERPARAMETERS['weight_decay']
        )
    
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode='min',
        factor=HYPERPARAMETERS['lr_scheduler_factor'],
        patience=HYPERPARAMETERS['lr_scheduler_patience'],
    )
    
    print(f"  Training final MLP...")
    warmup_epochs = HYPERPARAMETERS['lr_warmup_epochs']
    
    for epoch in range(HYPERPARAMETERS['max_epochs']):
        if epoch < warmup_epochs:
            lr_mult = get_lr_warmup_multiplier(epoch, warmup_epochs)
            for param_group in optimizer.param_groups:
                param_group['lr'] = HYPERPARAMETERS['learning_rate'] * lr_mult
        
        train_loss = train_epoch(model, train_loader, criterion, optimizer,
                                grad_clip_norm=HYPERPARAMETERS['grad_clip_norm'])
        
        if epoch >= warmup_epochs:
            scheduler.step(train_loss)
        
        if (epoch + 1) % 50 == 0:
            current_lr = optimizer.param_groups[0]['lr']
            print(f"    Epoch {epoch+1}/{HYPERPARAMETERS['max_epochs']}: "
                  f"Train Loss = {train_loss:.4f}, LR = {current_lr:.6f}")
    
    y_true, y_pred = evaluate(model, train_loader)
    metrics = calculate_metrics(y_true, y_pred)
    
    print(f"  Training performance:")
    print(f"    RMSE: {metrics['RMSE']:.4f}")
    print(f"    MAE: {metrics['MAE']:.4f}")
    print(f"    R²: {metrics['R2']:.4f}")
    print(f"    Pearson r: {metrics['Pearson_r']:.4f}")
    
    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    model_path = MODELS_DIR / "AMAX_MLP1/AMAX_MLP1.pt"
    torch.save({
        'model_state_dict': model.state_dict(),
        'scaler': scaler,
        'feature_cols': feature_cols,
        'hyperparameters': HYPERPARAMETERS,
        'input_dim': input_dim
    }, model_path)
    print(f"\n  Model saved to: {model_path}")
    
    return model, metrics, feature_cols


def save_outputs(summary_df, all_predictions, final_metrics, feature_cols):
    """
    Save all output files.
    """
    print("\n" + "="*70)
    print("SAVING OUTPUTS")
    print("="*70)
    
    PERFORMANCE_DIR.mkdir(parents=True, exist_ok=True)
    
    # summary table
    summary_path = PERFORMANCE_DIR / "summary_table.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"  Saved summary table: {summary_path}")
    
    # all predictions
    predictions_df = pd.concat(all_predictions, ignore_index=True)
    predictions_path = PERFORMANCE_DIR / "all_predictions.csv"
    predictions_df.to_csv(predictions_path, index=False)
    print(f"  Saved all predictions: {predictions_path}")
    print(f"    Total predictions: {len(predictions_df):,}")
    
    # 4. JSON summary
    json_summary = {
        'model': 'AMAX_MLP1',
        'model_type': 'PyTorch MLP',
        'hyperparameters': {k: str(v) if isinstance(v, (list, tuple)) else v for k, v in HYPERPARAMETERS.items()},
        'n_features': len(feature_cols),
        'n_samples': 40013,
        'hardware': {
            'compute_type': 'GPU',
            'gpu': {
                'model': 'NVIDIA RTX A5000',
                'architecture': 'Ampere',
                'memory_gb': 24,
                'compute_capability': '8.6',
                'cuda_cores': 8192,
                'num_gpus': 4
            }
        },
        'cross_validation': {
            'scaffold': {
                'RMSE': summary_df[summary_df['Split_Type'] == 'Scaffold (avg)']['RMSE'].values[0],
                'MAE': summary_df[summary_df['Split_Type'] == 'Scaffold (avg)']['MAE'].values[0],
                'R2': summary_df[summary_df['Split_Type'] == 'Scaffold (avg)']['R2'].values[0],
                'Pearson_r': summary_df[summary_df['Split_Type'] == 'Scaffold (avg)']['Pearson_r'].values[0]
            },
            'solvent': {
                'RMSE': summary_df[summary_df['Split_Type'] == 'Solvent (avg)']['RMSE'].values[0],
                'MAE': summary_df[summary_df['Split_Type'] == 'Solvent (avg)']['MAE'].values[0],
                'R2': summary_df[summary_df['Split_Type'] == 'Solvent (avg)']['R2'].values[0],
                'Pearson_r': summary_df[summary_df['Split_Type'] == 'Solvent (avg)']['Pearson_r'].values[0]
            },
            'cluster': {
                'RMSE': summary_df[summary_df['Split_Type'] == 'Cluster (avg)']['RMSE'].values[0],
                'MAE': summary_df[summary_df['Split_Type'] == 'Cluster (avg)']['MAE'].values[0],
                'R2': summary_df[summary_df['Split_Type'] == 'Cluster (avg)']['R2'].values[0],
                'Pearson_r': summary_df[summary_df['Split_Type'] == 'Cluster (avg)']['Pearson_r'].values[0]
            },
            'overall': {
                'RMSE': summary_df[summary_df['Split_Type'] == 'Overall (avg)']['RMSE'].values[0],
                'MAE': summary_df[summary_df['Split_Type'] == 'Overall (avg)']['MAE'].values[0],
                'R2': summary_df[summary_df['Split_Type'] == 'Overall (avg)']['R2'].values[0],
                'Pearson_r': summary_df[summary_df['Split_Type'] == 'Overall (avg)']['Pearson_r'].values[0]
            }
        },
        'final_model': {
            'trained_on': 'full_dataset',
            'n_samples': int(summary_df[summary_df['Split_Type'] == 'Overall (avg)']['N_test'].values[0]),
            'n_features': len(feature_cols),
            'training_metrics': {
                'RMSE': final_metrics['RMSE'],
                'MAE': final_metrics['MAE'],
                'R2': final_metrics['R2'],
                'Pearson_r': final_metrics['Pearson_r']
            }
        }
    }
    
    json_path = PERFORMANCE_DIR / "summary.json"
    with open(json_path, 'w') as f:
        json.dump(json_summary, f, indent=2)
    print(f"  Saved JSON summary: {json_path}")


def print_final_summary(summary_df):
    """
    Print final summary table to console.
    """
    print("\n" + "="*70)
    print("FINAL SUMMARY TABLE")
    print("="*70)
    print(summary_df.to_string(index=False))
    print("="*70)


# ===== Main ===== #    


def main():
    """
    Main training pipeline.
    """
    print("\n" + "="*70)
    print("AMAX_MLP1 TRAINING PIPELINE")
    print("="*70)
    print(f"Device: {DEVICE}")
    print(f"Hyperparameters:")
    for key, value in HYPERPARAMETERS.items():
        print(f"  {key}: {value}")
    
    # training all splits
    all_results, all_predictions, comp_desc, solv_desc = train_all_splits()
    
    # creating summary table
    summary_df = create_summary_table(all_results)
    
    # training final model
    final_model, final_metrics, feature_cols = train_final_model(comp_desc, solv_desc)
    
    # saving all outputs
    save_outputs(summary_df, all_predictions, final_metrics, feature_cols)
    
    # printing final summary
    print_final_summary(summary_df)
    
    print("\n" + "="*70)
    print("TRAINING COMPLETE!")
    print("="*70)
    print(f"All outputs saved to: {PERFORMANCE_DIR}")
    print(f"Final model saved to: {MODELS_DIR / 'AMAX_MLP1/AMAX_MLP1.pt'}")


if __name__ == "__main__":
    main()

