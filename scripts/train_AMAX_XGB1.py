#!/usr/bin/env python3
"""
train_AMAX_XGB1.py

Author: natelgrw
Date: 11/07/2025

Trains XGBoost model with GPU acceleration on compound and solvent descriptors
(excluding Morgan fingerprints) to predict lambda_max.
"""

import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
import xgboost as xgb
from scipy import stats
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from pathlib import Path

warnings.filterwarnings('ignore')

# ===== Configuration ===== #

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

DATA_DIR = Path("../data")
MODELS_DIR = Path("../models")
PERFORMANCE_DIR = Path("../performance/AMAX_XGB1")


HYPERPARAMETERS = {
    "n_estimators": 1000,
    "max_depth": 9,
    "learning_rate": 0.05,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "tree_method": "hist",
    "device": "cuda",
    "random_state": RANDOM_SEED,
    "verbosity": 1
}

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


def prepare_data(df, comp_desc, solv_desc):
    """
    Prepare features and target for training.
    """
    df_merged = merge_descriptors(df, comp_desc, solv_desc)
    
    exclude_cols = ['compound', 'solvent', 'lambda_max', 'scaffold', 'source']
    feature_cols = [col for col in df_merged.columns if col not in exclude_cols]
    
    X = df_merged[feature_cols]
    y = df_merged['lambda_max']
    
    X = X.fillna(X.mean())
    
    return X, y, feature_cols, df_merged


def calculate_metrics(y_true, y_pred):
    """Calculate regression metrics."""
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
    
    X_train, y_train, feature_cols, train_merged = prepare_data(train_df, comp_desc, solv_desc)
    X_test, y_test, _, test_merged = prepare_data(test_df, comp_desc, solv_desc)
    
    print(f"  Training XGBoost model...")
    model = xgb.XGBRegressor(**HYPERPARAMETERS)
    model.fit(X_train, y_train, verbose=False)
    
    y_pred = model.predict(X_test)
    
    metrics = calculate_metrics(y_test, y_pred)
    
    print(f"  Results:")
    print(f"    RMSE: {metrics['RMSE']:.4f}")
    print(f"    MAE: {metrics['MAE']:.4f}")
    print(f"    R²: {metrics['R2']:.4f}")
    print(f"    Pearson r: {metrics['Pearson_r']:.4f}")
    
    predictions_df = test_merged[['compound', 'solvent']].copy()
    predictions_df['split_type'] = split_type
    predictions_df['holdout_id'] = holdout_id
    predictions_df['lmax_true'] = y_test.values
    predictions_df['lmax_pred'] = y_pred

    return {
        'split_type': split_type,
        'holdout_id': holdout_id,
        'metrics': metrics,
        'predictions': predictions_df,
        'model': model,
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
    
    X, y, feature_cols, _ = prepare_data(df, comp_desc, solv_desc)
    
    print(f"  Training final XGBoost model...")
    model = xgb.XGBRegressor(**HYPERPARAMETERS)
    model.fit(X, y, verbose=False)
    
    y_pred = model.predict(X)
    metrics = calculate_metrics(y, y_pred)
    
    print(f"  Training performance:")
    print(f"    RMSE: {metrics['RMSE']:.4f}")
    print(f"    MAE: {metrics['MAE']:.4f}")
    print(f"    R²: {metrics['R2']:.4f}")
    print(f"    Pearson r: {metrics['Pearson_r']:.4f}")
    
    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    model_path = MODELS_DIR / "AMAX_XGB1/AMAX_XGB1.json"
    model.save_model(model_path)
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
    
    # JSON summary
    json_summary = {
        'model': 'AMAX_XGB1',
        'model_type': 'XGBoost',
        'hyperparameters': HYPERPARAMETERS,
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
    print("AMAX_XGB1 TRAINING PIPELINE")
    print("="*70)
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
    print(f"Final model saved to: {MODELS_DIR / 'AMAX_XGB1.json'}")


if __name__ == "__main__":
    main()

