#!/usr/bin/env python3
"""
extract_AMAX_XGB1.py

A comprehensive Python script to extract all graphs and data needed for publication
for the AMAX XGBoost model. Generates publication-ready figures and tables.
"""

import os
import sys
import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import xgboost as xgb
import warnings
warnings.filterwarnings('ignore')

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

def get_descriptor_features(df, comp_desc_df, solv_desc_df):
    """Merge compound and solvent descriptors with the dataset"""
    df_merged = df.merge(comp_desc_df, left_on='compound', right_on='smiles', how='left')
    df_merged = df_merged.merge(solv_desc_df, left_on='solvent', right_on='smiles', how='left', suffixes=('_comp', '_solv'))

    exclude_cols = {"lambda_max", "compound", "solvent", "source", "smiles", "smiles_comp", "smiles_solv"}
    fingerprint_cols = [c for c in df_merged.columns if "Morgan_Fingerprint" in c]
    
    feature_cols = [
        c for c in df_merged.columns
        if c not in exclude_cols and c not in fingerprint_cols
    ]
    
    return df_merged[feature_cols], df_merged["lambda_max"]

def load_xgb_model():
    """Load the trained XGBoost model"""
    model = joblib.load('models/AMAX_XGB1.pkl')
    return model

def evaluate_model(model, X, y):
    """Evaluate XGBoost model and return predictions"""
    predictions = model.predict(X)
    
    r2 = r2_score(y, predictions)
    mae = mean_absolute_error(y, predictions)
    rmse = np.sqrt(mean_squared_error(y, predictions))
    
    return r2, mae, rmse, predictions

def create_prediction_vs_actual_plot(y_true, y_pred, title, filename):
    """Create prediction vs actual scatter plot"""
    plt.figure(figsize=(8, 8))
    
    # calculate metrics
    r2 = r2_score(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    
    # create scatter plot
    plt.scatter(y_true, y_pred, alpha=0.6, s=20, color='darkorange')
    
    # add perfect prediction line
    min_val = min(min(y_true), min(y_pred))
    max_val = max(max(y_true), max(y_pred))
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, label='Perfect Prediction')
    
    # add regression line
    z = np.polyfit(y_true, y_pred, 1)
    p = np.poly1d(z)
    plt.plot(y_true, p(y_true), 'g-', linewidth=2, alpha=0.8, label=f'Regression Line (R²={r2:.4f})')
    
    plt.xlabel('Experimental λmax (nm)', fontsize=12)
    plt.ylabel('Predicted λmax (nm)', fontsize=12)
    plt.title(f'{title}\nR² = {r2:.4f}, MAE = {mae:.3f} nm, RMSE = {rmse:.3f} nm', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # set equal aspect ratio
    plt.axis('equal')
    
    plt.tight_layout()
    plt.savefig(f'performance/AMAX_XGB1/{filename}.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_residual_plot(y_true, y_pred, title, filename):
    """Create residual plot"""
    residuals = y_pred - y_true
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # residuals vs predicted
    ax1.scatter(y_pred, residuals, alpha=0.6, s=20, color='darkorange')
    ax1.axhline(y=0, color='red', linestyle='--', linewidth=2)
    ax1.set_xlabel('Predicted λmax (nm)', fontsize=12)
    ax1.set_ylabel('Residuals (nm)', fontsize=12)
    ax1.set_title('Residuals vs Predicted', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # histogram of residuals
    ax2.hist(residuals, bins=50, alpha=0.7, color='darkorange', edgecolor='black')
    ax2.axvline(x=0, color='red', linestyle='--', linewidth=2)
    ax2.set_xlabel('Residuals (nm)', fontsize=12)
    ax2.set_ylabel('Frequency', fontsize=12)
    ax2.set_title('Distribution of Residuals', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    plt.suptitle(f'{title} - Residual Analysis', fontsize=16)
    plt.tight_layout()
    plt.savefig(f'performance/AMAX_XGB1/{filename}.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_error_distribution_plot(y_true, y_pred, title, filename):
    """Create error distribution plot"""
    errors = np.abs(y_pred - y_true)
    
    plt.figure(figsize=(10, 6))
    
    # create histogram
    plt.hist(errors, bins=50, alpha=0.7, color='darkorange', edgecolor='black')
    
    # add statistics
    mean_error = np.mean(errors)
    median_error = np.median(errors)
    std_error = np.std(errors)
    
    plt.axvline(mean_error, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_error:.2f} nm')
    plt.axvline(median_error, color='blue', linestyle='--', linewidth=2, label=f'Median: {median_error:.2f} nm')
    
    plt.xlabel('Absolute Error (nm)', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title(f'{title} - Error Distribution\nMean: {mean_error:.2f} nm, Std: {std_error:.2f} nm', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'performance/AMAX_XGB1/{filename}.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_feature_importance_plot(model, feature_names, X_train, y_train, title, filename, top_n=20):
    """Create feature importance plot for top N features"""
    # use built-in XGB feature importance (much faster than permutation)
    importance = model.feature_importances_
    
    # get top N features
    top_indices = np.argsort(importance)[-top_n:]
    top_features = [feature_names[i] for i in top_indices]
    top_importance = importance[top_indices]
    
    # create the plot
    plt.figure(figsize=(12, 8))
    
    y_pos = np.arange(len(top_features))
    bars = plt.barh(y_pos, top_importance, color='darkorange', alpha=0.7)
    
    plt.yticks(y_pos, top_features)
    plt.xlabel('Feature Importance', fontsize=12)
    plt.title(f'{title} - Top {top_n} Feature Importance', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3, axis='x')
    
    # add value labels on bars
    for i, (bar, importance_val) in enumerate(zip(bars, top_importance)):
        width = bar.get_width()
        plt.text(width + width*0.01, bar.get_y() + bar.get_height()/2, 
                f'{importance_val:.4f}', ha='left', va='center', fontsize=9)
    
    plt.gca().invert_yaxis()
    
    plt.tight_layout()
    
    plt.savefig(f'performance/AMAX_XGB1/{filename}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Feature importance plot saved: {filename}.png")
    
    return top_features, top_importance

def get_descriptor_category(feature_name):
    """Categorize descriptor based on name"""
    feature_lower = feature_name.lower()
    
    if 'molecular' in feature_lower or 'mol' in feature_lower:
        return 'Molecular Properties'
    elif 'logp' in feature_lower or 'logp' in feature_lower:
        return 'Lipophilicity'
    elif 'tpsa' in feature_lower or 'polar' in feature_lower:
        return 'Polar Surface Area'
    elif 'molar' in feature_lower or 'mass' in feature_lower:
        return 'Molecular Weight'
    elif 'aromatic' in feature_lower or 'ring' in feature_lower:
        return 'Aromaticity'
    elif 'charge' in feature_lower or 'electro' in feature_lower:
        return 'Electronic Properties'
    elif 'solvent' in feature_lower or '_solv' in feature_lower:
        return 'Solvent Properties'
    elif 'compound' in feature_lower or '_comp' in feature_lower:
        return 'Compound Properties'
    else:
        return 'Other'

def main():
    """Main function to extract all publication data for XGBoost model"""
    print("AMAX XGBoost Model - Publication Data Extraction")
    print("=" * 60)
    
    # create performance directories
    os.makedirs('performance', exist_ok=True)
    os.makedirs('performance/AMAX_XGB1', exist_ok=True)
    
    # load data
    print("Loading data...")
    train_data = pd.read_csv('data/training/amax_training.csv')
    val_data = pd.read_csv('data/validation/amax_validation.csv')
    test_data = pd.read_csv('data/testing/amax_testing.csv')
    
    # load descriptors
    comp_descriptors = pd.read_csv('data/compounds/comp_descriptors.csv')
    solv_descriptors = pd.read_csv('data/solvents/solv_descriptors.csv')
    
    # get features and targets
    X_train, y_train = get_descriptor_features(train_data, comp_descriptors, solv_descriptors)
    X_val, y_val = get_descriptor_features(val_data, comp_descriptors, solv_descriptors)
    X_test, y_test = get_descriptor_features(test_data, comp_descriptors, solv_descriptors)
    
    print(f"Training samples: {len(X_train)}, Features: {X_train.shape[1]}")
    print(f"Validation samples: {len(X_val)}")
    print(f"Test samples: {len(X_test)}")
    
    # load XGBoost model
    print("\nLoading XGBoost model...")
    model = load_xgb_model()
    print(f"Model loaded: {model.get_params()['n_estimators']} boosting rounds")
    
    # evaluate on validation set
    print("\nEvaluating on validation set...")
    val_r2, val_mae, val_rmse, val_pred = evaluate_model(model, X_val, y_val)
    val_results = {'r2': val_r2, 'mae': val_mae, 'rmse': val_rmse}
    print(f"Validation: R²={val_r2:.4f}, MAE={val_mae:.3f}, RMSE={val_rmse:.3f}")
    
    # evaluate on test set
    print("\nEvaluating on test set...")
    test_r2, test_mae, test_rmse, test_pred = evaluate_model(model, X_test, y_test)
    test_results = {'r2': test_r2, 'mae': test_mae, 'rmse': test_rmse}
    print(f"Test: R²={test_r2:.4f}, MAE={test_mae:.3f}, RMSE={test_rmse:.3f}")
    
    # generate publication plots
    print("\nGenerating publication plots...")
    
    # 1. prediction vs actual plots
    create_prediction_vs_actual_plot(y_val, val_pred, "AMAX_XGB1 - Validation Set", "validation/pred_vs_actual")
    create_prediction_vs_actual_plot(y_test, test_pred, "AMAX_XGB1 - Test Set", "testing/pred_vs_actual")
    
    # 2. residual plots
    create_residual_plot(y_val, val_pred, "AMAX_XGB1 - Validation Set", "validation/residuals")
    create_residual_plot(y_test, test_pred, "AMAX_XGB1 - Test Set", "testing/residuals")
    
    # 3. error distribution plots
    create_error_distribution_plot(y_val, val_pred, "AMAX_XGB1 - Validation Set", "validation/error_dist")
    create_error_distribution_plot(y_test, test_pred, "AMAX_XGB1 - Test Set", "testing/error_dist")

    # 4. feature importance plot
    print("\nGenerating feature importance plot...")
    top_features, top_importance = create_feature_importance_plot(model, X_train.columns, X_train, y_train, "AMAX_XGB1", "feature_importance", top_n=20)
    
    # save model predictions for further analysis
    print("\nSaving predictions...")
    val_predictions_df = pd.DataFrame({
        'compound': val_data['compound'],
        'solvent': val_data['solvent'],
        'actual': y_val,
        'predicted': val_pred,
        'error': np.abs(val_pred - y_val)
    })
    val_predictions_df.to_csv('performance/AMAX_XGB1/validation/val_predictions.csv', index=False)
    
    test_predictions_df = pd.DataFrame({
        'compound': test_data['compound'],
        'solvent': test_data['solvent'],
        'actual': y_test,
        'predicted': test_pred,
        'error': np.abs(test_pred - y_test)
    })
    test_predictions_df.to_csv('performance/AMAX_XGB1/testing/test_predictions.csv', index=False)
    
    print("\n" + "="*50)
    print("AMAX_XGB1 PUBLICATION DATA EXTRACTION COMPLETE")
    print("="*50)

if __name__ == "__main__":
    main()
