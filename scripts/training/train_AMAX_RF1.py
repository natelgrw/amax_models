"""
train_AMAX_RF1.py

Python script to train the AMAX_RF1 prediction model. Uses a random forest model architecture
with 312 compound + solvent RDKitdescriptors and no molecular fingerprints.
"""

import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import joblib

# path configuration
TRAIN_PATH = "data/training/amax_training.csv"
VAL_PATH   = "data/validation/amax_validation.csv"
TEST_PATH  = "data/testing/amax_testing.csv"
COMP_DESC_PATH = "data/compounds/comp_descriptors.csv"
SOLV_DESC_PATH = "data/solvents/solv_descriptors.csv"
MODEL_DIR  = "models"
MODEL_PATH = os.path.join(MODEL_DIR, "AMAX_RF1.pkl")

# loading data
train_df = pd.read_csv(TRAIN_PATH)
val_df   = pd.read_csv(VAL_PATH)
test_df  = pd.read_csv(TEST_PATH)

# loading descriptor data
comp_desc_df = pd.read_csv(COMP_DESC_PATH)
solv_desc_df = pd.read_csv(SOLV_DESC_PATH)

def get_descriptor_features(df, comp_desc_df, solv_desc_df):
    """
    Merge compound and solvent descriptors with the dataset
    Exclude fingerprint columns and SMILES
    """
    # merging compound and solvent descriptors
    df_merged = df.merge(comp_desc_df, left_on='compound', right_on='smiles', how='left')
    df_merged = df_merged.merge(solv_desc_df, left_on='solvent', right_on='smiles', how='left', suffixes=('_comp', '_solv'))

    # Select descriptor features only (exclude SMILES, fingerprints, and target)
    exclude_cols = {"lambda_max", "compound", "solvent", "source", "smiles", "smiles_comp", "smiles_solv"}
    fingerprint_cols = [c for c in df_merged.columns if "Morgan_Fingerprint" in c]
    
    feature_cols = [
        c for c in df_merged.columns
        if c not in exclude_cols and c not in fingerprint_cols
    ]
    
    return df_merged[feature_cols], df_merged["lambda_max"]

X_train, y_train = get_descriptor_features(train_df, comp_desc_df, solv_desc_df)
X_val,   y_val   = get_descriptor_features(val_df, comp_desc_df, solv_desc_df)
X_test,  y_test  = get_descriptor_features(test_df, comp_desc_df, solv_desc_df)

print(f"Training samples: {len(X_train)}, Features: {X_train.shape[1]}")

# hyperparameter tuning
rf = RandomForestRegressor(random_state=42, n_jobs=32)
param_grid = {
    "n_estimators": [200, 500],      # 2 options
    "max_depth": [None, 20],         # 2 options  
    "min_samples_split": [2, 5],     # 2 options
    "min_samples_leaf": [1, 2],      # 2 options
    "max_features": [0.8, 1.0]      # 2 options
}

grid = GridSearchCV(
    rf,
    param_grid,
    cv=3,
    scoring="neg_mean_absolute_error",
    verbose=2,
    n_jobs=32
)

grid.fit(X_train, y_train)
print("Best parameters:", grid.best_params_)
best_rf = grid.best_estimator_

# model evaluation
def evaluate(model, X, y, label="Validation"):
    """
    Evaluate the model on the given data
    """
    pred = model.predict(X)
    r2 = r2_score(y, pred)
    mae = mean_absolute_error(y, pred)
    rmse = np.sqrt(mean_squared_error(y, pred))
    print(f"{label}  R²: {r2:.4f}, MAE: {mae:.3f}, RMSE: {rmse:.3f}")
    return r2, mae, rmse

evaluate(best_rf, X_val, y_val, "Validation")
evaluate(best_rf, X_test, y_test, "Test")

# training model on training data only
final_rf = RandomForestRegressor(**grid.best_params_, random_state=42, n_jobs=32)
final_rf.fit(X_train, y_train)

# save model
os.makedirs(MODEL_DIR, exist_ok=True)
joblib.dump(final_rf, MODEL_PATH)
print(f"Model saved to {MODEL_PATH}")