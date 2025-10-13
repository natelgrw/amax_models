#!/usr/bin/env python3
"""
train_AMAX_MLP1.py

Python script to train the AMAX_RF1 prediction model. Uses a random forest model architecture
with 312 compound + solvent RDKitdescriptors and no molecular fingerprints.
"""

import os
import sys
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import GridSearchCV, train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import pickle
import warnings
import random
warnings.filterwarnings('ignore')

# Set unbuffered output
os.environ['PYTHONUNBUFFERED'] = '1'

# Set all random seeds for full reproducibility
random.seed(42)
np.random.seed(42)
torch.manual_seed(42)
torch.cuda.manual_seed(42)
torch.cuda.manual_seed_all(42)  # For multi-GPU

class MLPRegressor(nn.Module):
    """PyTorch MLP Regressor"""
    
    def __init__(self, input_size, hidden_sizes, dropout_rate=0.2, activation='relu'):
        super(MLPRegressor, self).__init__()
        
        layers = []
        prev_size = input_size
        
        # Build hidden layers
        for hidden_size in hidden_sizes:
            layers.append(nn.Linear(prev_size, hidden_size))
            
            if activation == 'relu':
                layers.append(nn.ReLU())
            elif activation == 'tanh':
                layers.append(nn.Tanh())
            elif activation == 'leaky_relu':
                layers.append(nn.LeakyReLU())
            
            layers.append(nn.Dropout(dropout_rate))
            prev_size = hidden_size
        
        # Output layer
        layers.append(nn.Linear(prev_size, 1))
        
        self.network = nn.Sequential(*layers)
        
    def forward(self, x):
        return self.network(x)

class MLPTrainer:
    """MLP Training Class"""
    
    def __init__(self, device='cuda'):
        self.device = device
        self.scaler = StandardScaler()
        self.model = None
        self.best_params = None
        
    def load_data(self):
        """Load and preprocess data"""
        print("Loading data...")
        
        # Load training data
        train_data = pd.read_csv('data/training/amax_training.csv')
        val_data = pd.read_csv('data/validation/amax_validation.csv')
        test_data = pd.read_csv('data/testing/amax_testing.csv')
        
        # Load descriptors
        comp_descriptors = pd.read_csv('data/compounds/comp_descriptors.csv')
        solv_descriptors = pd.read_csv('data/solvents/solv_descriptors.csv')
        
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
        
        # Get features and targets
        X_train, y_train = get_descriptor_features(train_data, comp_descriptors, solv_descriptors)
        X_val, y_val = get_descriptor_features(val_data, comp_descriptors, solv_descriptors)
        X_test, y_test = get_descriptor_features(test_data, comp_descriptors, solv_descriptors)
        
        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_val_scaled = self.scaler.transform(X_val)
        X_test_scaled = self.scaler.transform(X_test)
        
        print(f"Training samples: {len(X_train_scaled)}, Features: {X_train_scaled.shape[1]}")
        
        return X_train_scaled, y_train, X_val_scaled, y_val, X_test_scaled, y_test
    
    def train_model(self, X_train, y_train, X_val, y_val, params):
        """Train MLP model with given parameters"""
        
        # Create model
        model = MLPRegressor(
            input_size=X_train.shape[1],
            hidden_sizes=params['hidden_sizes'],
            dropout_rate=params['dropout_rate'],
            activation=params['activation']
        ).to(self.device)
        
        # Create data loaders
        train_dataset = TensorDataset(
            torch.FloatTensor(X_train), 
            torch.FloatTensor(y_train)
        )
        val_dataset = TensorDataset(
            torch.FloatTensor(X_val), 
            torch.FloatTensor(y_val)
        )
        
        train_loader = DataLoader(train_dataset, batch_size=params['batch_size'], shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=params['batch_size'], shuffle=False)
        
        # Loss and optimizer
        criterion = nn.MSELoss()
        optimizer = optim.Adam(
            model.parameters(), 
            lr=params['learning_rate']
        )
        
        # Training loop
        best_val_loss = float('inf')
        patience_counter = 0
        
        for epoch in range(params['max_epochs']):
            # Training
            model.train()
            train_loss = 0.0
            
            for batch_X, batch_y in train_loader:
                batch_X, batch_y = batch_X.to(self.device), batch_y.to(self.device)
                
                optimizer.zero_grad()
                outputs = model(batch_X).squeeze()
                loss = criterion(outputs, batch_y)
                loss.backward()
                optimizer.step()
                
                train_loss += loss.item()
            
            # Validation
            model.eval()
            val_loss = 0.0
            
            with torch.no_grad():
                for batch_X, batch_y in val_loader:
                    batch_X, batch_y = batch_X.to(self.device), batch_y.to(self.device)
                    outputs = model(batch_X).squeeze()
                    loss = criterion(outputs, batch_y)
                    val_loss += loss.item()
            
            train_loss /= len(train_loader)
            val_loss /= len(val_loader)
            
            # Early stopping
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_counter = 0
            else:
                patience_counter += 1
                
            if patience_counter >= params['patience']:
                break
        
        return model
    
    def evaluate_model(self, model, X, y):
        """Evaluate model performance"""
        model.eval()
        
        with torch.no_grad():
            X_tensor = torch.FloatTensor(X).to(self.device)
            predictions = model(X_tensor).cpu().numpy().squeeze()
        
        mse = mean_squared_error(y, predictions)
        r2 = r2_score(y, predictions)
        mae = mean_absolute_error(y, predictions)
        
        return r2, mae, np.sqrt(mse)
    
    def hyperparameter_search(self, X_train, y_train, X_val, y_val):
        """Perform hyperparameter search"""
        print("Starting hyperparameter search...")
        
        # Parameter grid - COMPARABLE TO RF AND XGBOOST (32 combinations)
        param_grid = {
            'hidden_sizes': [
                (1024, 512),      # 2 hidden layers  
                (1024, 512, 256)  # 3 hidden layers
            ],
            'dropout_rate': [0.1, 0.2],           # 2 options
            'learning_rate': [0.01, 0.1],         # 2 options (higher for faster convergence)
            'batch_size': [2048, 4096],          # 2 options (much larger for maximum GPU utilization)
            'activation': ['relu', 'tanh']        # 2 options
        }
        
        # Calculate total combinations
        total_combinations = 1
        for key, values in param_grid.items():
            total_combinations *= len(values)
        print(f"Total parameter combinations: {total_combinations}")
        
        best_score = float('inf')
        best_params = None
        
        # Grid search
        from itertools import product
        
        for i, params in enumerate(product(*param_grid.values())):
            param_dict = dict(zip(param_grid.keys(), params))
            
            # Add fixed parameters (increased for better performance)
            param_dict.update({
                'max_epochs': 150,   # Increased from 50
                'patience': 25       # Increased from 10
            })
            
            print(f"Testing combination {i+1}/{total_combinations}: {param_dict}")
            sys.stdout.flush()  # Force flush
            
            try:
                model = self.train_model(X_train, y_train, X_val, y_val, param_dict)
                r2, mae, rmse = self.evaluate_model(model, X_val, y_val)
                
                # Use MAE as scoring metric (like RF and XGBoost)
                score = mae
                
                if score < best_score:
                    best_score = score
                    best_params = param_dict.copy()
                    self.model = model
                    
                print(f"  Validation MAE: {mae:.4f}, R²: {r2:.4f}")
                sys.stdout.flush()  # Force flush
                
            except Exception as e:
                print(f"  Error: {e}")
                sys.stdout.flush()  # Force flush
                continue
        
        self.best_params = best_params
        print(f"\nBest parameters: {best_params}")
        print(f"Best validation MAE: {best_score:.4f}")
        
        return best_params

def main():
    """Main training function"""
    print("AMAX MLP Training - PyTorch Implementation")
    print("=" * 50)
    
    # Check device - REQUIRE GPU
    if not torch.cuda.is_available():
        print("ERROR: CUDA is not available!")
        print("This script requires GPU training. Please run on a machine with CUDA support.")
        exit(1)
    
    device = 'cuda'
    print(f"Using device: {device}")
    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"CUDA version: {torch.version.cuda}")
    
    # Initialize trainer
    trainer = MLPTrainer(device=device)
    
    # Load data
    X_train, y_train, X_val, y_val, X_test, y_test = trainer.load_data()
    
    # Hyperparameter search
    best_params = trainer.hyperparameter_search(X_train, y_train, X_val, y_val)
    
    # Final evaluation
    print("\nFinal Model Evaluation:")
    print("-" * 30)
    
    # Validation performance
    val_r2, val_mae, val_rmse = trainer.evaluate_model(trainer.model, X_val, y_val)
    print(f"Validation  R²: {val_r2:.4f}, MAE: {val_mae:.3f}, RMSE: {val_rmse:.3f}")
    
    # Test performance
    test_r2, test_mae, test_rmse = trainer.evaluate_model(trainer.model, X_test, y_test)
    print(f"Test        R²: {test_r2:.4f}, MAE: {test_mae:.3f}, RMSE: {test_rmse:.3f}")
    
    # Save model
    os.makedirs('models', exist_ok=True)
    
    # Save PyTorch model
    torch.save({
        'model_state_dict': trainer.model.state_dict(),
        'best_params': best_params,
        'scaler': trainer.scaler,
        'input_size': X_train.shape[1]
    }, 'models/AMAX_MLP1.pth')
    
    print("Model saved to models/AMAX_MLP1.pth")

if __name__ == "__main__":
    main()
