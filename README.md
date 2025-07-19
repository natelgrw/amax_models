# AMAX-1: A Benchmark UV-Vis λmax Prediction Model for Aiding Compound Characterization in LC-MS

**AMAX-1** is a large scale open source dataset and machine learning framework for predicting the **UV-Vis absorption maxima (λmax)** of small molecules in various solvents. Designed to support **compound characterization workflows in LC-MS**, AMAX-1 includes both an extensive dataset and a pretrained ML model for fast, accurate λmax estimation.

## 📊 The AMAX-1 Dataset

The **AMAX-1 dataset** contains:

- **40,026 molecule–solvent combinations**, the largest λmax dataset of its kind to date
- Experimentally measured **λmax values** curated from public datasets, benchmark papers, and literature
- Hundreds of **compound and solvent descriptors** for **22,420** unique molecules and **356** unique solvents

The dataset is actively expanding with new experimental λmax values from the **Coley Research Group at MIT**, ensuring AMAX-1 remains a growing resource for optical property prediction.

## ⚙️ The AMAX-1 Model

AMAX-1 includes a fully trained, optimized machine learning model for λmax prediction. The model was selected after benchmarking over **210 algorithms** using cross-validation, hyperparameter optimization, and external test sets.

Key modeling features:

- **Input:** Concatenated vector of compound and solvent descriptors (900–1000 total features)
- **Output:** Predicted UV-Vis absorption maximum (λmax, in nm)
- **Best-performing models:** Gradient boosting (XGBoost, LightGBM), multilayer perceptrons, and ensemble regressors

The pretrained model can be used for:

- Estimating λmax for new compound–solvent pairs
- Aiding in peak assignment during LC-MS method development
- Supporting chromophore screening and spectral ML workflows