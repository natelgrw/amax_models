# AMAX: UV-Vis Lambda Max Prediction for LC-MS

AMAX is a collection of machine learning models under active development designed to predict the λ<sub>max</sub> of chemical compounds and aid LC-MS compound characterization workflows.

This repository contains scripts for model training and processing, along with model performance metrics.

## 🤖 Available Models

| Model | Framework | Architecture | R² Score | MAE (nm) | RMSE (nm) | Status |
|-------|-----------|--------------|----------|----------|-----------|---------|
| **AMAX_XGB1** | XGBoost | Gradient Boosting (500 estimators) | 0.9084 | 17.682 | 35.507 | Active |
| **AMAX_RF1** | Scikit-Learn | Random Forest (500 trees) | 0.9035 | 18.601 | 36.441 | Active |
| **AMAX_MLP1** | PyTorch | Sequential NN (1024 -> 512) | 0.8913 | 23.956 | 38.68 | Active |

Active AMAX prediction models are accessible at this [Hugging Face Model Repository](https://huggingface.co/datasets/natelgrw/AMAX). Depreciated models are available upon request at `natelgrw.tech@gmail.com`. 

All models utilize 312 RDKit molecular descriptors combining both compound and solvent features, trained on a random data split of 32,010 training samples with 4,001 validation and 4,002 test samples. Each model has been retrained to eliminate data leakage and ensure robust performance evaluation.

## 📊 The AMAX Dataset

The AMAX dataset is an open source dataset designed to assist machine learning models in small molecule UV-Vis absorption maxima (λ<sub>max</sub>) prediction and LC-MS compound characterization workflows.

The AMAX dataset contains:

- 40,016 unique molecule–environment combinations, the largest singular LC-MS retention time dataset of its kind to date
- Experimentally measured λ<sub>max</sub> values in nm, curated from public datasets, benchmark papers, and literature.

The AMAX dataset is actively expanding with new experimental retention time values from the Coley Research Group at MIT, ensuring it remains a growing resource for optical property prediction.

Additionally, AMAX includes ```.smi``` lists of 22,418 unique compounds and 356 unique solvents in the dataset for chemical descriptor calculations.

The full dataset is accessible at this [Hugging Face Dataset Repository](https://huggingface.co/datasets/natelgrw/AMAX).

AMAX is designed for use in:

- Estimating retention times for new compound–environment combinations
- Aiding in peak assignment in LC-MS method development
- Training ML models for retention time prediction under specific conditions

## 📋 Data Sources Used

Detailed information on the data sources comprising the AMAX dataset can be found in the Hugging Face repository linked above.

## ✒️ Citation

If you use the AMAX dataset in a project, please cite the following:

```
@dataset{natelgrwamaxdataset,
  title={AMAX: A Benchmark Dataset for UV-Vis Lambda Max Prediction in LC-MS},
  author={Leung, Nathan},
  institution={Coley Research Group @ MIT}
  year={2025},
  howpublished={\url{https://huggingface.co/datasets/natelgrw/AMAX}}
}
```
