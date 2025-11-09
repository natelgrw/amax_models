# AMAX: UV-Vis Lambda Max Prediction for LC-MS

AMAX is a collection of machine learning models under active development designed to predict the λ<sub>max</sub> of chemical compounds and aid LC-MS compound characterization workflows.

Current Version: **1.0.0**

Active AMAX prediction models are accessible at this [Hugging Face Model Repository](https://huggingface.co/natelgrw/AMAX-Models). Depreciated models are available upon request at `natelgrw.tech@gmail.com`. 

This repository contains scripts for model training and processing, along with model performance metrics.

## 🧠 Available Models

For compound λ<sub>max</sub> prediction, we recommend that you use AMAX_XGB1, the model with the greatest overall prediction accuracy.

| Model | Architecture | Overall RMSE (s) | Overall MAE (s) | Overall R<sup>2</sup> | Status |
|-----|-----|-----|-----|-----|-----|
| AMAX_XGB1 | XGBoost | 56.488 | 36.005 | 0.746 | Active |
| AMAX_MLP1 | PyTorch Sequential MLP | 64.152 | 44.388 | 0.669 | Active |

All models were evaluated across rigorous scaffold, cluster, and method splits.

## 📊 The AMAX Dataset

The AMAX dataset is an open source dataset designed to assist machine learning models in small molecule UV-Vis absorption maxima (λ<sub>max</sub>) prediction and LC-MS compound characterization workflows.

The AMAX dataset contains:

- 40,013 unique molecule–environment combinations, the largest singular LC-MS retention time dataset of its kind to date
- Experimentally measured λ<sub>max</sub> values in nm, curated from public datasets, benchmark papers, and literature.

The AMAX dataset is actively expanding with new experimental retention time values from the Coley Research Group at MIT, ensuring it remains a growing resource for optical property prediction.

Additionally, AMAX includes ```.smi``` lists of 22,415 unique compounds and 356 unique solvents in the dataset for chemical descriptor calculations.

The full dataset is accessible at this [Hugging Face Dataset Repository](https://huggingface.co/datasets/natelgrw/AMAX).

AMAX is designed for use in:

- Estimating retention times for new compound–environment combinations
- Aiding in peak assignment in LC-MS method development
- Training ML models for retention time prediction under specific conditions

## 📋 Data Sources Used

Detailed information on the data sources comprising the AMAX dataset can be found in the Hugging Face repository linked above.

## ✒️ Citation

If you use an AMAX prediction model in your research, please cite:

```
@modelcollection{amaxmodels,
  title={AMAX-Models: Machine Learning Models for Molecular Absorption Wavelength Prediction},
  author={Leung, Nathan},
  institution={Coley Research Group @ MIT}
  year={2025},
  howpublished={\url{https://huggingface.co/natelgrw/AMAX-Models}},
}
```

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
