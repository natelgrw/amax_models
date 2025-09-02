# AMAX-1: UV-Vis Lambda Max Prediction for LC-MS

AMAX-1 is an open source dataset and machine learning model in active development designed to predict the (λ<sub>max</sub>) of chemical compounds and aid LC-MS compound characterization workflows.

## 📊 The AMAX-1 Dataset

The AMAX-1 dataset is an open source dataset designed to assist machine learning models in small molecule UV-Vis absorption maxima (λ<sub>max</sub>) prediction and LC-MS compound characterization workflows.

It is actively expanding with new experimental retention time values from the Coley Research Group at MIT, ensuring it remains a growing resource for optical property prediction.

Additionally, AMAX-1 includes ```.smi``` lists of 22,418 unique compounds and 356 unique solvents in the dataset for chemical descriptor calculations.

The full dataset is accessible at this [HuggingFace Repository](https://huggingface.co/datasets/natelgrw/AMAX-1).

AMAX-1 is designed for use in:

- Estimating retention times for new compound–environment combinations
- Aiding in peak assignment in LC-MS method development
- Training ML models for retention time prediction under specific conditions

The AMAX-1 dataset contains:

- 40,016 unique molecule–environment combinations, the largest singular LC-MS retention time dataset of its kind to date
- Experimentally measured λ<sub>max</sub> values in nm, curated from public datasets, benchmark papers, and literature

## 📋 Data Sources Used

Detailed information on the data sources comprising AMAX-1 can be found in the linked HuggingFace repository.

## ✒️ Citation

If you use the AMAX-1 dataset in a project, please cite the following:

```
@dataset{natelgrwamax1dataset,
  title={AMAX-1: A Benchmark Dataset for UV-Vis Lambda Max Prediction in LC-MS},
  author={Leung, Nathan},
  institution={Coley Research Group @ MIT}
  year={2025},
  howpublished={\url{https://huggingface.co/datasets/natelgrw/AMAX-1}}
}
```
