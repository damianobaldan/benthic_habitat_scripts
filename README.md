# benthic_habitat_scripts

Correlative benthic habitat models for the NECCTON project. Species occurrence data are fitted with a set of statistical models. Models are then evaluated with a spatial blocks approach. An ensemble is generated and used for projecting the habitat suitability for the species in the geographic space.  
This example uses Posidonia oceanica occurrence data.

For a better visualization see: https://damianobaldan.github.io/benthic_habitat_scripts/


# Usage

To fit, evaluate, and project the species distribution models, open and run benthic_habitat_scripts.Rmd.

## Data Source

The data required to run the Notebook are:

- species occurrence data, formatted as georeferenced point shapefile (e.g. .shp file)
- envrionmental predictors, formatted as georeferenced rasters (e.g. .tif file)


## Metrics

The performance of the super-resolution model is evaluated using the following metrics:

### Area under the Receiver Operating Curve (AUC)
AUC is a metric used to evaluate the performance of a binary classification model. 
AUC represents represents the probability that the model, if given a randomly chosen positive and negative 
example, will rank the positive higher than the negative.
AUC varies between 0.5 and 1, with 0.5 representing a random classifier and 1 representing a perfect classifier.

### True Skill Statistics (TSS)
TSS is a metric used to evaluate the performance of species distribution models, 
particularly in the context of presence/absence predictions. 
It is calculated from a confusion matrix based on sensitivity (true positive rate) and specificity 
(true negative rate), and is defined as sensitivity + specificity - 1. 
TSS values range from -1 to +1, with higher values indicating better model performance.

TSS relies on a confusion matrix, which summarizes the model's predictions against actual observations. It includes:

- True Positives (TP): Correctly predicted presences.
- True Negatives (TN): Correctly predicted absences.
- False Positives (FP): Incorrectly predicted presences (also known as commission errors).
- False Negatives (FN): Incorrectly predicted absences (also known as omission errors)
- Sensitivity (also known as Recall or True Positive Rate): TP / (TP + FN). It represents the proportion of actual presences that were correctly predicted. 
- Specificity (also known as True Negative Rate): TN / (TN + FP). It represents the proportion of actual absences that were correctly predicted.

### Boyce's index
The Boyce index is a occurrence-only index and measures how much model predictions differ from random distribution 
of the observed occurrences across the prediction gradients.
It is thus an appropriate metric in the case of occurrence-only models. It is continuous and varies between -1 
and +1. Positive values indicate a model which present predictions are consistent with the distribution 
of occurrences in the evaluation dataset, values close to zero mean that the model is not different from a 
random model, negative values indicate counter predictions, i.e., predicting poor quality areas where 
occurrences are more frequent.

## List of Dependencies
- R version 4.4.1

## Citations and Links
[Baldan et al., 2024](https://doi.org/10.1111/ddi.13922)













