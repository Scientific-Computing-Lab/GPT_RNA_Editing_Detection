# RNA Editing Detection using GPT

## Overview
This repository contains the code and data associated with our NeurIPS 2024 submission: "Detection of RNA Editing Sites by GPT Fine-tuning." Our work focuses on leveraging GPT models, specifically GPT-3.5/4, to predict RNA editing sites, particularly involving ADAR enzymes and Alu elements. By framing the problem as both a generation and classification task, we demonstrate the transformative potential of large language models in understanding and predicting RNA editing events.

### Main Contributions:

- Generation Approach: Predict RNA edited structures given RNA sequences and secondary structures.
- Classification Approach: Predict whether specific sites within RNA sequences are edited, using a binary classification model.
- Data Augmentation: Develop robust strategies to enhance the diversity of training data.
- Threshold Adjustment: Implement optimization techniques to balance sensitivity and specificity in predictions.


## Repository Structure
- data_preparation/:
   - Code for preparing and augmenting RNA sequence data for model training.
- inference/:
   - Scripts for generating predictions using both the generative and classification models.
   - Instructions for fine-tuning GPT models on RNA editing tasks.
- figure/:
   - Code and data used to create key figures from the article.
   - Instructions for replicating the graphical representations and results shown in the paper.

   
## Getting Started
### Requirments

First, clone this repository. 

You may use the file  `environment.yml` to create anaconda environment (python 3.8) with the required packages.

### Steps to Use the environment.yml File:
#### Create the Environment:
1. Save the `environment.yml` file in your project directory, then run the following command:
   
```
conda env create -f environment.yml
```

2. Activate the Environment:
   
```
conda activate rna_editing
```

## Data Preparation
### Methodology

Data preparation differs depending on whether you are focusing on generation or classification tasks.

![methodology](https://github.com/user-attachments/assets/baa19f96-3e4b-44b8-8c44-0825420783de)

### 1. Generation
   
For the generation task, data preparation involves augmenting the sequences and then preparing them for model input.

Creation Data And Augmentation Script: This script generates the dsRNA structure and modifies RNA sequences by introducing variations while preserving secondary structures.

To run this script, navigate to the scripts/generation folder and use the following command:

```
python RNA_Structure_Data_Creation_Augmentation.py  [-h] --pair_region PAIR_REGION --output_dir OUTPUT_DIR
                      --editing_site_plus EDITING_SITE_PLUS --editing_site_minus
                      EDITING_SITE_MINUS --genome GENOME --processes PROCESSES --editing_level
                      EDITING_LEVEL --distance DISTANCE --num_augment_seq NUM_AUGMENT_SEQ
                      --min_percent_change MIN_PERCENT_CHANGE --max_percent_change
                      MAX_PERCENT_CHANGE --structural_identity PERCENT_STRUCTURAL_IDENTITY
                      --method method
```
* Example input files (pair_region,editing_site_plus,editing_site_minus) for the script can be found in the scripts/data_preparation/input_file folder
  
After running the creation and augmentation script, prepare the data for model input. This step converts the sequences into a format suitable for the model.
Run the preparation script:
```
python Model_Input_Preparation_Generation.py <input_csv_file>
```
- input_csv_file: The input file for this script should be the CSV file generated by the RNA_Structure_Generation_Augmentation.py script.

* Example input files for model (train and valid) can be found in the folder scripts/data_preparation/input_file folder

### 2. Classification

For the classification task, you have two options depending on whether you want to use full or partial context around RNA editing sites. After choosing the context, you will prepare the data and balance it using a selection script.

To run this script, navigate to the scripts/classification folder and use the following command:
```
python Classification_Data_Creation.py [-h] --pair_region PAIR_REGION --output_dir OUTPUT_DIR
                                        --editing_site_plus EDITING_SITE_PLUS
                                        --editing_site_minus EDITING_SITE_MINUS --genome GENOME
                                        --editing_level EDITING_LEVEL --context FULL\PARTIAL
```
* Example input files for the script can be found in the scripts/data_preparation/input_samples folder
  
Run the preparation script:
```
python Model_Input_Preparation_Classification.py <input_csv_file>
```
- input_csv_file: The input file is the output CSV from Classification_Data_Creation.py.

Finally, balance the prepared data using the selection script.
Run the selector script:
```
python Balanced_Data_Selector.py <input_file.jsonl> <output_file.jsonl>
```
input_file: The input JSONL file should be the one created by Model_Input_Preparation_Classification.py.


Note: Each of these scripts comes with various parameters that can be adjusted according to your needs. You can use --help with any script to see all available options and required parameters:

```
python <script_name>.py --help
```
### Baseline Non-AI Model
In the classification approach, we include a Baseline_Non_AI script that tests RNA editing site predictions based on established biological rules, without the use of AI. This script serves as a benchmark, allowing for a comparison between traditional rule-based predictions and AI-driven models.

Run the Baseline Non-AI Model script:

```
python Rule_Based_RNA_Editing_Prediction.py <input_file>
```
- The <input_file> should be the file generated from the Model_Input_Preparation_Classification script.


## Inference
Building on the methodology you've chosen (classification or generation), the inference process is straightforward. Depending on whether you're using a classification or generation approach, you'll find the relevant scripts in the scripts/inferencing folder.

### Running Inference
- Classification Methodology: If you opted for the classification-based approach, use the following command:
```
  python infer_classification.py <input_file> <output_file> <temperature>
```
Input: The <input_file> is the file created in the Balanced_Data_Selector.py step.

- Generation Methodology: For those utilizing the generation approach, run:
```
   python infer_generation.py <input_file> <output_file> <temperature>
```
Input: The <input_file> is the file generated by the generation/Model_Input_Preparation_Generation.py script during the data preparation step.

## Training and validation performance
Below is an illustration showing the performance of fine-tuned GPT-3.5 models in detecting RNA editing sites. This visual captures the model's accuracy and loss over the training epochs, providing insight into how effectively the model learns and generalizes.

![train_valid_trends (1)](https://github.com/user-attachments/assets/a91cbcc4-04c7-4105-890c-2f70d2e454bc)

To generate this illustration, navigate to the figure/Training_Validation_Performance directory and run the provided scripts:

```
python Training_Validation_Performance_Plot.py <Generation.csv> <Classification-PARTIAL.csv> <Classification-FULL.csv>
```
- The input files (Generation.csv, Classification-PARTIAL.csv, Classification-FULL.csv) are provided within the figure/fig5 directory.

## Threshold Adjustment
The following illustration depicts the threshold adjustment process, which optimizes the classification boundary for detecting RNA editing sites. This adjustment ensures that the model achieves the best balance between sensitivity and specificity by fine-tuning the decision threshold based on log probabilities.

![Threshold Adjustment (1)](https://github.com/user-attachments/assets/e0e06963-96d4-4d30-aeee-b7f18a95ca47)

To create this illustration, navigate to the figure/Threshold_Optimization directory and run the provided scripts:

```
python Threshold_Optimization_Visualization.py <input_file.csv>
```
- <input_file>:
   - full context -  figure/Threshold_Optimization/full_context.txt
   - partial context - figure/Threshold_Optimization/partial_context.txt

## Tissue-Specific Analysis

In this part of our study, we conducted tissue-specific analysis to investigate potential differences in RNA editing patterns based on varying ADAR expression levels across tissues. This analysis covers six tissues, each with distinct ADAR and RNA-binding protein profiles: Muscle Skeletal, Whole Blood, Artery Tibial, Brain Cerebellum, Esophagus Muscularis and Brain Spinal Cord.

For each tissue, we generated and analyzed RNA editing patterns and used classification-based models to predict whether adenosines were edited or not. The tissue-specific visualization graph illustrates the threshold optimization for each tissue.

![image](https://github.com/user-attachments/assets/43695c23-c9bb-4716-a5c2-383766b20ccc)

To replicate this analysis for a specific tissue, navigate to the figure/Tissue_Specific_Analysis directory and run the provided scripts:

```
python Threshold_Optimization_Tissues_Visualization.py <{TissueName}\{TissueName}_Thresholds.txt>
```

For example, to visualize Brain Cerebellum analysis:
```
python Threshold_Optimization_Tissues_Visualization.py BrainCerebellum\BrainCerebellum_Thresholds.txt
```

1. {TissueName}_Training: Contains training statistics for each tissue.
2. {TissueName}_Thresholds: Contains threshold data used for illustrates the threshold optimization for each tissue.
3. {TissueName}_output_meta_stats: Metadata statistics of the analysis for each tissue.
4. {TissueName}.png: Graphical representation of threshold adjustments.

The visualization script takes the threshold data as input and produces the tissue-specific plot demonstrating the effect of threshold adjustments.
