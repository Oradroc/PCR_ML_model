# PCR_ML_model

A random forest classifier (RFC) for predicting the outcome of a predesigned polymerase chain reaction (PCR) aiding wetlab researchers in successful PCR design.

## Description

A tool for researchers in the biological sciences to predict PCR outcome of a designed reaction prior to purchasing reagents and running PCRs. Researchers can enter information for designed reactions into the PCR_Reaction_Collection_Form_Name.xlsx and export the colected_data sheet as a csv. Researchers can then open the feature engineering notebook for reaction prediction (PREDICTION_PCR_Feature_Engineering_7-10-20_2-8-21.ipynb) and change the path, filename, and export name and run the notebook. Following feature engineering, researchers can alter the filename, path, and export name in the PCR_Predictor jupyter notebook. After executing the predictor notebook, the researcher can use the output spreadsheet to view the outcome of their reaction predcited in addition to biophysical parameter estimates used by the model that may help in guiding the redesign their PCR reaction.

## Getting Started

### Dependencies

* Everything was developed in jupyter notebook using Python version 3.6.13 on Windows 10
* Libraries: Pandas, numpy, os, matplotlib, sklearn, primer3, melting, Bio.SeqUtils, datetime

### Installing

* Modify file paths for input and output files and input/ouput filenames for user data.

### Executing program

For using the model as a reaction predictor: 
* First enter user input into the Data collection spreadsheet and save 'collected_data" sheet.
* Change filepath and filename in PREDICTION_PCR_Feature_Engineering_7-10-20_2-8-21.ipynb and run notebook for reaction processing.
* Change filepath, filename (processed PCR reactions name), and export name in PCR_Predictor.ipynb and run the notebook to obtain predictions and useful biophyscial parameters for reaction design.


## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors
Contributors names and contact info

Nicholas Cordaro - nicholas.cordaro@colorado.edu
Aaron Clauset - Aaron.clauset@colorado.edu
