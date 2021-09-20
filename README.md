# PCR_ML_model

A random forest classifier (RFC) for predicting the outcome of a predesigned polymerase chain reaction (PCR) aiding wetlab researchers in successful PCR design.

## Description

A tool for researchers in the biological sciences to predict PCR outcome of a designed reaction prior to purchasing reagents and running PCRs. Researchers can enter information for designed reactions into the PCR_Reaction_Collection_Form_Name.xlsx and export the colected_data sheet as a csv. Researchers can then open the feature engineering notebook for reaction prediction using the scripts in the "PredictMyPCRs" folder using the instructions below to first execute FE and then predict their reactions. After executing the predictor script, the researcher can use the output spreadsheet to view the outcome of their reaction predicted in addition to biophysical parameter estimates used by the model that may help in guiding the redesign of their PCR reactions.

## Getting Started

### Dependencies

* Everything was developed in jupyter notebook using Python version 3.6.13 on Windows 10
* Libraries: pandas, numpy, os, matplotlib, sklearn, primer3, melting, Bio.SeqUtils, datetime, xlsxwriter

### Executing program

For using the model as a reaction predictor: 
* First enter user input into the Data collection spreadsheet UI macro on the first sheet of "PCR_Reaction_Collection_Form_Name.xlsx" and save the output 'collected_data" sheet to YourRaw_collected_data.csv or whatever you would like.
* Second execute the feature engineering script using the following line: python Prediction_PCR_Feature_Engineering.py YourRaw_collected_data.csv predictMyPCRs.csv
* Third execute the prediction script on your processed PCR data: python PCR_Predictor.py predictMyPCRs.csv outcomes 
* You can label "outcomes" file will be output with todays date for tracking and contains the predicted outcomes of your PCRs and biophysical parameters calculated from your PCRs that may help guide your reaction redesign if necessary.

## Authors
Contributors names and contact info

Nicholas Cordaro - nicholas.cordaro@colorado.edu
Aaron Clauset - Aaron.clauset@colorado.edu
