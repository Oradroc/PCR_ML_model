#Adjust Data
#python PCR_ML_FE_Prediction.py Raw.csv EngineeredData Predictions
#import libraries
import pandas as pd
import numpy as np
import os
import sys

#import data

RawData = sys.argv[1]
outFile = sys.argv[2]

import datetime
now = datetime.datetime.now()
now = str(now)
date = now.split(" ")

Predictions = sys.argv[3]+date[0]+'.xlsx'
PCR_Data = pd.read_csv(RawData)

print("\nFeature Engineering...")
#Data Setup
#Drop unnecessary columns
remove=["timestamp","user","email","species","polymerase_catalog_number","pcrid",
        "template_type","polymerase_brand","dna_amt_units","band_size_units",'customized_salts']
PCR_Data.drop(remove,axis=1,inplace=True)

#Drop original output columns
remove=["wrong_bands","right_band","no_bands"]
PCR_Data.drop(remove,axis=1,inplace=True)

#Data Cleaning

#rename incorrect column labels to correct column labels
PCR_Data.rename(columns={'extension_t2_min.1':'extension_t2_sec','t5_sec':'t3_sec'},inplace=True)

#Convert time inputs to seconds
PCR_Data["t1(sec)"] = (PCR_Data["t1_min"]*60) + PCR_Data["t1_sec"]
PCR_Data["melting_t2(sec)"] = (PCR_Data["melting_t2_min"]*60) + PCR_Data["melting_t2_sec"]
PCR_Data["annealing_t3(sec)"] = (PCR_Data["annealing_t2_min"]*60) + PCR_Data["annealing_t2_sec"]
PCR_Data["extension_t4(sec)"] = (PCR_Data["extension_t2_min"]*60) + PCR_Data["extension_t2_sec"]
PCR_Data["t5(sec)"] = (PCR_Data["t3_min"]*60) + PCR_Data["t3_sec"]

#drop old time columns
remove=["t1_min","t1_sec","melting_t2_min","melting_t2_sec","annealing_t2_min",
        "annealing_t2_sec","extension_t2_min","extension_t2_sec","t3_min","t3_sec"]
PCR_Data.drop(remove,axis=1,inplace=True)

#Feature Engineering

from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
import primer3 as pr
import melting as melt
from NickFunctions import*


#Change polymerase name to generalized form for Buffer dictionary
PCR_Data["polymerase_name"]=PCR_Data['polymerase_name'].apply(lambda x: ChangeName(x))
PCR_Data.head(2)

#Strip blank spaces from primers
PCR_Data["primer1"] = PCR_Data["primer1"].apply(lambda x: x.replace(' ',''))# remove all white space including middle of string
PCR_Data["primer2"] = PCR_Data["primer2"].apply(lambda x: x.replace(' ',''))# remove all white space including middle of string
#Ensure all primers are upper case
PCR_Data["primer1"] = PCR_Data["primer1"].apply(lambda x: x.upper())
PCR_Data["primer2"] = PCR_Data["primer2"].apply(lambda x: x.upper())

#Add Primer GC content
PCR_Data["primer1_GC"] = PCR_Data["primer1"].apply(lambda x: GC(x))
PCR_Data["primer2_GC"] = PCR_Data["primer2"].apply(lambda x: GC(x))

#Add Primer Lengths
PCR_Data["P1_len"] = PCR_Data["primer1"].apply(lambda x: len(x))
PCR_Data["P2_len"] = PCR_Data["primer2"].apply(lambda x: len(x))
#Remove primers with length greater than 60 bioinformatics algorithms do not work well when greater than 60
longPs = PCR_Data[(PCR_Data["P1_len"]>60) | (PCR_Data["P2_len"]>60)].index 
PCR_Data.drop(longPs,inplace=True)

#Add Primer melting temps
#Primer 3 no correction original model
PCR_Data["Tm_p1_p3orig"]=PCR_Data.apply(lambda x: pr.calcTm(x.primer1),axis=1)
PCR_Data["Tm_p2_p3orig"]=PCR_Data.apply(lambda x: pr.calcTm(x.primer2),axis=1)
#Primer3
PCR_Data["Tm_p1_p3"]=PCR_Data.apply(lambda x: Tm(primer=x.primer1,polymerase=x.polymerase_name,primerconc=(x.amtprimer1+x.amtprimer2),option=1),axis=1)
PCR_Data["Tm_p2_p3"]=PCR_Data.apply(lambda x: Tm(primer=x.primer2,polymerase=x.polymerase_name,primerconc=(x.amtprimer1+x.amtprimer2),option=1),axis=1)
#BioseqUtils
PCR_Data["Tm_p1_BSU"]=PCR_Data.apply(lambda x: Tm(primer=x.primer1,polymerase=x.polymerase_name,primerconc=(x.amtprimer1+x.amtprimer2),option=2),axis=1)
PCR_Data["Tm_p2_BSU"]=PCR_Data.apply(lambda x: Tm(primer=x.primer2,polymerase=x.polymerase_name,primerconc=(x.amtprimer1+x.amtprimer2),option=2),axis=1)
#Melting (IDT analog)
PCR_Data["Tm_p1_IDT"]=PCR_Data.apply(lambda x: Tm(primer=x.primer1,polymerase=x.polymerase_name,primerconc=(x.amtprimer1+x.amtprimer2),option=3),axis=1)
PCR_Data["Tm_p2_IDT"]=PCR_Data.apply(lambda x: Tm(primer=x.primer2,polymerase=x.polymerase_name,primerconc=(x.amtprimer1+x.amtprimer2),option=3),axis=1)
#Average
#PCR_Data["Tm_p1_Avg"]=(PCR_Data["Tm_p1_p3"]+PCR_Data["Tm_p1_BSU"]+PCR_Data["Tm_p1_IDT"])/3
#PCR_Data["Tm_p2_Avg"]=(PCR_Data["Tm_p2_p3"]+PCR_Data["Tm_p2_BSU"]+PCR_Data["Tm_p2_IDT"])/3

#difference in primer melting temps
PCR_Data["Tm_dif_IDT"]=abs(PCR_Data["Tm_p1_IDT"]-PCR_Data["Tm_p2_IDT"])

#Add if hairpin structure is predicted
PCR_Data["HP_P1"] = PCR_Data.apply(lambda x: Tm(primer=x.primer1,polymerase=x.polymerase_name,primerconc=(x.amtprimer1+x.amtprimer2),option=4),axis=1)
PCR_Data["HP_P2"] = PCR_Data.apply(lambda x: Tm(primer=x.primer2,polymerase=x.polymerase_name,primerconc=(x.amtprimer1+x.amtprimer2),option=4),axis=1)

#Combine hairpins to one feature that has the highest hairpin temp
#if first statement is true return 2nd value first statement is false return third value.
PCR_Data["HPin_Tm_Max"] = np.where(PCR_Data["HP_P1"]>PCR_Data["HP_P2"],PCR_Data["HP_P1"],PCR_Data["HP_P2"])

#(PCR_Data.size)/len(PCR_Data.columns)

#Calculate Heterodimer
PCR_Data["Heterodimer_tm"] = PCR_Data.apply(lambda x : Tm(primer=x.primer1,primer2=x.primer2,polymerase=x.polymerase_name,primerconc=(x.amtprimer1+x.amtprimer2),option=6),axis=1)

#Calculate Homodimer
PCR_Data["Homodimer_tm1"] = PCR_Data.apply(lambda x : Tm(primer=x.primer1,polymerase=x.polymerase_name,primerconc=(x.amtprimer1+x.amtprimer2),option=5),axis=1)
PCR_Data["Homodimer_tm2"] = PCR_Data.apply(lambda x : Tm(primer=x.primer2,polymerase=x.polymerase_name,primerconc=(x.amtprimer1+x.amtprimer2),option=5),axis=1)
PCR_Data["Homodimer_maxTm"] = np.where(PCR_Data["Homodimer_tm1"]>PCR_Data["Homodimer_tm2"],PCR_Data["Homodimer_tm1"],PCR_Data["Homodimer_tm2"])
PCR_Data["Homodimer_minTm"] = np.where(PCR_Data["Homodimer_tm1"]<PCR_Data["Homodimer_tm2"],PCR_Data["Homodimer_tm1"],PCR_Data["Homodimer_tm2"])

def changeTm(Tm):
    '''
    Tms of less than or equal to zero mean there was no significant structure. To help the RF make these temps all 0
    '''
    if Tm <0:
        return 0
    else:
        return Tm
#Homodimer and Heterodimer calculation
PCR_Data["Homodimer_maxTm"]=PCR_Data["Homodimer_maxTm"].apply(lambda x: changeTm(x))
PCR_Data["Heterodimer_tm"]=PCR_Data["Heterodimer_tm"].apply(lambda x: changeTm(x))

#Separate variables by minimum and maximum
PCR_Data['MaxTm_IDT']=np.where(PCR_Data["Tm_p1_IDT"]>PCR_Data["Tm_p2_IDT"],PCR_Data["Tm_p1_IDT"],PCR_Data["Tm_p2_IDT"])
PCR_Data['MinTm_IDT']=np.where(PCR_Data["Tm_p1_IDT"]<PCR_Data["Tm_p2_IDT"],PCR_Data["Tm_p1_IDT"],PCR_Data["Tm_p2_IDT"])
PCR_Data['MaxGC']=np.where(PCR_Data["primer1_GC"]>PCR_Data["primer2_GC"],PCR_Data["primer1_GC"],PCR_Data["primer2_GC"])
PCR_Data['MinGC']=np.where(PCR_Data["primer1_GC"]<PCR_Data["primer2_GC"],PCR_Data["primer1_GC"],PCR_Data["primer2_GC"])
PCR_Data['MaxPlen']=np.where(PCR_Data["P1_len"]>PCR_Data["P2_len"],PCR_Data["P1_len"],PCR_Data["P2_len"])
PCR_Data['MinPlen']=np.where(PCR_Data["P1_len"]<PCR_Data["P2_len"],PCR_Data["P1_len"],PCR_Data["P2_len"])
PCR_Data['TotPrimer']=PCR_Data['amtprimer1']+PCR_Data['amtprimer2']

#Create GC_Clamp_features Score
PCR_Data["GC_Clamp1"] = PCR_Data["primer1"].apply(lambda x: GC_Clamp_features(x).Score)
PCR_Data["GC_Clamp2"] = PCR_Data["primer2"].apply(lambda x: GC_Clamp_features(x).Score)


# Finalize GC_Clamp features
#Take the lowest GC in clamp and positional weight for p1 and p2
PCR_Data["GC_ClampMin"] = np.where(PCR_Data["GC_Clamp1"]<PCR_Data["GC_Clamp2"],PCR_Data["GC_Clamp1"],PCR_Data["GC_Clamp2"])
PCR_Data["GC_ClampMax"] = np.where(PCR_Data["GC_Clamp1"]>PCR_Data["GC_Clamp2"],PCR_Data["GC_Clamp1"],PCR_Data["GC_Clamp2"])

#Drop GC_CLamp columns used for analysis and primers because we are done pulling data from them
remove=["GC_Clamp1","GC_Clamp2","primer1","primer2","polymerase_name","cycles1","cycles3","HP_P1","HP_P2","Homodimer_tm1","Homodimer_tm2","Tm_p1_IDT","Tm_p2_IDT","primer1_GC","primer2_GC","P1_len","P2_len",'amtprimer1','amtprimer2']
PCR_Data.drop(remove,axis=1,inplace=True)

#Move outcome to the first position
PCR_Data=PCR_Data[['temp1', 'melting_temp2',
       'annealing_temp2', 'extension_temp2', 'cycles2', 'temp3',
       'band_size', 'dna_amt','t1(sec)', 'melting_t2(sec)',
       'annealing_t3(sec)', 'extension_t4(sec)', 't5(sec)', 'Tm_p1_p3orig', 'Tm_p2_p3orig', 'Tm_p1_p3',
       'Tm_p2_p3', 'Tm_p1_BSU', 'Tm_p2_BSU', 'Tm_dif_IDT',
       'HPin_Tm_Max', 'Heterodimer_tm', 'Homodimer_maxTm', 'GC_ClampMin',
       'GC_ClampMax', 'MaxTm_IDT','MinTm_IDT','MaxGC','MinGC','MaxPlen','MinPlen','TotPrimer']]

#print('before dropping duplicates',len(PCR_Data))
#PCR_Data.drop_duplicates(keep='first',inplace=True)
#print('after dropping duplicates',len(PCR_Data))

engfname = outFile+date[0]+".csv"
print("Writing Engineered Features to csv:", engfname)
PCR_Data.to_csv(engfname,index=False)

#PCR Predictor
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from sklearn import preprocessing   
from sklearn.utils import shuffle
from sklearn.ensemble import RandomForestClassifier

print("Predicting Reaction Outcomes...")

PCR_Data_M = pd.read_csv('ProcessedTrainingData.csv')

Test_Data_M = pd.read_csv(engfname)

initialRemove=['outcome','outcome_clean','Tm_p1_p3orig', 'Tm_p2_p3orig', 'Tm_p1_p3','Tm_p2_p3','Tm_p1_BSU','Tm_p2_BSU','temp1','MaxPlen','MinPlen']#features to remove from trainingdata
initialRemove2=['Tm_p1_p3orig', 'Tm_p2_p3orig', 'Tm_p1_p3', 'Tm_p2_p3','Tm_p1_BSU', 'Tm_p2_BSU','temp1','MaxPlen','MinPlen']#features to remove from TestData
PCR_Data=PCR_Data_M.drop(initialRemove,axis=1)
Test_Data=Test_Data_M.drop(initialRemove2,axis=1)

def name(df,df_norm):
    names={}
    for i in range(len(df.columns.values)):
        names[i]=df.columns.values[i]
    df_norm.rename(columns=names,inplace=True)
    return df_norm
#Normalize PCR_Data and normalize test PCR_Data
#normalize columns
feat = PCR_Data.values #returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()#create scaler object
feat_scaled = min_max_scaler.fit_transform(feat)#scale each feature between 0 and 1
PCR_Data_norm = pd.DataFrame(feat_scaled)#create PCR_Dataframe from normalized PCR_Data
PCR_Data_norm = name(PCR_Data,PCR_Data_norm)

#Fit predict data
feat_pred = Test_Data.values
pred_scaled = min_max_scaler.transform(feat_pred)
Test_norm = pd.DataFrame(pred_scaled)
Test_norm = name(Test_Data,Test_norm)

finalremove=['outcome_clean','HPin_Tm_Max','t1(sec)','annealing_t3(sec)','MinGC','Heterodimer_tm','Tm_dif_IDT','MaxGC','extension_t4(sec)']#remove extreaneous features
rfc = RandomForestClassifier(n_estimators=1000,oob_score=True,class_weight="balanced",max_features='auto',n_jobs=-1,random_state=42)#initialize rfc 
PCR_Data_norm['outcome_clean']=PCR_Data_M['outcome_clean'].copy(deep=True)#add outcome column back
PCR_Data_norm_t = shuffle(PCR_Data_norm,random_state=42)#shuffle data
X = PCR_Data_norm_t.drop(finalremove,axis=1) #(axis zero index 1 column)
y = PCR_Data_norm_t['outcome_clean']

rfc.fit(X,y)#fit data

finalremove_p=['HPin_Tm_Max','t1(sec)','annealing_t3(sec)','MinGC','Heterodimer_tm','Tm_dif_IDT','MaxGC','extension_t4(sec)']

X_pred=Test_norm.drop(finalremove_p,axis=1)
outcome=rfc.predict(X_pred)

print("Reaction Outcomes Predicted: ",outcome)
print("Exporting Predictions: ", Predictions)
#Export outcomes to excel file with important features
lbs1=pd.read_csv('FeatureRelabel.csv')
lbs1[lbs1['Old Title']=='Hpin_Tm_Max']='HPin_Tm_Max'
lbs=dict(zip(lbs1['Old Title'],lbs1['New Title']))
s=['t1(sec)', 'annealing_t3(sec)', 'MinGC', 'Heterodimer_tm', 'Tm_dif_IDT', 'MaxGC', 'extension_t4(sec)']
mapper={}
for l in lbs.keys():
  if l not in s:
      mapper[l]=lbs[l]
Test_Data.rename(mapper=mapper,axis=1,inplace=True)
Test_Data.drop(finalremove_p,axis=1,inplace=True)
Test_Data['MODEL_PREDICTION']=outcome
#write file
excelwriter =pd.ExcelWriter(Predictions)
Test_Data[['MODEL_PREDICTION','Melting Cycle Temp. (C)', 'Annealing Cycle Temp (C)',
       'Extension  Temp (C)', 'Amplification Cycles (cycles)',
       'Final Extension Cycle Temp. (C)', 'Product Band Size (bp)',
       'Template Amount (ng)', 'Melting Cycle Time (s)',
       'Final Extension Time (s)', 'Max Homodimer Tm (C)',
       'Min GC Clamp Strength (score)', 'Max GC Clamp Strength (score)',
       'Max Primer Tm (C)', 'Min Primer Tm (C)', 'Total Primer (uM)']].to_excel(excelwriter,sheet_name='Outcomes',index=False)
for c in Test_Data: #fix column width
  width = max(Test_Data[c].astype(str).map(len).max(), len(c))
  index = Test_Data.columns.get_loc(c)
  excelwriter.sheets['Outcomes'].set_column(index,index,width)
excelwriter.save()