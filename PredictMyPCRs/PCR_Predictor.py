
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from sklearn import preprocessing   
from sklearn.utils import shuffle
from sklearn.ensemble import RandomForestClassifier

#UserInput
TestFile = sys.argv[1]
outFile = sys.argv[2]+'.xlsx'

PCR_Data_M = pd.read_csv('ProcessedTrainingData.csv')

Test_Data_M = pd.read_csv(TestFile)

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

print(outcome)

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
excelwriter =pd.ExcelWriter(outFile)
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
       
