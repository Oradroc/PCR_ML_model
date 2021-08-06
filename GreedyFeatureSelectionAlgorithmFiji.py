#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd
import numpy as np
import os
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import shuffle
from sklearn.base import clone
from sklearn.model_selection import cross_val_score
from sklearn import preprocessing
from scipy import stats
from scipy.stats import sem
import copy
from joblib import Parallel, delayed
import pickle


# In[4]:


def ZScore(data_nor,drop,axis):
    '''
    Description:
    ------------
    zscore normalizes rows of dataframe 
    
    Parameters:
    ------------
    data_nor, DataFrame: dataframe with data to be normalized rowise
    drop, list: list of columns to be left out of the normalization
    axis, int: 1 for columns 0 for rows
    Return:
    ------------
    zscore normalized dataframe
    '''
    #calculate zscore for each column
    if axis==1:#normalize column
        features=set(data_nor.columns.values)-set(drop)
        for feat in features:
            mew=data_nor[feat].mean()
            std=data_nor[feat].std(ddof=1)
            data_nor[feat]=data_nor[feat].apply(lambda x: ((x-mew)/std))
        return data_nor
    else: #normalize row axis==0
        vector=data_nor.drop(drop,axis=1).values #get data frame rows as array of numpyy arrays
        columns=data_nor.drop(drop,axis=1).columns
        scaledrows=[]
        vlen=len(vector)
        for i in range(vlen):
            mew=vector[i].mean()
            std=vector[i].std(ddof=1)
            srow=[]
            row=vector[i]
            lrow=len(row)
            for j in range(lrow):
                srow.append((row[j]-mew)/std)#append zscore
            scaledrows.append(np.array(srow))#append scaled row
        scaledrows=np.array(scaledrows)
        #print(len(scaledrows),len(columns))
        ZSCORED=pd.DataFrame(data=scaledrows,columns=columns)
        ZSCORED['outcome_clean']=data_nor['outcome_clean']
        return ZSCORED
def NormalizeDf(data):
    '''
    MinMaxNormalization
    '''
    #normalize columns
    feat = data.values #returns a numpy array
    min_max_scaler = preprocessing.MinMaxScaler()#create scaler object
    feat_scaled = min_max_scaler.fit_transform(feat)#scale each feature between 0 and 1
    data_norm = pd.DataFrame(feat_scaled)#create dataframe from normalized data
    names={}
    for i in range(len(data.columns.values)):
        names[i]=data.columns.values[i]
    data_norm.rename(columns=names,inplace=True)
    data_norm['random']=np.random.random(size=len(data))#add random feature
    return data_norm
    
def Best_or_Worst_Score(scoreA,scoreF,ScoreType,best=True):
    '''
    Description:
    ------------
    Determines best new feature for model from dictionary data
    
    Parameters:
    ------------
    scoreA, dict: dictionary with feature as key and numpy array of scores
    scoreF, dict: dictionary with feature as key and numpy array of scores
    ScoreType, string: determines if 'best feature' is decided with accuracy (scoreA) or f1 (scoreF)
    best, bool: If evaluation should be done by best value (True) or worst value (False)
    Return:
    ------------
    Function returns feature addition that lead to greatest improvement in score
    '''
    if best:
        bestF=''
        bestV=0
        bestValues='None'
        if ScoreType=='accuracy':
            k=scoreA.keys()
            for f in k:
                S=np.array(scoreA[f]).flatten().mean()
                if S>bestV:
                    bestF=f
                    bestValues=np.array(scoreA[f]).flatten()
                    bestV=bestValues.mean()
                    
        else:
            k=scoreF.keys()
            for f in k:
                S=np.array(scoreF[f]).flatten().mean()
                if S>bestV:
                    bestF=f
                    bestValues=np.array(scoreA[f]).flatten()
                    bestV=bestValues.mean()
        return bestF,bestV
    else:
        worstF=''
        worstV=np.inf
        worstValues='None'
        if ScoreType=='accuracy':
            k=scoreA.keys()
            for f in k:
                S=np.array(scoreA[f]).flatten().mean()
                if S<worstV:
                    worstF=f
                    worstValues=np.array(scoreA[f]).flatten()
                    worstV=bestValues.mean()
        else:
            k=scoreF.keys()
            for f in k:
                S=np.array(scoreF[f]).flatten().mean()
                if S<worstV:
                    worstF=f
                    worstValues=np.array(scoreA[f]).flatten()
                    worstV=bestValues.mean()
                    
        return worstF,worstV
def Forests(data_norm,feature,featuresAdded,rfc,runs,cv):
    '''
    Description:
    ------------
    Parallelized model building. Normalizes rows (zscore) when 3 or more features are present.
    
    Parameters:
    ------------
    data, DataFrame: data for processing
    featuresAdded, list: features already in the model
    rfc, object: random forest classifier object 
    runs, int: number of cross validation iterations to run on the data
    cv, int: number of folds for cross validation
    
    Return:
    ------------
    Tuple: (accuracy,f1) ({feature:np.array(accuracy)},{feature:np.array(f1)})
    '''
    data_nor=copy.deepcopy(data_norm)#create deep copy
    tempFeats=copy.deepcopy(featuresAdded)
    tempFeats.append(feature)
    drop=(set(data_nor.columns.values)-set(tempFeats))#drop all features except temp feats
    numfeats=len(featuresAdded)
    #Initialize score savers
    scoreA={feature:[]}#{feature:np.array(scores)} storing reference scores
    scoreF={feature:[]}#{feature:np.array(scores)} storing reference scores
    #Run runs times for this model add to scores to dict for best
    for i in range(runs):
        #reshuffle folds 6times
        data_nor = shuffle(data_nor,random_state=i) #Use different random state shuffle for each cv calc
        X = data_nor.drop(drop,axis=1) #(axis zero index 1 column)
        y = data_nor['outcome_clean']
        #record 4 fold cross validation
        accuracy = cross_val_score(rfc, X, y, cv=cv, scoring='accuracy')
        f1 = cross_val_score(rfc, X, y, cv=cv, scoring='f1')
        scoreA[feature].append(accuracy)
        scoreF[feature].append(f1)
    scoreA[feature]=np.array(scoreA[feature]).flatten()
    scoreF[feature]=np.array(scoreF[feature]).flatten()
    return scoreA,scoreF
                
def FwdFeatureImportance(data_nor,initialRemove,numfeats,cores,ScoreType='accuracy',cv=10,runs=5):
    '''
    Description:
    ------------
    Adds features to model based on best accuracy or f1 score. 
    Records reference score and score not used as a reference.
    
    Parameters:
    ------------
    data_nor, DataFrame: Normalized data for model
    numfeats, int: Number of final features for the test
    initialRemove, list: List of features to be excluded from the test
    ScoreType, string: Determines 'best' feature to add based on model 'accuracy' or 'f1' from cross validation
    cv, int: number of crossvalidation folds (default=10)
    
    Return:
    ------------
    Function returns tuple containing the following:
    [0] Dictionary containing feature name, score type, and value {'Features Dropped (count)':[],'Score':[],'Value':[]}
    [1] List of features added into the model
    '''
    num=0 
    #create nested dictionary
    modelTracker={'Accuracy':{},'F1':{}}
    remove=list(initialRemove)
    data_nor = shuffle(data_nor,random_state=42)#shuffle data to remove bias in intial structure
    featuresAdded=[]#list of features added so they can be removed from the list
    #intialize instance of rfc
    MostFeatures=set(data_nor.columns.values)-set(remove)#features without the initial remove
    rfc = RandomForestClassifier(n_estimators=1000,oob_score=True,class_weight="balanced",max_features='auto',n_jobs=1)#,max_leaf_nodes=3
    while(len(featuresAdded)<numfeats): #while the number of features added back into the model is less than the number of features
        featuresRemaining=list(MostFeatures-set(featuresAdded))
        #Parallel process optimal model building list of tuples of nested dictionaries
        AllScores=Parallel(n_jobs=cores)(delayed(Forests)(data_nor,feature,featuresAdded,rfc,runs,cv) for feature in featuresRemaining)
        #process parallelized data for finding best score
        scoreA={}
        scoreF={}
        for run in AllScores:
            A=run[0]#Grab Accuaracy
            F=run[1]#Grab Accuaracy
            scoreA[list(A.keys())[0]]=np.array(list(A.values())[0])#feature:list of values
            scoreF[list(F.keys())[0]]=np.array(list(F.values())[0])
        #Determine which model had the best accuracy
        if ScoreType=='accuracy':
            bestFeat,bestScore=Best_or_Worst_Score(scoreA,scoreF,'accuracy')
        elif ScoreType=='F1':
            bestFeat,bestScore=Best_or_Worst_Score(scoreA,scoreF,'F1')
        #Add to model builder
        modelTracker['Accuracy'][bestFeat]=scoreA[bestFeat]
        modelTracker['F1'][bestFeat]=scoreF[bestFeat]
        featuresAdded.append(bestFeat)
    return modelTracker,featuresAdded


# In[5]:


def main():
    #Variables to change
    initialRemove=['outcome','outcome_clean','Tm_p1_p3orig', 'Tm_p2_p3orig', 'Tm_p1_p3', 'Tm_p2_p3','Tm_p1_BSU', 'Tm_p2_BSU','random','temp1','MaxPlen','MinPlen']
    numfeats=23 #number of features
    cores=32 #number of cores to use
    cv=10 #cross validation
    runs=6 #number of tests done for each cv
    exportname='PCR_BestModel23Features_outcomeclean.txt'
    
    #Load Data
    cleanedF = os.path.join("/Users/nico3390/AC_proj/", "Combined_2017-19-20_RawData2020-11-01_myeng110120_outcome_clean.csv")
    PCR=pd.read_csv(cleanedF)
    #drop any rows with na (should already be good)
    PCR.dropna(axis=0,inplace=True)
    PCR_norm=NormalizeDf(PCR) #normalize data
    d,feats=FwdFeatureImportance(PCR_norm,initialRemove,numfeats,cores,ScoreType='accuracy',cv=cv,runs=runs)
    
    #Export nested Dictionary with pickle
    file = open(exportname, "wb") 
    pickle.dump(d,file)
    file.close()
    


# In[6]:


if __name__ == "__main__":
    main()


# In[8]:


'''with open('BestModel2Features.txt', 'rb') as handle:
    b = pickle.load(handle)
'''


# In[ ]:




