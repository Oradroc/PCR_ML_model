#Adjust Data
#python Training_PCR_Feature_Engineering.py TrainingDataRaw.csv TrainingDataProcessed.csv
#import libraries
import pandas as pd
import numpy as np
import os
import sys

#import data

combinedF = sys.argv[1]
outFile = sys.argv[2]
PCR_Data = pd.read_csv(combinedF)
#Data Setup
#Drop duplicate data 10/26/20
print('before dropping duplicates',len(PCR_Data))
PCR_Data.drop_duplicates(keep='first',inplace=True)
print('after dropping duplicates',len(PCR_Data))
#Drop unnecessary columns
remove=["timestamp","user","email","species","polymerase_catalog_number","pcrid",
        "template_type","polymerase_brand","dna_amt_units","band_size_units",'customized_salts']
PCR_Data.drop(remove,axis=1,inplace=True)
PCR_Data.dropna(inplace=True)
#print('before dropping duplicates',len(PCR_Data))
#PCR_Data.drop_duplicates(keep='first',inplace=True)
#print('after dropping duplicates',len(PCR_Data))
#Set Outcomes

def outcome_bool(a):
    '''
    Description:
    ------------
    Converts Truth values for inputs right_band, wrong_band, no_bands into 1 or 0 for success or fail respectively.
    Classifies "dirty" PCRs as success by only considering the "right_band" input.
    
    Parameters:
    ------------
    right_band; bool: True or False
    
    Return:
    ------------
    int: 1 or 0
    '''
    if a:
        return 1
    else:
        return 0
def outcome_bool_clean(a,b,c):
    '''
    Description:
    ------------
    Converts Truth values for inputs right_band, wrong_band, no_bands into 1 or 0 for success or fail.
    Classifies only "clean" PCRs as successes by considering both right_band and wrong_band inputs.
    
    Parameters:
    ------------
    a,b,c ; bool: True or False #right_band, wrong_bands, no_bands
    
    Return:
    ------------
    int: 1 or 0
    '''
    if c: #no bands
        return 0
    elif (a and not b): #right bands only
        return 1
    else: #wrong bands and right bands
        return 0

#outcome
PCR_Data["outcome"] = PCR_Data["right_band"].apply(lambda x:outcome_bool(x))

#Only 100%successful no nonspecific products
PCR_Data['outcome_clean'] = PCR_Data.apply(lambda x: outcome_bool_clean(x.right_band,x.wrong_bands,x.no_bands),axis=1)
PCR_Data[(PCR_Data['right_band']==True)&(PCR_Data['wrong_bands']==True)].head(3)

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

def Tm(primer,polymerase,primerconc, option,primer2='NA'):
    '''
    
    
    Description:
    ------------
    Calculates primer melting temperature (Tm) based on a variety of options. Uilizes Buffer function to lookup salt concentrations for each function. Multiple options for testing.
     NOTE: primer2 only necessary when calculating heterodimer
    Options (default 3): 
    1 - Primer3 melting temperature
    2 - Tm_NN BioSeqUtils
    3 - melting (Virtually the same as IDT)
    4 - Hairpin
    5 - Homodimer
    6 - Heterodimer need to input primer2
    
    Differences:

    IDT unsure what NN method is used for calculation. Should use Owczarzy (2008) salt correction.

    Primer3 uses santaLucia (1997) correction. Uses Owczarzy (2008) salt correction

    BioSeqUtils uses santaLucia (2004) review correction. Uses Owczarzy (2008) salt correction
    
    Parameters:
    ------------
    primer, str: primer to calculate melting temperature for
    polymerase, str: simplified polyerase name from "ChangeName" function. Used for lookup table of salt concentrations.
    primerconc, float: from user input
    option, int: melting temp prediction option
    primer2, str: only needed for heterodimer calculation option.
    
    Return:
    ------------
    calculated melting temperature as float.
    '''
    
    b=Buffer(polymerase)#salts in each polymerase buffer
    Tms=[]
    if (option==1):
        return pr.calcTm(seq=primer,
                                    mv_conc=b['MV_tot'],
                                    dv_conc=b['DV_tot'],
                                    dntp_conc=b['dNTPs'],
                                    dna_conc=primerconc,
                                    max_nn_length=60,
                                    tm_method='santalucia',
                                    salt_corrections_method=2)
    elif (option==2):
        return mt.Tm_NN(seq=primer,
                                    nn_table=mt.DNA_NN4,
                                    dnac1=primerconc,
                                    dnac2=0,
                                    Na=b['Na'], 
                                    K=b['K'],
                                    Tris=b['Tris'],
                                    Mg=b['Mg'],
                                    dNTPs=b['dNTPs'],
                                    saltcorr=7)
    elif (option==3):
        return melt.temp(primer,
                                    DNA_c=primerconc,
                                    Na_c=b['MV_tot'],
                                    Mg_c=b['DV_tot'],
                                    dNTPs_c=b['dNTPs'])
    elif (option==4):
        return pr.calcHairpin(seq=primer,
                                    mv_conc=b['MV_tot'],
                                    dv_conc=b['DV_tot'],
                                    dntp_conc=b['dNTPs'],
                                    dna_conc=primerconc,
                                    temp_c=37,
                                    max_loop=30,
                                    output_structure=False).tm
    elif (option==5):
        return pr.calcHomodimer(seq=primer,
                                    mv_conc=b['MV_tot'],
                                    dv_conc=b['DV_tot'],
                                    dntp_conc=b['dNTPs'],
                                    dna_conc=primerconc,
                                    temp_c=37,
                                    max_loop=30,
                                    output_structure=False).tm
    elif (option==6):
        return pr.calcHeterodimer(seq1=primer, seq2=primer2,
                                    mv_conc=b['MV_tot'],
                                    dv_conc=b['DV_tot'],
                                    dntp_conc=b['dNTPs'],
                                    dna_conc=primerconc,
                                    temp_c=37,
                                    max_loop=30,
                                    output_structure=False).tm
    else:
        print("Not a valid option. Please input 1,2,3,4,5, or 6")
        return 0

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

# **GC Clamp:**
# 
# According to Top Tip Bio (TTP) (https://toptipbio.com/gc-clamp-pcr/) GC Clamp is defined by a G or C in the last 5 nucleotides on the 3' end of the primer.
# 
# Strength table for NN model taken from Khandelwahl (2010)

class GC_Clamp_features(object):
    '''
    Description:
    ------------
    Class obtains features from GC clamp such as strength and num_GC 
    NN model strength table taken from 
    Khandelwal G, Bhyravabhotla J. A Phenomenological Model for Predicting Melting Temperatures of DNA Sequences. PLOS ONE. 2010;5: e12433. doi:10.1371/journal.pone.0012433
    '''
    #initialize class variables
    clamp = ""
    num_GC = 0
    Strength={'GC': 13,'CC': 11,'GG': 11,'CG': 10,'AC': 10,'TC': 8,'AG': 8,'TG': 7,'GT': 10,'CT': 8,'GA': 8,'CA': 7,'AT': 7,'TT': 5,'AA': 5,'TA': 4}
    Score=0
    #GC_Clamp class constructor
    def __init__(self,primer):
        self.clamp = primer[-5:]
        for nt in range(len(self.clamp)-1):#inrement by 1 nucleotide at a time
            NN=self.clamp[nt]+self.clamp[nt+1]#grab dinucleotides
            self.Score+=self.Strength[NN]#GC clamp strength score
        for nt in range(len(self.clamp)):
            if self.clamp[nt]=="G" or self.clamp[nt]=="C":
                self.num_GC+=1

'''
#Test GC_Clamp_features class
clamp = GC_Clamp_features("ATACGTACGATAGCA")
print(clamp.clamp)
print(clamp.num_GC)
print(clamp.Score)'''

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
PCR_Data=PCR_Data[['outcome', 'outcome_clean','temp1', 'melting_temp2',
       'annealing_temp2', 'extension_temp2', 'cycles2', 'temp3',
       'band_size', 'dna_amt','t1(sec)', 'melting_t2(sec)',
       'annealing_t3(sec)', 'extension_t4(sec)', 't5(sec)', 'Tm_p1_p3orig', 'Tm_p2_p3orig', 'Tm_p1_p3',
       'Tm_p2_p3', 'Tm_p1_BSU', 'Tm_p2_BSU', 'Tm_dif_IDT',
       'HPin_Tm_Max', 'Heterodimer_tm', 'Homodimer_maxTm', 'GC_ClampMin',
       'GC_ClampMax', 'MaxTm_IDT','MinTm_IDT','MaxGC','MinGC','MaxPlen','MinPlen','TotPrimer']]
PCR_Data[(PCR_Data['outcome']==1)&(PCR_Data['outcome_clean']==0)]


#print('before dropping duplicates',len(PCR_Data))
#PCR_Data.drop_duplicates(keep='first',inplace=True)
#print('after dropping duplicates',len(PCR_Data))

import datetime
now = datetime.datetime.now()
now = str(now)
date = now.split(" ")

engfname = outFile+date[0]+".csv"
print(engfname)

PCR_Data.to_csv(engfname,index=False)
