#NickFunctions.py>
from Bio.SeqUtils import MeltingTemp as mt
import primer3 as pr
import melting as melt
def Buffer(polymerase):
    '''
    Description:
    ------------
    Takes ChangeName simiplified polymerase name and returns salt concentrations in mM for reaction conditions for that polymerase.
    
    Parameters:
    ------------
    polymerase, str: simplified polymerase name
    
    Return:
    ------------
    list of salt types with concentrations in mM for primer feature functions.
    
    '''
    Salts={'Pfu':{"Tris":20.0,"Na":0.0,"K":10.0,"Mg":1.5,"dNTPs":0.8, "MV_tot":20.0,"DV_tot":1.5},
           'LongAmp':{"Tris":60.0,"Na":0.0,"K":0.0,"Mg":2.0,"dNTPs":1.2, "MV_tot":20.0,"DV_tot":2.0},
           'PhusionHF':{"Tris":25.0,"Na":0.0,"K":50.0,"Mg":1.5,"dNTPs":0.8, "MV_tot":50,"DV_tot":1.5},
           'SuperMix':{"Tris":66.0,"Na":0.0,"K":0.0,"Mg":2.4,"dNTPs":0.88, "MV_tot":19.8,"DV_tot":2.4},
           'Vent':{"Tris":20.0,"Na":0.0,"K":10.0,"Mg":2.0,"dNTPs":0.8, "MV_tot":20.0,"DV_tot":2.0},
           'PrimeSTAR':{"Tris":0,"Na":0.0,"K":0,"Mg":1,"dNTPs":1.0, "MV_tot":50,"DV_tot":1.0},
           'Custom':{"Tris":10,"Na":0.0,"K":50.0,"Mg":1.5,"dNTPs":0.8, "MV_tot":50,"DV_tot":1.5},
           'Taq':{"Tris":0.0,"Na":0.0,"K":50.5,"Mg":1.5,"dNTPs":0.8, "MV_tot":50.5,"DV_tot":1.5},
           'Q5':{"Tris":25.0,"Na":0.0,"K":50.0,"Mg":2.0,"dNTPs":0.2, "MV_tot":50.0,"DV_tot":2.0}
            }
    return Salts[polymerase]
def ChangeName(polymerase):
    '''
    Description:
    ------------
    Crude function that takes polymerase name from common user input and changes to simplified polymerase name for key use in "Buffer"  function. Can also be written using polymerase catalog number, but was easier to manually lookup and confirm concentrations by name.
     
    Parameters:
    ------------
    polymerase, str: original polymerase name from user input. Users typically input the same name, but this function will have to be manually modified based on new polymerase values. Can be simplified with the polymerase catalog number.
    
    Return:
    ------------
    Simplified polymerase name as a string.
    
   
    '''
    Taq=['Taq','Taq DNA Polymerase']
    for Pol in Taq:
        if Pol==polymerase:
            return "Taq"
    Custom=['Custom','Custom Pfu']
    for Pol in Custom:
        if Pol==polymerase:
            return  "Custom"
    SuperMix=['SuperMix','PCR SuperMix High Fidelity','PCR SuperMix High Fidelity - Invitrogen/10790020','Platinum PCR SuperMix High Fidelity','Platinum SuperMix High Fidelity','Platinum PCR SuperMix High Fidelity - Invitrogen/12532016']
    for Pol in SuperMix:
        if Pol==polymerase:
            return "SuperMix"
    # NEED TO REMOVE LINES WITH QUICKCHANGE LIGHTNING TEMPORARY FIX
    Pfu=['Pfu','QuikChange Lightning Enzyme','Pfu Turbo','PfuTurbo DNA Polymerase - Agilent/600250 ','PfuUltra High-Fidelity DNA Polymerase','PFU Ultra HF polymerase','PfuTurbo Cx Hotstart DNA Polymerase']
    for Pol in Pfu:
        if Pol==polymerase:
            return "Pfu"
    PhusionHF=['PhusionHF','Phusion High-Fidelity PCR Master Mix with HF Buffer','Phusion High-Fidelity PCR Master Mix with HF Buffer - NEB/M0531S','Phusion DNA polymerase','Phusion Hot Start II High-Fidelity DNA Polymerase','Phusion HF master mix']
    for Pol in PhusionHF:
        if Pol==polymerase:
            return "PhusionHF"
    Vent=['Vent','Vent DNA Polymerase']
    for Pol in Vent:
        if Pol==polymerase:
            return "Vent"
    PrimeSTAR=['PrimeSTAR','PrimeSTAR HS DNA Polymerase - Takara/R010A']
    for Pol in PrimeSTAR:
        if Pol==polymerase:
            return "PrimeSTAR"
    LongAmp=['LongAmp','LongAmp Taq DNA Polymerase - NEB/M0323S','LongAmp Taq DNA Polymerase']
    for Pol in LongAmp:
        if Pol==polymerase:
            return "LongAmp"
    #if has not been found
    Q5=['Q5','Q5 High-Fidelity DNA Polymerase']
    for Pol in Q5:
        if Pol==polymerase:
            return 'Q5'
    return None
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

class GC_Clamp_features(object):
    '''
    Description:
    ------------
    Class obtains features from GC clamp such as strength and num_GC 
    NN model strength table taken from 
    Khandelwal G, Bhyravabhotla J. A Phenomenological Model for Predicting Melting Temperatures of DNA Sequences. PLOS ONE. 2010;5: e12433. doi:10.1371/journal.pone.0012433
    GC Clamp is defined by a G or C in the last 5 nucleotides on the 3' end of the primer.
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