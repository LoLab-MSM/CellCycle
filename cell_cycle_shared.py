"""An implementation of the model from:

Mathematical modeling of cell cycle regulation in response to DNA Damage: Exploring mechanisms of 
cell-fate determination. Kuzunari Iwamoto, Hiroyuki Hamada, Yukihiro Eguchi, Masahiro Okamoto. 
Biosystems. 2011 March, 103(2011), pp. 384-391
DOI: 10.1016/j.biosystems.2010.11.011

http://www.sciencedirect.com/science/article/pii/S0303264710002133

Implemented by: Corey Hayford
"""

from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify

# Model()

def set_dna_damage(damage):
    model.parameters['DDS_0'].value = damage

def declare_monomers():
    """Declare the monomers in the Iwamoto model
    'state' is the activity state between the monomers
    'a' is the active conformation binding state
    'i' is the inactive conformation binding state
    'phos' is the phosphorylation state
    'u' is the unphosphorylated monomer
    'p' is the phosphorylated monomer
    'b' is the binding site between p## and other proteins"""
    
    ## G2/M Monomers
    
    Monomer("x14_3_3", ['b'])

#     Monomer("p53")
#     Monomer("p21", ['b'])
#     Monomer("G2_M_Mdm2")
#     Monomer('G2_M_Cdc25', ['b', 'state', 'phos'], {'state':['i','a'], 'phos':['u','p']})
#     Monomer("CycB",['c'])    
#     Monomer("CDK1",['phos','b','c'], {'phos':['u','p']})
    
    ## G1/S Monomers
    
    Monomer("CycD", [        'c'])
    Monomer("CycE", [        'c'])
    Monomer("CycA", [        'c'])
    Monomer("CDK4_6",['phos','b','c'], {'phos':['u','p']})
    Monomer("CDK2",    ['phos','b','c'], {'phos':['u','p']})
    Monomer("p27",     [    'b'])
    Monomer("p16",     [    'b'])
    Monomer("Rb",     ['phos','b'], {'phos':['u','pp','pppp']})
    Monomer("E2F",     [    'b'])
    Monomer("X")
#     Monomer("I")
    
    ## New Monomers
    
    Monomer('NF_Y')
    Monomer('B_Myb',['state'], {'state': ['i', 'a']})
    Monomer('APC', ['state', 'b'], {'state': ['i', 'a']})
    Monomer('Ccdc20', ['b'])
    Monomer('Ccdh1', ['b'])
    
    # **Shared Monomers**

    Monomer("Signal")       #State can be on or off
    Monomer("SignalDamp")   #Dampens signal in Deg(t) function    
    Monomer('CycB', ['c'])
    Monomer('Cdc25', ['b', 'state','state1', 'phos'], {'state':['i','a'], 'state1':['A','C'], 'phos':['u','p']})
    Monomer('CDK1_cyto',['phos','b','c'], {'phos':['u','p']})
    Monomer('CDK1_nuc',['phos','b','c'], {'phos':['u','p']})
    Monomer("p53",    [    'b'])
    Monomer("p21",     [    'b']) 
    Monomer("ATM_ATR", ['b'])      
    Monomer("Mdm2", ['b'])
    Monomer('Wee1', ['phos'], {'phos':['u','p']})
    Monomer('Chk1', ['phos'], {'phos':['u','p']}) #u=inactive p=active
    Monomer("I")
#     Monomer('Im')
    
def declare_parameters():
    
# ***Declare Initial Conditions***
    
    Parameter("shared_Y5_0", 1.50)         #CDK2
    Parameter("shared_Y6_0", 2.0)          #CycD/CDK4_6
    Parameter("shared_Y11_0", 1.4)         #p27
    Parameter("shared_Y25_0", 2.65e-2)     #p53
    Parameter("shared_Y27_0", 0)           #ATM/ATR
    Parameter("shared_Y28_0", 1.0e-3)      #iCdc25A
    Parameter("shared_Y29_0", 1.0e-4)      #aCdc25A
    Parameter("shared_Y30_0", 0.99)        #iChk1
    Parameter("shared_Y31_0", 1.0e-2)      #aChk1
    Parameter("shared_Y32_0", 0)           #NF_Y
    Parameter("shared_Y33_0", 0)           #CycB
    Parameter("shared_Y34_0", 1.0)         #CDK1
    Parameter("shared_Y35_0", 1.0e-4)      #iCycB/CDK1_cyto
    Parameter("shared_Y36_0", 1.0e-4)      #aCycB/CDK1_cyto
    Parameter("shared_Y40_0", 0)           #iB-Myb
    Parameter("shared_Y41_0", 0)           #aB-Myb
    Parameter("shared_Y42_0", 1.0e-6)      #iCdc25C
    Parameter("shared_Y43_0", 1.0e-6)      #aCdc25C
    Parameter("shared_Y44_0", 3.0e-2)      #iCdc25CPs216
    Parameter("shared_Y45_0", 0)           #aCdc25CPs216
    Parameter("shared_Y48_0", 0.9)         #iAPC/Ccdc20
    Parameter("shared_Y49_0", 1.0e-1)      #aAPC/Ccdc20
    Parameter("shared_Y50_0", 1.0e-1)      #iAPC/Ccdh1
    Parameter("shared_Y51_0", 0.9)         #aAPC/Ccdh1
    Parameter("shared_Y52_0", 0)           #iCycB/CDK1_nuc
    Parameter("shared_Y53_0", 0)           #iCycB/CDK1_nuc
    
# *** Declare Kinetic Parameters ***

    Parameter("shared_k1", 5.00e-4)
    Parameter("shared_k5", 1.00e-1)
    Parameter("shared_k7", 2.50e-3)
    Parameter("shared_k8", 2.50e-5)
    Parameter("shared_k9", 3.00e-4)
    Parameter("shared_k11", 5.00e-4)
    Parameter("shared_k14", 7.50e-3)
    Parameter("shared_k15", 5.00e-3)
    Parameter("shared_k16", 5.00e-3)
    Parameter("shared_k17", 5.00e-2)
    Parameter("shared_k22", 6.00e-3)
    Parameter("shared_k26", 2.25e-2)
    Parameter("shared_k28", 9.00e-4)
    Parameter("shared_k29", 5.00e-5)
    Parameter("shared_k35", 5.00e-2)
    Parameter("shared_k38", 1.00e-3)
    Parameter("shared_k61", 7.00e-2)
    Parameter("shared_k81", 1.00e-3)
    Parameter("shared_k84", 1.00e-3)
    Parameter("shared_k85", 5.00e-3)
    Parameter("shared_k86", 5.00e-4)
    Parameter("shared_k89", 1.00e-3)
    Parameter("shared_k90", 5.00e-4)
    Parameter("shared_k91", 2.00e-2)
    Parameter("shared_k92", 5.00e-3)
    Parameter("shared_k93", 1.25e-3)
    Parameter("shared_k94", 2.50e-4)
    Parameter("shared_k95", 5.00e-2)
    Parameter("shared_k96", 1.00e-4)
    Parameter("shared_k97", 5.00e-3)
    Parameter("shared_k98", 5.00e-3)
    Parameter("shared_k103", 2.25e-2)
    Parameter("shared_k104", 1.75e-4)
    Parameter("shared_k105", 5.00e-2)
    Parameter("shared_k106", 5.00e-2)
    Parameter("shared_k107", 2.00e-3)
    Parameter("shared_k109", 1.00e-2)
    Parameter("shared_k111", 1.00e-3)
    Parameter("shared_k113", 1.00e-3)
    Parameter("shared_k114", 1.00e-4)
    Parameter("shared_k116", 1.00e0)
    Parameter("shared_k122", 5.00e-3)
    Parameter("shared_k123", 1.00e-2)
    Parameter("shared_k124", 1.00e-2)
    Parameter("shared_k125", 5.00e-3)
    Parameter("shared_k128", 1.00e-3)
    Parameter("shared_k129", 3.00e-1)
    Parameter("shared_k131", 1.00e-2)
    Parameter("shared_k132", 5.00e-5)
    Parameter("shared_k133", 5.00e-4)
    Parameter("shared_k134", 1.00e-2)
    Parameter("shared_k135", 5.00e-3)
    Parameter("shared_k136", 5.00e-3)
    Parameter("shared_k137", 3.00e-2)
    Parameter("n", 50)
    
    alias_model_components()
    
### ** Initial Conditions **
def declare_initial_conditions():
      
    ## G2/M

    Initial(x14_3_3(b=None), X13_0)                                     #13
    
    ## G1/S
     
    Initial(CycD(c=None), Y0_0)    #0
    Initial(CycE(c=None), Y1_0)    #1
    Initial(CycA(c=None), Y2_0)    #2
    Initial(CDK4_6(phos='u',b=None,c=None), Y3_0)    #3
#     Initial(p27(b=None), Y10_0)    #10
    Initial(p21(b=None), Y14_0)    #14
    Initial(p16(b=None), Y18_0)    #18
    Initial(Rb(b=1,phos='u') % E2F(b=1), Y19_0)    #19
    Initial(Rb(b=1,phos='pp') % E2F(b=1), Y20_0)    #20
    Initial(E2F(b=None), Y21_0)    #21
    Initial(Rb(b=None,phos='pppp'), Y22_0)    #22
    Initial(Rb(b=None,phos='u'), Y23_0)    #23
    Initial(X(), Y26_0)    #26
    
    ##Shared
    
    Initial(Signal(), DDS_0)                                            #18
    Initial(SignalDamp(), DDS_0)            
    
#     Initial(Chk1(phos='p'), X1_0)                                       #0
#     Initial(Chk1(phos='u'), X1pre_0)                                    #1
    Initial(CycB(c=2) % CDK1_nuc(phos='p',b=1,c=2) % p21(b=1), X7_0)        #7    Inactive MPF (sequestered by p21)
    Initial(Cdc25(b=1, state= 'i', state1='C', phos= 'p') % x14_3_3(b=1), X10_0)    #10
    Initial(Wee1(phos= 'u'), X14_0)                                     #14
    Initial(Wee1(phos= 'p'), X15_0)                                     #15
    
    Initial(CycE(c=1) % CDK2(phos='u',b=None,c=1), Y6_0)    #6
    Initial(CycE(c=1) % CDK2(phos='p',b=None,c=1), Y7_0)    #7
    Initial(CycA(c=1) % CDK2(phos='u',b=None,c=1), Y8_0)    #8
    Initial(CycA(c=1) % CDK2(phos='p',b=None,c=1), Y9_0)    #9
    Initial(CycD(c=2) % CDK4_6(phos='p',b=1,c=2) % p27(b=1), Y11_0)    #11
    Initial(CycE(c=2) % CDK2(phos='p',b=1,c=2) % p27(b=1), Y12_0)    #12
    Initial(CycA(c=2) % CDK2(phos='p',b=1,c=2) % p27(b=1), Y13_0)    #13
    Initial(CycD(c=2) % CDK4_6(phos='p',b=1,c=2) % p21(b=1), Y15_0)    #15
    Initial(CycE(c=2) % CDK2(phos='p',b=1,c=2) % p21(b=1), Y16_0)    #16
    Initial(CycA(c=2) % CDK2(phos='p',b=1,c=2) % p21(b=1), Y17_0)    #17
    Initial(Mdm2(b=None), Y25_0)    #25
    Initial(I(), Y27_0)    #27

    Initial(CDK2(phos='u',b=None,c=None), shared_Y5_0)
    Initial(CycD(c=1) % CDK4_6(phos='p', b=None, c=1), shared_Y6_0)
    Initial(p27(b=None), shared_Y11_0)
    Initial(p53(b=None), shared_Y25_0)
    Initial(ATM_ATR(b=None), shared_Y27_0)
    Initial(Cdc25(b=None, state='i', state1='A', phos= 'u'), shared_Y28_0) # added state1 = 'A'
    Initial(Cdc25(b=None, state='a', state1='A', phos= 'u'), shared_Y29_0) # added state1 = "A"
    Initial(Chk1(phos='u'), shared_Y30_0) #'u'= inactive
    Initial(Chk1(phos='p'), shared_Y31_0) #'p' = active
    Initial(CycB(c=None), shared_Y33_0)
#     Initial(CDK1(phos='u', b=None, c=None), shared_Y34_0)
    Initial(CycB(c=1) % CDK1_cyto(phos='u', b=None, c=1), shared_Y35_0)
    Initial(CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1), shared_Y36_0)
    Initial(Cdc25(b=None, state='i', state1='C', phos= 'u'), shared_Y42_0) # added state1 = 'C'
    Initial(Cdc25(b=None, state='a', state1='C', phos= 'u'), shared_Y43_0) # added state1 = 'C'
    Initial(Cdc25(b=None, state='i', state1='C', phos= 'p'), shared_Y44_0) # added state1 = 'C'
    Initial(Cdc25(b=None, state='a', state1='C', phos= 'p'), shared_Y45_0) # added state1 = 'C'
    Initial(CycB(c=1) % CDK1_nuc(phos='u', b=None, c=1), shared_Y52_0)
    Initial(CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1), shared_Y53_0)
    
    ## New
    
    Initial(NF_Y(), shared_Y32_0) # Not sure
    Initial(B_Myb(state='i'), shared_Y40_0)
    Initial(B_Myb(state='a'), shared_Y41_0)
    Initial(APC(state='i', b=1) % Ccdc20(b=1), shared_Y48_0)
    Initial(APC(state='a', b=1) % Ccdc20(b=1), shared_Y49_0)
    Initial(APC(state='i', b=1) % Ccdh1(b=1), shared_Y50_0)
    Initial(APC(state='a', b=1) % Ccdh1(b=1), shared_Y51_0)
    

### ** Declare Observables **
def declare_observables():
    
    ## G2/M
    
    Observable("OBS_aCdc25C", Cdc25(b=None, state= 'a', state1='C', phos= 'u'))
    Observable("OBS_MPF", CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1))
    Observable("OBS_p53", p53(b=None))
    Observable("OBS_Wee1", Wee1(phos= 'u'))
    
    ## G1/S
    
    Observable("OBS_E2F", E2F(b=None))
    Observable("OBS_p21", p21(b=None))
    Observable("OBS_p16", p16(b=None))
    Observable("OBS_Rb", Rb(b=None, phos='u'))
    Observable("OBS_CycD_CDK4_6", CycD(c=1) % CDK4_6(phos='p',b=None,c=1))

    ## Shared
    
    Observable("signal", Signal())
    Observable("signal_damp", SignalDamp())
    Observable("OBS_Mdm2", Mdm2(b=None))
    Observable("OBS_p27", p27(b=None))
    Observable("OBS_CycE", CycE(c=None))
    Observable("OBS_aCycE_CDK2", CycE(c=1) % CDK2(phos='p', b=None, c=1))
    Observable("OBS_CycA", CycA(c=None))
    Observable("OBS_CycB", CycB(c=None))
    Observable("OBS_Int", I())
    
    ## New
    
    Observable("OBS_APC_Ccdc20", APC(state='a', b=1) % Ccdc20(b=1)) # state= 'a', 'i', both?

    
### ** Functions **
def declare_functions():
    
    ##  G2/M
    
    Expression("create_preMPF", sympify("G2_M_k9/(1 + G2_M_k31*OBS_p53)"))
    
    ## G1/S
    
    Expression("create_p16", sympify("G1_S_k41/((One + G1_S_k42*OBS_Rb) - (G1_S_k43*OBS_p16 + G1_S_k44*OBS_p16*OBS_CycD_CDK4_6))"))
    Expression("create_Rb", sympify("G1_S_k58/(One + G1_S_k59*OBS_p16)"))
    
    ## Shared
  
    Expression("create_Mdm2", sympify("G1_S_k66*OBS_Int**n/(G1_S_k65**9 + OBS_Int**n)"))
    Expression("create_Int", sympify("(G1_S_k70*OBS_p53*signal)/(One + G1_S_k71*OBS_p53*OBS_Mdm2)"))
    Expression("sig_deg", sympify("G1_S_k74 - G1_S_k73*(signal-signal_damp)"))
    Expression("kdamp_DDS0", sympify("G1_S_k75*DDS_0"))
    
   
### ** Rules **
def declare_rules():
    
## *** New Terms *** ##    
    Rule('E2F_Create_CycA', E2F(b=None) >> E2F(b=None) + CycA(c=None), shared_k130)
    Rule('NF_Y_Create_CycA', NF_Y() >> NF_Y() + CycA(c=None), shared_k75)
    Rule('aAPC_Ccdc20_Create_CycA', CycA(c=None) + APC(state='a', b=1) % Ccdc20(b=1) >> APC(state='a', b=1) % Ccdc20(b=1), shared_k126)
    Rule('aAPC_Ccdh1_Create_CycA', CycA(c=None) + APC(state='a', b=1) % Ccdh1(b=1) >> APC(state='a', b=1) % Ccdh1(b=1), shared_k127)
    Rule('iAPC_Ccdc20_Degrade_CycA', CycA(c=1) % CDK2(phos='u', b=None, c=1) + APC(state='i', b=1) % Ccdc20(b=1) >> CDK2(b=None, c=None) + APC(state='i', b=1) % Ccdc20(b=1), shared_k15)
    Rule('aAPC_Ccdc20_Degrade_CycA', CycA(c=1) % CDK2(phos='p', b=None, c=1) + APC(state='a', b=1) % Ccdc20(b=1) >> CDK2(b=None, c=None) + APC(state='a', b=1) % Ccdc20(b=1), shared_k14)
    Rule('iAPC_Ccdh1_Degrade_CycA', CycA(c=1) % CDK2(phos='u', b=None, c=1) + APC(state='i', b=1) % Ccdh1(b=1) >> CDK2(b=None, c=None) + APC(state='i', b=1) % Ccdh1(b=1), shared_k15)
    Rule('aAPC_Ccdh1_Degrade_CycA', CycA(c=1) % CDK2(phos='p', b=None, c=1) + APC(state='a', b=1) % Ccdh1(b=1) >> CDK2(b=None, c=None) + APC(state='a', b=1) % Ccdh1(b=1), shared_k14)
    Rule('aCdc25A_Activate_CycE_CDK2', CycE(c=1) % CDK2(phos='u', b=None, c=1) + Cdc25(b=None, state='a', state1='A', phos= 'u') >> CycE(c=1) % CDK2(phos='p', b=None, c=1) + Cdc25(b=None, state='a', state1='A', phos= 'u'), shared_k22)
    Rule('aCdc25A_Activate_CycA_CDK2', CycA(c=1) % CDK2(phos='u', b=None, c=1) + Cdc25(b=None, state='a', state1='A', phos= 'u') >> CycA(c=1) % CDK2(phos='p', b=None, c=1) + Cdc25(b=None, state='a', state1='A', phos= 'u'), shared_k28)
    Rule('p27_CycD_CDK4_6_Degrade_CycD_CDK4_6', CycD(c=2) % CDK4_6(phos='p',b=1,c=2) % p27(b=1) >> p27(b=None), shared_k21)
    Rule('Create_p21_CycB_CDK1', p21(b=None) + CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) >> CycB(c=2) % CDK1(phos='p',b=1,c=2) % p21(b=1), shared_k103) #Rate Constant
    Rule('ATM_ATR_Create_p53', ATM_ATR(b=None) >> p53(b=None), shared_k61)
    Rule('Create_Mdm2', None >> Mdm2(b=None), create_Mdm2)
    Rule('Degrade_ATM_ATR', ATM_ATR(b=None) >> None, shared_k79) #Rate Constant
    Rule('E2F_Create_iCdc25A', E2F(b=None) >> E2F(b=None) + Cdc25(b=None, state='i', state1='A', phos= 'u'), shared_k80)
    Rule('Deactivate_aCdc25A', Cdc25(b=None, state='a', state1='A', phos= 'u') >> Cdc25(b=None, state='i', state1='A', phos= 'u'), shared_k55) #Rate Constant
    Rule('aChk1_Degrade_iCdc25A', Chk1(phos='p') + Cdc25(b=None, state='i', state1='A', phos= 'u') >> Chk1(phos='p'), shared_k81) #Rate Constant
    Rule('aCycE_CDK2_Activate_Cdc25A', CycE(c=1) % CDK2(phos='p',b=None,c=1) + Cdc25(b=None, state='i', state1='A', phos= 'u') >> CycE(c=1) % CDK2(phos='p',b=None,c=1) + Cdc25(b=None, state='a', state1='A', phos= 'u'), shared_k82)
    Rule('aCycA_CDK2_Activate_Cdc25A', CycA(c=1) % CDK2(phos='p',b=None,c=1) + Cdc25(b=None, state='i', state1='A', phos= 'u') >> CycA(c=1) % CDK2(phos='p',b=None,c=1) + Cdc25(b=None, state='a', state1='A', phos= 'u'), shared_k82)
    Rule('aChk1_Degrade_aCdc25A', Chk1(phos='p') + Cdc25(b=None, state='a', state1='A', phos= 'u') >> Chk1(phos='p'), shared_k84) #Rate Constant
    Rule('Degrade_aCdc25', Cdc25(b=None, state='a', state1='A', phos= 'u') >> None, shared_k85) #Rate Constant
    Rule('Degrade1_aCdc25', Cdc25(b=None, state='a', state1='A', phos= 'u') >> None, shared_k86) #Rate Constant
    Rule('p21_CycB_CDK1_Degrade_p21', CycB(c=2) % CDK1(phos='p',b=1,c=2) % p21(b=1) >> CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1), shared_k104) #Rate Constant
    Rule('aCycB_CDK1_cyto_Activate_Cdc25C', CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) + Cdc25(b=None, state='i', state1='C', phos='u') >> CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='u'), shared_k110)
    Rule('aCycB_CDK1_nuc_Activate_Cdc25C', CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) + Cdc25(b=None, state='i', state1='C', phos='u') >> CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='u'), shared_k110)
    Rule('aChk1_Phos_iCdc25C', Chk1(phos='p') + Cdc25(b=None, state='i', state1='C', phos='u') >> Chk1(phos='p') + Cdc25(b=None, state='a', state1='C', phos='p'), shared_k111) #Rate Constant
    Rule('Degrade1_aCdc25C', Cdc25(b=None, state='a', state1='C', phos='u') + Cdc25(b=None, state='a', state1='C', phos='u') >> None, shared_k109)
    Rule('aChk1_Degrade_aCdc25C', Cdc25(b=None, state='a', state1='C', phos='u') + Chk1(phos='u') >> Chk1(phos='u'), shared_k113)
    Rule('Degrade2_aCdc25C', Cdc25(b=None, state='a', state1='C', phos='u') >> None, shared_k114)
    Rule('aCycB_CDK1_cyto_Activate_Cdc25CPs216', Cdc25(b=None, state='i', state1='C', phos='p') + CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) >> CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='p'), shared_k116)
    Rule('aCycB_CDK1_nuc_Activate_Cdc25CPs216', Cdc25(b=None, state='i', state1='C', phos='p') + CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) >> CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='p'), shared_k116)
    Rule('Wee1_Deactivate_CycB_CDK1_nuc', CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) + Wee1(phos='u') >> CycB(c=1) % CDK1_nuc(phos='u', b=None, c=1) + Wee1(phos='u'), shared_k133) #Rate Constant
    Rule('aCdc25C_Activate_CycB_CDK1_nuc', CycB(c=1) % CDK1_nuc(phos='u', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='u') >> CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='u'), shared_k134) #Rate Constant
    Rule('aCdc25CPs216_Activate_CycB_CDK1_nuc', CycB(c=1) % CDK1_nuc(phos='u', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='p') >> CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='p'), shared_k134) #Rate Constant
    Rule('aAPC_Ccdc20_Degrade_iCycB_CDK1_nuc', APC(state='a', b=1) % Ccdc20(b=1) + CycB(c=1) % CDK1_nuc(phos='u', b=None, c=1) >> APC(state='a', b=1) % Ccdc20(b=1), shared_k136)
    Rule('aAPC_Ccdc20_Degrade_aCycB_CDK1_nuc', APC(state='a', b=1) % Ccdh1(b=1) + CycB(c=1) % CDK1_nuc(phos='u', b=None, c=1) >> APC(state='a', b=1) % Ccdh1(b=1), shared_k136)
    Rule('aAPC_Ccdh1_Degrade_iCycB_CDK1_nuc', CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) + CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) >> CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1), shared_k131)
    Rule('aAPC_Ccdh1_Degrade_aCycB_CDK1_nuc', APC(state='a', b=1) % Ccdc20(b=1) + CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) >> APC(state='a', b=1) % Ccdc20(b=1), shared_k135)
    Rule('aCycB_CDK1_cyto_Create_aCycB_CDK1_nuc', APC(state='a', b=1) % Ccdh1(b=1) + CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) >> APC(state='a', b=1) % Ccdh1(b=1), shared_k137)
    
## *** New Species *** ##
    Rule('aCycA_CDK2_Create_NF_Y', CycA(c=1) % CDK2(phos='p', b=None, c=1) >> CycA(c=1) % CDK2(phos='p', b=None, c=1) + NF_Y(), shared_k89)
    Rule('Degrade_NF_Y', NF_Y() >> None, shared_k90)
    Rule('NF_Y_Create_CycB', NF_Y() >> NF_Y() + CycB(c=None), shared_k91)
    Rule('iCycB_CDK1_cyto_Degrade_CDK1', CycB(c=1) % CDK1_cyto(phos='u', b=None, c=1) >> CycB(c=None), shared_k95)
    Rule('Degrade_CycB', CycB(c=None) >> None, shared_k92)
    Rule('CycB_CDK1_Create CycB_CDK1_cyto', CycB(c=None) + CDK1(phos='u', b=None, c=None) >> CycB(c=1) % CDK1_cyto(phos='u', b=None, c=1), shared_k93)
    Rule('aAPC_Ccdc20_Degrade_CycB', CycB(c=None) + APC(state='a', b=1) % Ccdc20(b=1) >> APC(state='a', b=1) % Ccdc20(b=1), shared_k128)
    Rule('aAPC_Ccdh1_Degrade_CycB', CycB(c=None) + APC(state='a', b=1) % Ccdh1(b=1) >> APC(state='a', b=1) % Ccdh1(b=1), shared_k129)
    Rule('iCycB_CDK1_cyto_Degrade_CycB', CycB(c=1) % CDK1_cyto(phos='u', b=None, c=1) >> CDK1(phos='u', b=None, c=None), shared_k94)
    Rule('aAPC_Ccdc20_Degrade_iCycB_CDK1_cyto', CycB(c=1) % CDK1_cyto(phos='u', b=None, c=1) + CycB(c=None) + APC(state='a', b=1) % Ccdc20(b=1) >> CDK1(phos='u', b=None, c=None) + CycB(c=None) + APC(state='a', b=1) % Ccdc20(b=1), shared_k97)
    Rule('aAPC_Ccdh1_Degrade_iCycB_CDK1_cyto', CycB(c=1) % CDK1_cyto(phos='u', b=None, c=1) + CycB(c=None) + APC(state='a', b=1) % Ccdh1(b=1) >> CDK1(phos='u', b=None, c=None) + CycB(c=None) + APC(state='a', b=1) % Ccdh1(b=1), shared_k97)
    Rule('aAPC_Ccdc20_Degrade_aCycB_CDK1_cyto', CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) + CycB(c=None) + APC(state='a', b=1) % Ccdc20(b=1) >> CDK1(phos='u', b=None, c=None) + CycB(c=None) + APC(state='a', b=1) % Ccdc20(b=1), shared_k98)
    Rule('aAPC_Ccdh1_Degrade_aCycB_CDK1_cyto', CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) + CycB(c=None) + APC(state='a', b=1) % Ccdh1(b=1) >> CDK1(phos='u', b=None, c=None) + CycB(c=None) + APC(state='a', b=1) % Ccdh1(b=1), shared_k98)
    Rule('Deactivate_CycB_CDK1_cyto', CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) >> CycB(c=1) % CDK1_cyto(phos='u', b=None, c=1), shared_k96)
    Rule('aCdc25C_Activate_CycB_CDK1_cyto', CycB(c=1) % CDK1_cyto(phos='u', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='u') >> CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='u'), shared_k95)
    Rule('aCdc25CPs216_Activate_CycB_CDK1_cyto', CycB(c=1) % CDK1_cyto(phos='u', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='p') >> CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) + Cdc25(b=None, state='a', state1='C', phos='p'), shared_k95)
    Rule('aCycB_CDK1_nuc_to_cyto', CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) >> CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1), shared_k132)
    Rule('Degrade_aCycB_CDK1_cyto', CycB(c=1) % CDK1_cyto(phos='p', b=None, c=1) >> None, shared_k96)
    Rule('E2F_Create_iB_Myb', E2F(b=None) >> B_Myb(state='i'), shared_k105)
    Rule('aCycA_CDK2_Activate_B_Myb', CycA(c=1) % CDK2(phos='p', b=None, c=1) + B_Myb(state='i') >> CycA(c=1) % CDK2(phos='p', b=None, c=1) + B_Myb(state='a'), shared_k106)
    Rule('Degrade_aB_Myb', B_Myb(state='a') >> None, shared_k107)
    Rule('aAPC_Ccdh1_Deactivate_APC_Ccdc20', APC(state='a', b=1) % Ccdc20(b=1) + APC(state='a', b=1) % Ccdh1(b=1) >> APC(state='i', b=1) % Ccdc20(b=1) + APC(state='a', b=1) % Ccdh1(b=1), shared_k122)
    Rule('aCycB_CDK1_nuc_Activate_APC_Ccdc20', APC(state='i', b=1) % Ccdc20(b=1) + CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) >> APC(state='a', b=1) % Ccdc20(b=1) + CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1), shared_k123)
    Rule('aCycB_CDK1_nuc_Deactivate_APC_Ccdh1', APC(state='a', b=1) % Ccdh1(b=1) + CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1) >> APC(state='i', b=1) % Ccdh1(b=1) + CycB(c=1) % CDK1_nuc(phos='p', b=None, c=1), shared_k124)
    Rule('aCycA_CDK2_Deactivate_APC_Ccdh1', APC(state='a', b=1) % Ccdh1(b=1) + CycA(c=1) % CDK2(phos='p', b=None, c=1) >> APC(state='i', b=1) % Ccdh1(b=1) + CycA(c=1) % CDK2(phos='p', b=None, c=1), shared_k124)
    Rule('Activate_APC_Ccdh1', APC(state='i', b=1) % Ccdh1(b=1) >> APC(state='a', b=1) % Ccdh1(b=1), shared_k125)