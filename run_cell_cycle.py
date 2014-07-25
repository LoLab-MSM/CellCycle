from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify
from scipy import constants

from G1_S_v2 import *
from G2_M_v2 import *
from Cell_Cycle_Shared import *


declare_monomers()
declare_parameters()
declare_initial_conditions()
declare_observables()
declare_functions()
declare_rules()
     
generate_equations(model, verbose=True)

### *** Checking and Printing Everything to Screen ***

# print len(model.rules)
# print len(model.initial_conditions)
# print len(model.reactions)
# print len(model.species)
# quit()
  
# for monomers in model.monomers:
#     print monomers
# print
#   
# for parameters in model.parameters:
#     print parameters
# print
#  
# for initial_conditions in model.initial_conditions:
#     print initial_conditions
# print
#  
# for obs in model.observables:
#     print obs, ":", obs.species, ",", obs.coefficients
#     obs_string = ''
#     for i in range(len(obs.coefficients)):
#         if i > 0: obs_string += " + "
#         obs_string += "__s"+str(obs.species[i])
#         if obs.coefficients[i] > 1:
#              obs_string += "*"+str(obs.coefficients[i])
#     print obs_string
#   
# for rules in model.rules:
#     print rules
# print
#  
# for i in range(len(model.species)):
#     print str(i)+":", model.species[i]
# print
#  
# for i in range(len(model.odes)):
#     print str(i)+":", model.odes[i]
# print
# 
# for x in model.parameters_initial_conditions():
#     print x, ":", x.value
# print 
# 
# for x in model.parameters_unused():
#     print x, ":", x.value
# print 
# 
# for x in model.parameters_rules():
#     print x
#    
# quit()
#  
# from pysb.generator.bng import BngGenerator
# print BngGenerator(model).get_content()

t = linspace(0,6000,600)

## ** Set DNA Damage (None) **
set_dna_damage(0.0)
y = odesolve(model,t,verbose=True)
 
pl.figure()
for obs in ["OBS_p27", "OBS_CycE", "OBS_CycA", "OBS_CycB", "OBS_CycE_CDK2"]:
    pl.plot(t, y[obs], label=re.match(r"OBS_(\w+)", obs).group(1), linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (DNA Damage = 0)", fontsize=22)
pl.savefig("Whole Cell Cycle No DNA Damage1.png", format= "png")

pl.figure()
pl.plot(t, y["OBS_APC_Ccdc20"], label="APC_Cdc20", linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level (APC_Cdc20)", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (DNA Damage = 0)", fontsize=22)
pl.savefig("Whole Cell Cycle No DNA Damage2.png", format= "png")

## ** Set DNA Damage (Low) **
set_dna_damage(0.002)
y = odesolve(model,t,verbose=True)
 
pl.figure()
for obs in ["OBS_p27", "OBS_CycE", "OBS_CycA", "OBS_CycB", "OBS_CycE_CDK2"]:
    pl.plot(t, y[obs], label=re.match(r"OBS_(\w+)", obs).group(1), linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (DNA Damage = 0.002)", fontsize=22)
pl.savefig("Whole Cell Cycle Low DNA Damage1.png", format= "png")

pl.figure()
pl.plot(t, y["OBS_APC_Ccdc20"], label="APC_Cdc20", linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level (APC_Cdc20)", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (DNA Damage = 0.002)", fontsize=22)
pl.savefig("Whole Cell Cycle Low DNA Damage2.png", format= "png")

## ** Set DNA Damage (Medium) **
set_dna_damage(0.004)
y = odesolve(model,t,verbose=True)
 
pl.figure()
for obs in ["OBS_p27", "OBS_CycE", "OBS_CycA", "OBS_CycB", "OBS_CycE_CDK2"]:
    pl.plot(t, y[obs], label=re.match(r"OBS_(\w+)", obs).group(1), linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (DNA Damage = 0.004)", fontsize=22)
pl.savefig("Whole Cell Cycle Medium DNA Damage1.png", format= "png")

pl.figure()
pl.plot(t, y["OBS_APC_Ccdc20"], label="APC_Cdc20", linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level (APC_Cdc20)", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (DNA Damage = 0.004)", fontsize=22)
pl.savefig("Whole Cell Cycle Medium DNA Damage2.png", format= "png")

## ** Set DNA Damage (High) **
set_dna_damage(0.008)
y = odesolve(model,t,verbose=True)
 
pl.figure()
for obs in ["OBS_p27", "OBS_CycE", "OBS_CycA", "OBS_CycB", "OBS_CycE_CDK2"]:
    pl.plot(t, y[obs], label=re.match(r"OBS_(\w+)", obs).group(1), linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (DNA Damage = 0.008)", fontsize=22)
pl.savefig("Whole Cell Cycle High DNA Damage1.png", format= "png")

pl.figure()
pl.plot(t, y["OBS_APC_Ccdc20"], label="APC_Cdc20", linewidth=3)
pl.legend(loc='upper left', prop={'size': 16})
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level (APC_Cdc20)", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (DNA Damage = 0.008)", fontsize=22)
pl.savefig("Whole Cell Cycle High DNA Damage2.png", format= "png")

## ** Set DNA Damage (Extreme) **
set_dna_damage(0.016)
y = odesolve(model,t,verbose=True)
 
pl.figure()
for obs in ["OBS_p27", "OBS_CycE", "OBS_CycA", "OBS_CycB", "OBS_CycE_CDK2"]:
    pl.plot(t, y[obs], label=re.match(r"OBS_(\w+)", obs).group(1), linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (DNA Damage = 0.016)", fontsize=22)
pl.savefig("Whole Cell Cycle Extreme DNA Damage1.png", format= "png")

pl.figure()
pl.plot(t, y["OBS_APC_Ccdc20"], label="APC_Cdc20", linewidth=3)
pl.legend(loc='upper left', prop={'size': 16})
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level (APC_Cdc20)", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (DNA Damage = 0.016)", fontsize=22)
pl.savefig("Whole Cell Cycle Extreme DNA Damage2.png", format= "png")

pl.show()