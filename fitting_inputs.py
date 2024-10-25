import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit, minimize
import math as math
import os 
import warnings
from sklearn.metrics import r2_score 
import copy
import subprocess
from PIL import Image
import scipy.stats as stats
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from SALib.sample import saltelli
from SALib.analyze import sobol
import json
import datetime
from scipy.stats import t
from pyDOE import lhs
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit, minimize
import math as math
import os 
import warnings
import copy
import subprocess
from PIL import Image
from scipy.optimize import approx_fprime
import importlib
import textwrap
import sys

# /\_/\  
#( o.o ) Hello, This is the input file for the data fitting! 
# > ^ <


# ╔═════════════════════════════════════════════════════════════════╗
# ║         Information on Data Structures & Some Parameters        ║
# ╚═════════════════════════════════════════════════════════════════╝

#In the demo the LAOS data is stored in the following format
#  data_exp_avg =[ freq1-strain1, freq1-strain2, ...,
#                  freq2-strain1, freq2-strain2, ..., 
#                  freq3-strain1, freq2-strain2, ..., 
#                   ...
#                  freqN-strainN, freqN-strainN+1, ... ]

# where each freq*-strain* [ [time][stress][strain][shear rate] ]
# the following parameters describe the index of the time, strain, stress and shear rate
col_index_time = 0
col_index_stress = 1
col_index_strain = 2
col_index_shearrate = 3
#The sequential order of the data_exp_avg must follow the order of the strain and frequency values provided
#in the following arrays
strain_values_perc_nonlin = [120,100,50,30,15,10,5]
strain_values_nonlin = np.divide(strain_values_perc_nonlin,100) #strain should be actual not percentage
freqvals = [2,4,8]
#The demo set is [ freq2Hz+Strain120% , freq2Hz+Strain100%, freq2Hz+Strain50%, ... , 
#                  freq4Hz+Strain120% , freq4Hz+Strain100%, freq4Hz+Strain50%, ... , 
#                  freq8Hz+Strain120% , freq8Hz+Strain100%, freq8Hz+Strain50%, ... ]



#The stress growth data used in the demo
# data_exp_stgrw = [ shear_rate1, shear_rate2, ..., shear_rateN]
# shear_rateN = [  [time][strain{%}][shear stress]   ]
i_strain = 1
i_stress = 2
#The shear rates applied in the stress growth testing. It must follow the order in the data sheet. 
str_grw_shear_rates = [0.063, 0.1, 0.17, 0.28, 0.45, 0.73, 1.2,\
                1.96, 3.21, 5.25, 8.58, 14,23,37,61]
#While building a flow curve, which strain (in percentage) you want to collect the data?
target_strain =120 # to extract the stress at


#The SAOS data udes in the demo
# it is a frequency sweep data at an amplitude an various frequencies. 
#the dataset array is structureed as follows
#The demo set is 
# data_exp_frswp_waveform_avg= [ [ time , stress, strain , shear_rate ] ]
# time  = [ freq0_StrainLin , freq1_StrainLin, freq2_StrainLin, ... , freqN_StrainLin, ]
#for instance: data_exp_frswp_waveform_avg[0][0][1] has the stress data from all the frequencies sequentially. 
strain_lin = 0.001
freqvals_frfit= [ 0.1, 0.1585, 0.2512, 0.3981, 0.631, 1, 1.585, 2.512, 3.981, 6.31, 10]
         

#There are also datasets from strain sweeps 

# ╔═════════════════════════════════════════════════════════════════╗
# ║                         General Inputs                          ║
# ╚═════════════════════════════════════════════════════════════════╝

#The general color palette and plot details used for all the fitting in the notebook
color_palette = plt.cm.tab20(np.linspace(0, 1, 20))
markersize=10

#Naming for the input files
current_time = datetime.datetime.now()
formatted_time = current_time.strftime('%Y_%m_%d_%H_%M')

#How many cycles will the model be running for? the fitting is always by using the last cycle. 
cycles = 5

#The highest shear rate viscosy to prevent errors in calcualtion. 
eta_inf = 1e-9

#Number of elements and time
nelements_cycle = 257 #from 0 to final data
ntimesteps_cycle = nelements_cycle -1


# ╔═════════════════════════════════════════════════════════════════╗
# ║         Running the HB fitting for HB informed fitting          ║
# ╚═════════════════════════════════════════════════════════════════╝

#Option 1: Read a previous run case
run_HB = 'no'
readHB = 'HBfit_2024_10_24_18_12' #read this file if you do not run a new case

#Option2: Run the HB fitting from scratch
#run_HB = 'yes'
extract_HB_param_from_laos = 'yes' #use LAOS data?
extract_HB_param_from_external = 'yes' #provide an external result?
extract_HB_param_from_stress_growth = 'yes' #use stress growth data?
#There is one general output file for all cases in this section

#Which data to take from LAOS set while building the flow curves from LAOS data?
#Provide the index for the data in the data_exp_all. 
#These indicate the highest strains from the each 3 amplitudes in the demo set. 
data_index_for_laos_built_flow_curves = [0,1,2,7,8,9,14,15,16]
skip_first = 3 #skip the first ** of laos extracted flow data while evaluating

#External parameters
HB_external_params = [1000, 10, 0.2] #sigma_y, dot_gamma_cr, n
error_HBext = [0,0,0] #sigma_y, dot_gamma_cr, n


# ╔═════════════════════════════════════════════════════════════════╗
# ║                       Running the SAOS Fitting                  ║
# ╚═════════════════════════════════════════════════════════════════╝

#Option 1:
#runSAOS = 'yes'

#Option 2:
runSAOS = 'no'
readSAOS = 'SAOSfit_2024_10_24_18_12'

#These inputs are used while performing the fitting
#The initial guess for the SAOS fitting
modulus = [0.5] 
t_rel=[0.5] 

#Bounds
bounds_mod0 = [2,6] #10^2 - 10^6
bounds_trel0 = [-2,2] #10^-2 - 10^2
bounds_mod0_n = (0,1)  #the values change between 0 and 1 after normalization. 
bounds_trel0_n = (0,1) #the values change between 0 and 1 after normalization. 

#Form the input arrays using the inputs provided above
bounds_lin_norm = [bounds_mod0, bounds_trel0]
bounds_lin = [bounds_mod0_n, bounds_trel0_n] 
init_guess_lin =  [modulus[0], t_rel[0]]



# ╔═════════════════════════════════════════════════════════════════╗
# ║                       Running the LAOS Fitting                  ║
# ╚═════════════════════════════════════════════════════════════════╝

#Option 1:
#runLAOS = 'yes'

#Option 2:
runLAOS = 'no'
readLAOS = 'LAOSfit_all_2024_10_24_18_12'

#The initial guess for the LAOS fitting
tau_y_log=[0.5] 
nexp=[0.5]   
alpha=[0.5] 
gammadot_cr_log=[0.5] 
I1c_log=[0.5] 

#Bounds
bounds_tauy = [2,4] #log
bounds_nexp = [0.1, 0.8] #linear
bounds_alpha =  [0.1, 0.8] #linear
bounds_gammadotcr = [-2,3] #log
bounds_I1c = [-4,3] #log
bounds_tauy_n = (0,1)
bounds_nexp_n = (0,1)
bounds_alpha_n =  (0,1)
bounds_gammadotcr_n = (0,1)
bounds_I1c_n = (0,1)

#Form the input arrays using the inputs provided above
bounds_nonlin_norm = [ bounds_tauy, bounds_nexp, bounds_alpha, bounds_gammadotcr, bounds_I1c ] 
bounds_nonlin = [ bounds_tauy_n, bounds_nexp_n, bounds_alpha_n, bounds_gammadotcr_n, bounds_I1c_n ] 
init_guess_nonlin = [tau_y_log[0], nexp[0], alpha[0], gammadot_cr_log[0],I1c_log[0] ]
 
## Herschel-Bulkley Parameters infromed fitting
init_guess_hb = [alpha[0], I1c_log[0] ]
bounds_nonlin_norm_hb = [  bounds_alpha, bounds_I1c ] 
bounds_nonlin_hb = [ bounds_alpha_n, bounds_I1c_n ] 


# ╔═════════════════════════════════════════════════════════════════╗
# ║         Prediction by Using the Fitted Parameters               ║
# ╚═════════════════════════════════════════════════════════════════╝

shear_rate_vals = np.logspace(-2,3,30)
ntimesteps_steady = 1e5










