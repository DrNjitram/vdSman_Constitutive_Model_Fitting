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

color_palette = plt.cm.tab20(np.linspace(0, 1, 20))
markersize=10

#Stress Growth Data if exsits
target_strain =120 # to extract the stress at
skip_first = 3 #skip the first * of laos extracted flow data

## Running the HB fitting for HB informed fitting ##

#Option 1: Read a previous run case
run_HB = 'no'
readHB = 'HBfit_2024_10_22_14_22' #read this file

# #Option2: Run the HB fitting from scratch
# run_HB = 'yes'
current_time = datetime.datetime.now()
formatted_time = current_time.strftime('%Y_%m_%d_%H_%M')

extract_HB_param_from_laos = 'yes' #use LAOS data?
data_index_for_laos_built_flow_curves = [0,1,2,7,8,9,14,15,16]

extract_HB_param_from_stress_growth = 'no' #use stress growth data?
str_grw_shear_rates = [0.063, 0.1, 0.17, 0.28, 0.45, 0.73, 1.2,\
                1.96, 3.21, 5.25, 8.58, 14,23,37,61]

extract_HB_param_from_external = 'no' #provide an external result?
HB_external_params = [1000, 10, 0.2] #sigma_y, dot_gamma_cr, n
error_HBext = [0,0,0] #sigma_y, dot_gamma_cr, n



## Running the SAOS Fitting ##
strain_lin = 0.001
freqvals_frfit= [ 0.1, 0.1585, 0.2512, 0.3981, 0.631, 1, 1.585, 2.512, 3.981, 6.31, 10]
cycles = 5
eta_inf = 1e-9
ntimesteps_cycle = 256
#Inputs
modulus = [0.5] 
t_rel=[0.5] 
#Bounds
bounds_mod0 = [2,6]
bounds_trel0 = [-2,2]
bounds_mod0_n = (0,1)
bounds_trel0_n = (0,1)
bounds_lin_norm = [bounds_mod0, bounds_trel0]
bounds_lin = [bounds_mod0_n, bounds_trel0_n] 

#Option 1:
#runSAOS = 'yes'

#Option 2:
runSAOS = 'no'
readSAOS = 'SAOSfit_2024_10_23_09_26'



















