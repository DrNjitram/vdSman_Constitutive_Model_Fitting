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

path_read="./experimental_data/2024-08-06/20240905_P60.xlsx"
freqvals_frfit= [ 0.1, 0.1585, 0.2512, 0.3981, 0.631, 1, 1.585, 2.512, 3.981, 6.31, 10]
color_palette = plt.cm.hsv(np.linspace(0, 1, 16))
str_grw_shear_rates = [0.063, 0.1, 0.17, 0.28, 0.45, 0.73, 1.2,\
                1.96, 3.21, 5.25, 8.58, 14,23,37,61]
saos_strain_perc = 0.1
ndata_cycle =257 
font = {'family' : 'arial',
        'weight' : 'bold',
        'size'   : 24}

plt.rcParams.update({'font.size': 24})

cyclesdata=5
ndata_percycle =257