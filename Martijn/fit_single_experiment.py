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

import functions_fitting
from functions_fitting import *
importlib.reload(functions_fitting)

import fitting_inputs
from fitting_inputs import *
importlib.reload(fitting_inputs)

from sklearn.neighbors import LocalOutlierFactor
ndata_cycle = 90

def show_data(data_exp):
    fig, axs = plt.subplots()

    axs.set_xlabel('Index/ Time [s]')

    axs.set_ylabel('Shear Rate [1/s]')

    axs.set_title('1 Hz')

    for i in range(len(data_exp)):
        x1 = np.arange(ndata_cycle)  # data_exp[i][:, 0]
        y1 = data_exp[i][:, 1][-ndata_cycle:]

        # Plot data in the first subplot
        axs.plot(x1, y1, marker='o', color=color_palette[i], label='Data 1', markersize=2, linewidth=0)

    # Adjust layout for better spacing
    plt.tight_layout()

    freq = 1
    ampl = max([np.max(data[:, 1]) for data in data_exp]) * 15
    t = np.linspace(0, ndata_cycle / freq, ndata_cycle)
    strains = []
    shrates = []
    for elem in t:
        shrate = ampl / 100 * freq * 2 * math.pi * np.cos(freq * elem * 2 * math.pi)

        strain = ampl / 100 * np.sin(freq * elem * 2 * math.pi)
        strains.append(strain)
        shrates.append(shrate)
    axs.plot(t, shrates, marker='o', color='blue', markersize=0, linewidth=1)

    # Display the plots
    plt.show()

def reject_outliers(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else np.zeros(len(d))
    return data[s<m]

def LAOS(data_exp_avg, exp):
    col_index_stress = 1
    col_index_shearrate = 0

    # Number of elements and time
    nelements_cycle = 90  # from 0 to final data


    # Figure settings
    fig, axs = plt.subplots(figsize=(10, 8))
    axs.set_xscale('log')
    axs.set_yscale('log')
    #axs.set_xlim(1E-4, 1E-1)
    #axs.set_ylim(1E1, 1E5)
    axs.set_xlabel("$\\dot{{\\gamma}}\\;\\;[s^{-1}]$")
    axs.set_ylabel("$\\sigma_{xy}\\;\\;[Pa]$")
    #Extract the portion of the data for the steady shear curve build and plot it
    x_shr,y_str = extract_steady_from_laos(data_exp_avg, nelements_cycle,
                                           col_index_shearrate=col_index_shearrate, col_index_stress=col_index_stress)
    lof = LocalOutlierFactor(n_neighbors=20, contamination=0.01)
    hbfit_x = []
    hbfit_y = []
    skip_first = 1
    for i in range(len(x_shr)):
        if x_shr[i][skip_first:]:
            x = x_shr[i][skip_first:] ; y=y_str[i][skip_first:]
            xy = np.array([x, y]).T
            if len(x)>5:
                good = lof.fit_predict(xy) == 1
                xy = xy[good, :]

            axs.errorbar(xy[:, 0], xy[:, 1], marker='o', linewidth=0, markersize = markersize, color=color_palette[i],
                        markeredgecolor='black' )

            for elem in xy[:, 0]: hbfit_x.append(elem)
            for elem in xy[:, 1]: hbfit_y.append(elem)


    try:
        #fit the data to the HB model
        popt_laos, pconv_laos= curve_fit(HB, hbfit_x, hbfit_y, p0=[1000, 5, 0.2], bounds=[[1E-10, 1E-2, 1E-4],[1E6, 1E6, 1.5]], maxfev=2000)
        dot_gammacr_laos = (popt_laos[0] / (popt_laos[1]))**(1/popt_laos[2])
        model_data_HB_laos = HB(np.logspace(-3,3,50), popt_laos[0],popt_laos[1],popt_laos[2])
        plt.plot(np.logspace(-3,3,50), model_data_HB_laos,  linewidth=2, color='black', linestyle='--')
        HBparam_laos = [popt_laos[0], popt_laos[2],  dot_gammacr_laos]
        HBparam_laos = list(map(float, HBparam_laos))
        # calcualte errors
        perr_laos = np.sqrt(np.diag(pconv_laos))
        HBlaos_dotgam_pluserr = (popt_laos[0] / (popt_laos[1] + perr_laos[1])) ** (1 / popt_laos[2])
        HBlaos_dotgam_minuserr = (popt_laos[0] / (popt_laos[1] - perr_laos[1])) ** (1 / popt_laos[2])
        error_HBlaos = [perr_laos[0], perr_laos[2], HBlaos_dotgam_pluserr, HBlaos_dotgam_minuserr]
        error_HBlaos = list(map(float, error_HBlaos))

        print(f'LAOS fitted HB: TAU - n - dotgammacr : {HBparam_laos} error {error_HBlaos}')

        plt.title(f"{exp} {error_HBlaos}")

        plt.show()
        return HBparam_laos+error_HBlaos
    except Exception as e:
        print(e)

        plt.title(f"{exp}")


        plt.show()
        return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]



def read_data(path, block, i):
    path_read = rf"{path}{block+1}"
    data_exp, amplitudes = read_data_CCR_cycle(path_read, i)
    return data_exp

with open("results.csv", "w") as f:
    # measurements = [10, 10, 10, 13, 8]
    # path = r"D:\Data\Rheology\CCR\202411 DoE pH T M Alph 8\Block "
    # for block, length in enumerate(measurements):
    #     for exp in range(length):
    #         data_exp_avg = read_data(path, block, exp)
    #         #show_data(data_exp_avg)
    #         result = LAOS(data_exp_avg, exp)
    #         f.write(f"{block},{exp},Alpha8,{",".join(map(str, result))}\n")
    #         print(f"{block},{exp},Alpha8,{",".join(map(str, result))}")

    measurements = [0, 11, 0, 2, 2]
    path = r"D:\Data\Rheology\CCR\202409 DOE pH T M Arcon F\Block"
    for block, length in enumerate(measurements):
        if length:
            for exp in range(length):
                data_exp_avg = read_data(path, block, exp)
                show_data(data_exp_avg)
                result = LAOS(data_exp_avg, exp)
                f.write(f"{block},{exp},ArconF,{",".join(map(str, result))}\n")
                print(f"{block},{exp},ArconF,{",".join(map(str, result))}")




