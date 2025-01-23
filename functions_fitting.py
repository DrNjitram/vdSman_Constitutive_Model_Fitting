from collections import namedtuple

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


# ───────────────────────────────────────────────────────────────────
#   Data reading scripts
# ───────────────────────────────────────────────────────────────────
def read_data(n_header, n_rows, header_list, what_to_collect_time,\
               what_to_collect_single,path_read, sheet_name,\
               sample_list, rows_btw_data, cyclesdata, ndata_percycle):
    
    data = []
    data_single = []
    i=0
    for sample in sample_list:
        sample_strain_single = []
        sample_strain = []

        var = pd.read_excel(io=path_read, sheet_name=sheet_name, header=n_header, names=header_list, nrows=n_rows)
        

        for elem1 in what_to_collect_time:
            
            elem = var[elem1].tolist()
            sample_strain.append(elem)

        for elem2 in what_to_collect_single:
            a= var[elem2].tolist()
            
            for k in range(0,int(cyclesdata)):
                sample_strain_single.append(a[k*int(ndata_percycle)])
                
        n_header = n_header + n_rows + rows_btw_data

        data.append(sample_strain)
        data_single.append(sample_strain_single)
        i+=1

    return data, data_single

def read_data1(n_header, n_rows, header_list, what_to_collect_time, path_read, sheet_name, sample_list, rows_btw_data):
    
    data = []
    i=0
    for sample in sample_list:
        sample_strain = []
        var = pd.read_excel(io=path_read, sheet_name=sheet_name, header=n_header, names=header_list, nrows=n_rows)
    
        for elem in what_to_collect_time:
            sample_strain.append(var[elem].tolist())
                
        n_header = n_header + n_rows + rows_btw_data

        data.append(sample_strain)
        i+=1

    return data



def read_data_CCR(n_header, header_list, what_to_collect, path_read, cols_btw_data):
    index = 0
    measurements = []
    Info = namedtuple("Measurement",
                             "Mixture, Customer, Instrument, Material, Comments, Info, Batch_NB, Order_NB, N3, N4, N5, Frequency, Temperature, InterID")
    infos = []
    while True:
        try:
            measurement = []
            info = pd.read_excel(io=path_read, header=None, skiprows=1, nrows=14,
                                 usecols=[(1 if not index else 0) + index * cols_btw_data], names=["Value"])
            data = pd.read_excel(io=path_read, header=None, names=header_list, skiprows=n_header,
                                 usecols=list(range(index * cols_btw_data, (index + 1) * cols_btw_data)))

            for header in what_to_collect:
                if header == "index":
                    measurement.append(np.arange(len(data["Strain %"])))
                else:
                    measurement.append(np.array(data[header]))

            measurements.append(measurement)
            infos.append(Info(*info["Value"]))

            index += 1
        except Exception as e:
            if not measurements:
                print(e)
                raise e
            break


    return infos, measurements

# ───────────────────────────────────────────────────────────────────
#   HB fitting for HB informed fitting
# ───────────────────────────────────────────────────────────────────    
def shear_rate_sweep(target_strain, data_exp_stgrw, i_strain, i_stress ):
    
    stress = []
    indexes = []
    for rate in data_exp_stgrw:
        diff_list = np.abs(np.subtract(rate[i_strain],[target_strain]*len(rate[i_strain])))
        min_diff = min(diff_list)

        for i in range(0,len(diff_list)):
            if diff_list[i] == min_diff:
                index = i

        stress.append(rate[i_stress][index])
        indexes.append(index)

    return stress, indexes

def extract_steady_from_laos(data_exp_avg,nelements_cycle, col_index_shearrate, col_index_stress):

    ultimate_shrate = []
    ultimate_stress = []

    for k in range(0,len(data_exp_avg)):

        x = data_exp_avg[k][col_index_shearrate][-nelements_cycle:] #shear rate [1/s] -Martijn: Strain rate gamma dot [1/s]
        y = data_exp_avg[k][col_index_stress][-nelements_cycle:] #shear stress [Pa] -Martijn: Stress Sigma [Pa]

        #extract the positive values of both sher stress and shearrate
        xpos = []
        ypos = []
        for i in range(0,len(x)):
            if y[i]>0 and x[i]>0: xpos.append(x[i]) ; ypos.append(y[i])

        #find the min of the shear rate
        minshr = min(xpos)
        for i in range(0,len(ypos)):
            if xpos[i] == minshr: minindex_x=i ; minindex_y=i

        #Stress corresponding to min shear rate --> increase shear rate to spot where a permanent decay in stress starts
        ydat = []
        xdat = []
        for i in np.linspace(minindex_x,0,minindex_x+1):
            i=int(i)
            if ypos[i-1]< ypos[i] and ypos[i-2] < ypos[i]: break
            else: ydat.append(ypos[i]) ; xdat.append(xpos[i])

        ultimate_shrate.append(xdat)  
        ultimate_stress.append(ydat)  

    return ultimate_shrate, ultimate_stress

# ───────────────────────────────────────────────────────────────────
#   Some general functions
# ───────────────────────────────────────────────────────────────────
def HB(x, tau0, K, n):
    return tau0 + K*np.power(x,n)

def power_law(x,k,n): 
    return np.multiply(k,np.power(x,n))

def linear(x,a,b):
    return np.multiply(x,a)+b

# ───────────────────────────────────────────────────────────────────
#   Constitutive model
# ───────────────────────────────────────────────────────────────────
def Ruud_Single_Cycle_Oldroydb(strain_lin,freq,params, ncycles, ntimesteps_percycle, gammadot_data):

    A_11 = 1 
    A_12 = 0
    A_22 = 1

    G,tau_M = params
    
    dt = 1/(freq*ntimesteps_percycle)
    G1 = []
    G2 = []

    sigma = []
    time = []
    for cycle in range(0,ncycles):

        sum_G1 = 0
        sum_G2 = 0 
        for t in range (0,ntimesteps_percycle):
     
            gammadot = gammadot_data[t]

            #Solve for the structure tensor
            A_11 += dt*(2*gammadot*A_12 - (A_11 - 1)/tau_M )
            A_12 += dt*(  gammadot*A_22 - A_12/tau_M )
            
            sigma_elem = G * A_12 
            
            sigma.append(sigma_elem)
            time.append(t*dt)
            sum_G1 = sum_G1 + sigma_elem * np.sin(freq*t*dt*2*math.pi)*dt
            sum_G2 = sum_G2 + sigma_elem * np.cos(freq*t*dt*2*math.pi)*dt
        
        
        G1.append( sum_G1 * freq * 2*math.pi  /  strain_lin*math.pi)
        G2.append( sum_G2 * freq * 2*math.pi /  strain_lin*math.pi)

    return G1, G2, sigma, time

def Ruud_Single_Cycle_Full( G, tau_M, freq, strain, params, ncycles, eta_w,\
                            ntimesteps_percycle, gammadot_data, gammadot_howto):

    A_11 = 1 
    A_12 = 0
    A_22 = 1

    sigma_Y, nexp, alpha, gammadot_cr, I1c = params
    
    dt = 1/(freq*ntimesteps_percycle)
    G1 = []
    G2 = []

    sigma = []
    time = []
    for cycle in range(0,ncycles):

        sum_G1 = 0
        sum_G2 = 0 
        for t in range (1,ntimesteps_percycle+1):
     
            if gammadot_howto == 'from data':
                gammadot = gammadot_data[t]
            
            if gammadot_howto == 'from calculation':
                gammadot = strain * freq * 2*math.pi * np.cos(t*dt*freq*2*math.pi)

            #define the shearrate dependent viscosity
            eta_eff = ( sigma_Y/ np.abs(gammadot) ) * ( 1 + (np.abs(gammadot)/gammadot_cr)**nexp ) +eta_w
                    
            #define the damping function to use in strain dependent elastic modulus
            I1 = A_11 + A_22 -2 
            G_updated = G / ( 1 + np.abs(I1/I1c)**alpha)


            #Shear rate dependent viscosity + dtrain dependent elstic modulus = relaxation time
            tau_eff = (eta_eff)/(G_updated)


            #Finalize relaxation time by adding the solvent relaxation time
            itau = 1/tau_M + 1/tau_eff;  # Cf. Miyazaki ExtraMaxwell mode
            tau = 1.0/itau  

            #Solve for the structure tensor
            A_11 += dt*(2*gammadot*A_12 - (A_11 - 1)/tau )
            A_12 += dt*(  gammadot*A_22 - A_12/tau )

            #calcuate sigma by eqn 7 from Ruud's paper    
            sigma_elem = G_updated * A_12 
            
            sigma.append(sigma_elem)
            time.append(t*dt)
            sum_G1 = sum_G1 + sigma_elem * np.sin(freq*t*dt*2*math.pi)*dt
            sum_G2 = sum_G2 + sigma_elem * np.cos(freq*t*dt*2*math.pi)*dt
        
        
        G1.append( sum_G1 * (freq * 2*math.pi)  /  (strain*math.pi) )
        G2.append( sum_G2 * (freq * 2*math.pi)  /  (strain*math.pi) )

    return G1, G2, sigma, time

def Ruud_Steady( G, tau_M, params, eta_w, ntimesteps, gammadot):

    A_11 = 1 
    A_12 = 0
    A_22 = 1

    sigma_Y, nexp, alpha, gammadot_cr, I1c = params
    
    dt_list = np.logspace(-6,-3,int(ntimesteps))
    sigma = []
    time = []
    for t in range (0,int(ntimesteps)):
            
            dt = dt_list[t]

            #define the shearrate dependent viscosity
            eta_eff = ( sigma_Y/ np.abs(gammadot) ) * ( 1 + (np.abs(gammadot)/gammadot_cr)**nexp ) +eta_w
                    
            #define the damping function to use in strain dependent elastic modulus
            I1 = A_11 + A_22 -2 
            G_updated = G / ( 1 + np.abs(I1/I1c)**alpha)


            #Shear rate dependent viscosity + dtrain dependent elstic modulus = relaxation time
            tau_eff = (eta_eff)/(G_updated)


            #Finalize relaxation time by adding the solvent relaxation time
            itau = 1/tau_M + 1/tau_eff;  # Cf. Miyazaki ExtraMaxwell mode
            tau = 1.0/itau  

            #Solve for the structure tensor
            A_11 += dt*(2*gammadot*A_12 - (A_11 - 1)/tau )
            A_12 += dt*(  gammadot*A_22 - A_12/tau )

            #calcuate sigma by eqn 7 from Ruud's paper    
            sigma_elem = G_updated * A_12 
            
            sigma.append(sigma_elem)
            time.append(t*dt)

    return sigma, time

# ───────────────────────────────────────────────────────────────────
#   Optimization Functions
# ───────────────────────────────────────────────────────────────────
def optimization_lin_perfreq_ruud(initial_guess_lin, bounds_lin_norm, strain_val_lin, 
                current_freq, wholefreq, cycles,ntimesteps_cycle, col_index_stress, all_data, nelements_cycle ):   

        
    def objective_linear(params):

        ## Experimental ##

        #find the right index forthe freq sweep data
        extracted_ydata = []
        max_val = []
        l=0
        for freqval in wholefreq:
            if current_freq == freqval:
                index=l
            l+=1

        #extract the data
        last2cycles = all_data[0][0][col_index_stress][index*nelements_cycle:(index+1)*nelements_cycle] #waveform data at given frequency
        gammadotlist = all_data[0][0][3][index*nelements_cycle:(index+1)*nelements_cycle] #gammadot for inputtin into simulation
        for elem in last2cycles:
            extracted_ydata.append(elem)
        maxelem = max(last2cycles) 
        arraymaxelem = [maxelem]*nelements_cycle 
        for elem in arraymaxelem:
            max_val.append(elem)

        ## Predicted ##

        #Convert the normalized params to actual params
        params_actual = []
        for i in range(0,len(params)):
            param = params[i]
            bounds = bounds_lin_norm[i]
            param_actual1 = param * (bounds[1]-bounds[0]) + bounds[0]
            param_actual = 10**param_actual1
            params_actual.append(param_actual)

        #generate model data for given strains in linear region
        G1, G2, sigma, time = Ruud_Single_Cycle_Oldroydb(strain_val_lin,current_freq,\
                                                           params_actual, cycles, ntimesteps_cycle, gammadotlist)
        
        #read model data 
        pred_extracted_ydata = sigma[-nelements_cycle:]


        ## Difference ##        
        
        #take the difference
        diff = []
        i=0
        for elemdata in extracted_ydata:
            diff_elem = ( extracted_ydata[i] -  pred_extracted_ydata[i])
            
            diff_elem_norm = diff_elem / max_val[i] #Normalize with the max
            diff.append(diff_elem_norm)
            i+=1
        #take the square of difference
        diff_sqr = []
        for elemdiff in diff:
            diff_sqr.append(elemdiff**2)
        #sum up the squared values
        diff_sqr_sum = 0 
        for elemdiffsqr in diff_sqr:
                diff_sqr_sum  += elemdiffsqr 

        return diff_sqr_sum 


    #Minimize the objective function for linear region
    result = minimize(fun=objective_linear, x0=initial_guess_lin, method='Nelder-Mead') 
                     
    optimized_params_lin = result.x
    min_obj_lin = result.fun

    
    optim_actual = []
    for i in range(0,len(optimized_params_lin)):
            param = optimized_params_lin[i]
            bounds = bounds_lin_norm[i]
            param_actual1 = param * (bounds[1]-bounds[0]) + bounds[0]
            param_actual = 10**param_actual1
            optim_actual.append(param_actual)

    print(result)

    return optim_actual, min_obj_lin, result

def optimization_nonlin_together_ruud_HB_given(initial_guess_nonlin, HBparam, bounds_nonlin_norm,
                 strain_values_nonlin, freq,cycles,eta_inf,ntimesteps_cycle,\
                 col_index_stress, col_index_shearrate, all_data, optimized_params_lin, nelements_cycle):

    
    def objective_maxnorm(params):

        penalty =0 
        for paramval in params:
            if paramval < 0 or paramval > 1: 
                penalty= 1e8

        #calculate the params with units
        params_actual = [0]*5 #number of parameters is 5

        #input the HB params (TAU, n, dotgammacr)
        params_actual[0] = HBparam[0] #tau
        params_actual[1] = HBparam[1] #n
        params_actual[3] = HBparam[2] #dotgammacr

        #calculate the params with units for alpha and I1c
        for i in range(0,len(params)):
            param = params[i]
            bounds = bounds_nonlin_norm[i]
            param_actual = param * (bounds[1]-bounds[0]) + bounds[0]

            if i==0: 
                params_actual[2] = param_actual
            
            if i==1: #log for I1c
                param_actual = 10**(param_actual)
                params_actual[4] = param_actual

        
        #Extract predicted data
        #doing everything per frequency
        pred_extracted_ydata = [] #all data will be here
        freqindex=0
        for freqvals in freq:
            modulus = optimized_params_lin[0][0]*freqvals**optimized_params_lin[0][1]
            t_rel = optimized_params_lin[1][0]*freqvals**optimized_params_lin[1][1]
            
            strain_index=0
            for strainval in strain_values_nonlin:

                gammadot_data =all_data[strain_index+len(strain_values_nonlin)*freqindex][col_index_shearrate][-nelements_cycle:]
                G1, G2, sigma, time = Ruud_Single_Cycle_Full(modulus, t_rel, freqvals, strainval,\
                                                              params_actual, cycles, eta_inf,\
                                                              ntimesteps_cycle, gammadot_data, gammadot_howto='from data')
                
                for elem in sigma[-nelements_cycle:]:
                    pred_extracted_ydata.append(elem)
                
                strain_index+=1

            freqindex  += 1


        #Read Experimental Data
        extracted_ydata = []
        max_val = []
        for i_fq in range(0,len(freq)):
            for i_st in range(0,len(strain_values_nonlin)): 
                last2cycles = all_data[i_fq*len(strain_values_nonlin)+i_st][col_index_stress][-nelements_cycle:]
                for elem in last2cycles:
                    extracted_ydata.append(elem)

                maxelem = max(all_data[i_fq*len(strain_values_nonlin)+i_st][col_index_stress][-nelements_cycle:]) 
                arraymaxelem = [maxelem]*nelements_cycle 
                        
                for elem in arraymaxelem:
                    max_val.append(elem)

                                
        #calculate the difference btw Model and Experiments
        #difference between exp and model
        diff = []
        i=0
        for elemdata in extracted_ydata:
            diff_elem =  extracted_ydata[i]  - pred_extracted_ydata[i]

            diff_elem_norm = diff_elem / max_val[i]
            diff.append(diff_elem_norm**2)
            i+=1

        #sum up the diff values
        diff_sqr_sum = 0 
        for elemdiffsqr in diff:
                diff_sqr_sum  += elemdiffsqr 
    
        return diff_sqr_sum + penalty


    #Minimize the objective function for non-linear region
    result = minimize(fun=objective_maxnorm, x0=initial_guess_nonlin, method='Nelder-Mead') #, options={'xatol': 1e-4, 'fatol': 1e-4})  
    optimized_params_nonlin = result.x
    min_obj_nonlin = result.fun

    if result.success:
        print('Minimization is finished successfully for the non-linear region') 
    else: 
        print('Minimization is not successful')

    return optimized_params_nonlin, min_obj_nonlin, result

def optimization_nonlin_together_ruud(initial_guess_nonlin, bounds_nonlin_norm,
                 strain_values_nonlin, freq,cycles,eta_inf,ntimesteps_cycle,\
                 col_index_stress,col_index_shearrate, all_data, optimized_params_lin, nelements_cycle):

    
    def objective_maxnorm(params):

        penalty =0 
        for paramval in params:
            if paramval < 0 or paramval > 1: 
                penalty= 1e8

        params_actual = []
        for i in range(0,len(params)):
            param = params[i]
            bounds = bounds_nonlin_norm[i]
            param_actual = param * (bounds[1]-bounds[0]) + bounds[0]

            if i==0 or i == 3 or i == 4:
                param_actual = 10**(param_actual)
            params_actual.append(param_actual)
        
        #Extract predicted data
        #doing everything per frequency
        pred_extracted_ydata = [] #all data will be here
        freqindex=0
        for freqvals in freq:
            modulus = optimized_params_lin[0][0]*freqvals**optimized_params_lin[0][1]
            t_rel = optimized_params_lin[1][0]*freqvals**optimized_params_lin[1][1]
            
            strain_index=0
            for strainval in strain_values_nonlin:

                gammadot_data =all_data[strain_index+len(strain_values_nonlin)*freqindex][col_index_shearrate][-nelements_cycle:]
                G1, G2, sigma, time = Ruud_Single_Cycle_Full(modulus, t_rel, freqvals, strainval,\
                                                              params_actual, cycles, eta_inf,\
                                                              ntimesteps_cycle, gammadot_data, gammadot_howto='from data')
                
                for elem in sigma[-nelements_cycle:]:
                    pred_extracted_ydata.append(elem)
                
                strain_index+=1

            freqindex  += 1


        #Read Experimental Data
        extracted_ydata = []
        max_val = []
        for i_fq in range(0,len(freq)):
            for i_st in range(0,len(strain_values_nonlin)): 
                last2cycles = all_data[i_fq*len(strain_values_nonlin)+i_st][col_index_stress][-nelements_cycle:]
                for elem in last2cycles:
                    extracted_ydata.append(elem)

                maxelem = max(all_data[i_fq*len(strain_values_nonlin)+i_st][col_index_stress][-nelements_cycle:]) 
                arraymaxelem = [maxelem]*nelements_cycle #spesific
                        
                for elem in arraymaxelem:
                    max_val.append(elem)

                                
        #calculate the difference btw Model and Experiments
        #difference between exp and model
        diff = []
        i=0
        for elemdata in extracted_ydata:
            diff_elem =  extracted_ydata[i]  - pred_extracted_ydata[i]

            diff_elem_norm = diff_elem / max_val[i]
            diff.append(diff_elem_norm**2)
            i+=1

        #sum up the diff values
        diff_sqr_sum = 0 
        for elemdiffsqr in diff:
                diff_sqr_sum  += elemdiffsqr 
    
        return diff_sqr_sum + penalty


    #Minimize the objective function for non-linear region
    result = minimize(fun=objective_maxnorm, x0=initial_guess_nonlin, method='Nelder-Mead') #, options={'xatol': 1e-4, 'fatol': 1e-4})  
    optimized_params_nonlin = result.x
    min_obj_nonlin = result.fun

    if result.success:
        print('Minimization is finished successfully for the non-linear region') 
    else: 
        print('Minimization is not successful')

    return optimized_params_nonlin, min_obj_nonlin, result

# ───────────────────────────────────────────────────────────────────
#   Sum of squared for error calcualtions
# ───────────────────────────────────────────────────────────────────
def sum_of_squared_residuals_ruud(params_actual, current_freq, strain_val_lin,\
                                     wholefreq, cycles,ntimesteps_cycle, col_index_stress,col_index_shearrate,\
                                     all_data, nelements_cycle):
    
        ## Experimental ##

        #find the right index forthe freq sweep data
        extracted_ydata = []
        max_val = []
        l=0
        for freqval in wholefreq:
            if current_freq == freqval:
                index=l
            l+=1

        #extract the data
        last2cycles = all_data[0][0][col_index_stress][index*nelements_cycle:(index+1)*nelements_cycle] #waveform data at given frequency
        gammadotlist = all_data[0][0][col_index_shearrate][index*nelements_cycle:(index+1)*nelements_cycle] #gammadot for inputtin into simulation
        for elem in last2cycles:
            extracted_ydata.append(elem)
        maxelem = max(last2cycles) 
        arraymaxelem = [maxelem]*nelements_cycle 
        for elem in arraymaxelem:
            max_val.append(elem)

        ## Predicted ##

        #generate model data for given strains in linear region
        G1, G2, sigma, time = Ruud_Single_Cycle_Oldroydb(strain_val_lin,current_freq,\
                                                           params_actual, cycles, ntimesteps_cycle, gammadotlist)
        
        #read model data 
        pred_extracted_ydata = sigma[-nelements_cycle:]


        ## Difference ##        
        
        #take the difference
        diff = []
        i=0
        for elemdata in extracted_ydata:
            diff_elem = ( extracted_ydata[i] -  pred_extracted_ydata[i])
            
            diff_elem_norm = diff_elem / max_val[i] #Normalize with the max
            diff.append(diff_elem_norm)
            i+=1
        #take the square of difference
        diff_sqr = []
        for elemdiff in diff:
            diff_sqr.append(elemdiff**2)
        #sum up the squared values
        diff_sqr_sum = 0 
        for elemdiffsqr in diff_sqr:
                diff_sqr_sum  += elemdiffsqr 

        return diff_sqr_sum

def sum_of_squared_residuals_laos_ruud(params_actual, strain_values_nonlin, freq, cycles, eta_inf, ntimesteps_cycle,\
                                   optimized_params_lin, col_index_stress, col_index_shearrate, all_data, nelements_cycle ):
        
         #Extract predicted data
        #doing everything per frequency
        pred_extracted_ydata = [] #all data will be here
        freqindex=0
        for freqvals in freq:
            modulus = optimized_params_lin[0][0]*freqvals**optimized_params_lin[0][1]
            t_rel = optimized_params_lin[1][0]*freqvals**optimized_params_lin[1][1]
            
            strain_index=0
            for strainval in strain_values_nonlin:

                gammadot_data =all_data[strain_index+len(strain_values_nonlin)*freqindex][col_index_shearrate][-nelements_cycle:]
                G1, G2, sigma, time = Ruud_Single_Cycle_Full(modulus, t_rel, freqvals, strainval,\
                                                              params_actual, cycles, eta_inf,\
                                                              ntimesteps_cycle, gammadot_data, gammadot_howto='from data')
                
                for elem in sigma[-nelements_cycle:]:
                    pred_extracted_ydata.append(elem)
                
                strain_index+=1

            freqindex  += 1


        #Read Experimental Data
        extracted_ydata = []
        max_val = []
        for i_fq in range(0,len(freq)):
            for i_st in range(0,len(strain_values_nonlin)): 
                last2cycles = all_data[i_fq*len(strain_values_nonlin)+i_st][col_index_stress][-nelements_cycle:]
                for elem in last2cycles:
                    extracted_ydata.append(elem)

                maxelem = max(all_data[i_fq*len(strain_values_nonlin)+i_st][col_index_stress][-nelements_cycle:]) 
                arraymaxelem = [maxelem]*nelements_cycle #spesific
                        
                for elem in arraymaxelem:
                    max_val.append(elem)
                                
        #calculate the difference btw Model and Experiments
        #difference between exp and model
        diff = []
        i=0
        for elemdata in extracted_ydata:
            diff_elem =  extracted_ydata[i]  - pred_extracted_ydata[i]

            diff_elem_norm = diff_elem / max_val[i]
            diff.append(diff_elem_norm**2)
            i+=1

        #sum up the diff values
        diff_sqr_sum = 0 
        for elemdiffsqr in diff:
                diff_sqr_sum  += elemdiffsqr 
    
        return diff_sqr_sum 

