
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

def read_data(n_header, n_rows, header_list, what_to_collect_time, what_to_collect_single,path_read, sheet_name, sample_list, rows_btw_data):
    
    data = []
    data_single = []
    i=0
    for sample in sample_list:
        sample_strain_single = []
        sample_strain = []

        ##print(sample)
        var = pd.read_excel(io=path_read, sheet_name=sheet_name[0], header=n_header, names=header_list, nrows=n_rows)
        

        for elem1 in what_to_collect_time:
            
            elem = var[elem1].tolist()
            sample_strain.append(elem)

        for elem2 in what_to_collect_single:
            a= var[elem2].tolist()
            #print(a)

            
            for k in range(0,5):
                sample_strain_single.append(a[k*257])
                
        n_header = n_header + n_rows + rows_btw_data

        data.append(sample_strain)
        data_single.append(sample_strain_single)
        i+=1

    return data, data_single

def read_data1(n_header, n_rows, header_list, what_to_collect_time, what_to_collect_single,path_read, sheet_name, sample_list, rows_btw_data):
    

    data = []
    data_single = []
    i=0
    for sample in sample_list:
        sample_strain = []

        ##print(sample)
        var = pd.read_excel(io=path_read, sheet_name=sheet_name[0], header=n_header, names=header_list, nrows=n_rows)
        
        # df = pd.DataFrame(var)

        # #print(df.colums)
        ##print(var)

        for elem1 in what_to_collect_time:
            sample_strain.append(var[elem1].tolist())
                
        n_header = n_header + n_rows + rows_btw_data

        data.append(sample_strain)
        i+=1

    return data

def fitting_func_SAOS_freq_and_strain(params, strain,\
                      freq,cycles,ntimesteps_cycle):
    #takes strain as unitless not percnetage!!
    generate_input=1
    run_the_cases=1
    read_the_cases=1

    #print(params)

    # Part 1: generate input 
    #here the gammadot is not taken from the data
    if generate_input == 1:
    
        os.system('rm -r inputs_lin')
        os.system('rm -r results_lin')
        os.system('mkdir -p inputs_lin')
        os.system('mkdir -p results_lin')
        os.system('mkdir results_lin/terminal_output')
        os.system('mkdir results_lin/temporal_output')

        modulus,t_rel = params
        #modify according to the frequency
        ii = 1
        f = open('inputs_lin/input{:04}.dat'.format(ii), 'w')       #write the input file
        f.write('&comppar\n\
        modulus = {}       \n\
        lambda = {}        \n\
        gamma_0 = {}       \n\
        freq = {}          \n\
        cycles = {}        \n\
        numtimesteps_percycle = {}  \n\
        filename = "fitout_{:04}"   \n\
    /\n\n\n\n'.format(modulus,t_rel, strain, freq, cycles, ntimesteps_cycle,ii))
        f.close
            
                
        #generate one extra set of input for the last one since very last iput file has an issue with reading
                    
        f = open('inputs_lin/input{:04}.dat'.format(ii), 'w')       #write the input file
        f.write('&comppar\n\
    modulus = {}       \n\
    lambda = {}        \n\
    gamma_0 = {}       \n\
    freq = {}          \n\
    cycles = {}        \n\
    numtimesteps_percycle = {}  \n\
    filename = "fitout_{:04}"   \n\
/\n\n\n\n'.format(modulus,t_rel, strain, freq, cycles, ntimesteps_cycle,ii))
        f.close

    # Part 2: run the cases   
    if run_the_cases == 1:
        
            os.system('touch strain_sweep_lin.txt')
            #os.system('make')
            os.system("./startup_fitting_linear < inputs_lin/input{:04}.dat > results_lin/terminal_output/{:04}".format(ii,ii))
            os.system('mv fitout_* results_lin/temporal_output')
            os.system('mv strain_sweep_lin.txt results_lin/strain_sweep_{:04}.txt'.format(ii))


    # Part 3: read the cases
    if read_the_cases == 1:

        col_names=['freq', 'G1', 'G2']
        
        sampleG1=[]   
        sampleG2=[]
        path = f'results_lin/strain_sweep_{1:04}.txt'
        var1 = pd.read_fwf(filepath_or_buffer=path , header=None, names=col_names,\
                                widths=[12,12,12])
                
        G1 = var1['G1'].tolist()
        G2 = var1['G2'].tolist()
            
        sampleG1.append(G1[-1])
        sampleG2.append(G2[-1])
            
    return sampleG1, sampleG2

def fitting_func_LAOS(params, strain_values,\
                      freq,cycles,eta_inf,ntimesteps_cycle, modulus, t_rel, data_exp, freqindex):
    #takes strain as unitless not percnetage!!
    generate_input=1
    run_the_cases=1
    read_the_cases=1

    # Part 1: generate input 
    if generate_input == 1:
    
        os.system('rm -r inputs')
        os.system('rm -r results')
        os.system('mkdir -p inputs')
        os.system('mkdir -p results')
        os.system('mkdir results/terminal_output')
        os.system('mkdir results/temporal_output')

        #print('calculating for params:{}'.format(params))
        tau_y,nexp,alpha,gammadot_cr,I1c = params

        ii=0
        for strain in strain_values:
                
                gammadotlist = data_exp[ii+len(strain_values)*freqindex][3][-257:]
                gammadotlist_check = replace_nan_with_avg(gammadotlist)
                gdlist =  ', '.join(map(str, np.array(gammadotlist_check)))
                

                f = open('inputs/input{:04}.dat'.format(ii), 'w')       #write the input file
                f.write('&comppar\n\
        modulus = {}       \n\
        lambda = {}        \n\
        tau_y = {}         \n\
        nexp = {}          \n\
        alpha = {}         \n\
        gammadot_cr = {}   \n\
        eta_inf = {}       \n\
        I1c = {}           \n\
        gamma_0 = {}       \n\
        freq = {}          \n\
        cycles = {}        \n\
        gammadotlist = {}  \n\
        numtimesteps_percycle = {}  \n\
        filename = "fitout_{:04}"   \n\
/\n\n\n\n'.format(modulus,t_rel,tau_y,nexp, alpha, gammadot_cr, eta_inf, I1c, strain,\
                        freq, cycles, gdlist, ntimesteps_cycle,ii))
                f.close
                ii+=1
    
                
        #generate one extra set of input for the last one since very last iput file has an issue with reading
                    
        f = open('inputs/input{:04}.dat'.format(ii), 'w')       #write the input file
        f.write('&comppar\n\
    modulus = {}       \n\
    lambda = {}        \n\
    tau_y = {}         \n\
    nexp = {}          \n\
    alpha = {}         \n\
    gammadot_cr = {}   \n\
    eta_inf = {}       \n\
    I1c = {}           \n\
    gamma_0 = {}       \n\
    freq = {}          \n\
    cycles = {}        \n\
    gammadotlist = {}  \n\
    numtimesteps_percycle = {}  \n\
    filename = "fitout_{:04}"   \n\
/\n\n\n\n'.format(modulus,t_rel,tau_y,nexp, alpha, gammadot_cr, eta_inf, I1c, strain,\
            freq, cycles, gdlist, ntimesteps_cycle,ii))
        f.close

    # Part 2: run the cases   
    if run_the_cases == 1:
        
        
        ii=0
        for strain in strain_values:

                os.system('touch strain_sweep.txt')
                os.system("./startup_fitting_datainput < inputs/input{:04}.dat > results/terminal_output/{:04}".format(ii,ii))
                os.system('mv fitout_* results/temporal_output')
                os.system('mv strain_sweep.txt results/strain_sweep_{:04}.txt'.format(ii))

                ii+=1

    # Part 3: read the cases
    if read_the_cases == 1:

        col_names=['cycle', 'pertime', 'tauxx', 'tauxy', 'tauyy',\
                    'vonmis','gammadot', 'gammadot_eff', 'eta_eff', 'G_mod', 'I1']
        

        res_list=sorted(os.listdir('results/temporal_output'))
        sample=[]
            
        for res in res_list:

            path = 'results/temporal_output/{}'.format(res)
            var1 = pd.read_fwf(filepath_or_buffer=path , header=None, names=col_names,\
                                widths=[12,12,12,12,12,12,12,12,12,12,12])
                
            time = var1['pertime'].tolist()
            xy = var1['tauxy'].tolist()

                
            element = [time,xy]
            sample.append(element)
            
    return sample

def fitting_func_LAOS_nodatainput(params, strain_values,\
                      freq,cycles,eta_inf,ntimesteps_cycle, modulus, t_rel):
    #takes strain as unitless not percnetage!!
    generate_input=1
    run_the_cases=1
    read_the_cases=1

    # Part 1: generate input 
    if generate_input == 1:
    
        os.system('rm -r inputs')
        os.system('rm -r results')
        os.system('mkdir -p inputs')
        os.system('mkdir -p results')
        os.system('mkdir results/terminal_output')
        os.system('mkdir results/temporal_output')

        #print('calculating for params:{}'.format(params))
        tau_y,nexp,alpha,gammadot_cr,I1c = params

        ii=0
        for strain in strain_values:                

                f = open('inputs/input{:04}.dat'.format(ii), 'w')       #write the input file
                f.write('&comppar\n\
        modulus = {}       \n\
        lambda = {}        \n\
        tau_y = {}         \n\
        nexp = {}          \n\
        alpha = {}         \n\
        gammadot_cr = {}   \n\
        eta_inf = {}       \n\
        I1c = {}           \n\
        gamma_0 = {}       \n\
        freq = {}          \n\
        cycles = {}        \n\
        numtimesteps_percycle = {}  \n\
        filename = "fitout_{:04}"   \n\
/\n\n\n\n'.format(modulus,t_rel,tau_y,nexp, alpha, gammadot_cr, eta_inf, I1c, strain,\
                        freq, cycles,  ntimesteps_cycle,ii))
                f.close
                ii+=1
    
                
        #generate one extra set of input for the last one since very last iput file has an issue with reading
                    
        f = open('inputs/input{:04}.dat'.format(ii), 'w')       #write the input file
        f.write('&comppar\n\
    modulus = {}       \n\
    lambda = {}        \n\
    tau_y = {}         \n\
    nexp = {}          \n\
    alpha = {}         \n\
    gammadot_cr = {}   \n\
    eta_inf = {}       \n\
    I1c = {}           \n\
    gamma_0 = {}       \n\
    freq = {}          \n\
    cycles = {}        \n\
    numtimesteps_percycle = {}  \n\
    filename = "fitout_{:04}"   \n\
/\n\n\n\n'.format(modulus,t_rel,tau_y,nexp, alpha, gammadot_cr, eta_inf, I1c, strain,\
            freq, cycles, ntimesteps_cycle,ii))
        f.close

    # Part 2: run the cases   
    if run_the_cases == 1:
        
        
        ii=0
        for strain in strain_values:
                os.system('touch strain_sweep.txt')
                os.system("./startup_fitting < inputs/input{:04}.dat > results/terminal_output/{:04}".format(ii,ii))
                os.system('mv fitout_* results/temporal_output')
                os.system('mv strain_sweep.txt results/strain_sweep_{:04}.txt'.format(ii))

                ii+=1

    # Part 3: read the cases
    if read_the_cases == 1:

        col_names=['cycle', 'pertime', 'tauxx', 'tauxy', 'tauyy',\
                    'vonmis','gammadot', 'gammadot_eff', 'eta_eff', 'G_mod', 'I1']
        

        res_list=sorted(os.listdir('results/temporal_output'))
        sample=[]
            
        for res in res_list:

            path = 'results/temporal_output/{}'.format(res)
            var1 = pd.read_fwf(filepath_or_buffer=path , header=None, names=col_names,\
                                widths=[12,12,12,12,12,12,12,12,12,12,12])
                
            time = var1['pertime'].tolist()
            xy = var1['tauxy'].tolist()

                
            element = [time,xy]
            sample.append(element)
            
    return sample

def optimization_lin_modulus(initial_guess_lin, bounds_lin_norm, bounds_lin, strain_lin, 
                freq, freqvals, cycles,ntimesteps_cycle, data_exp_frswp ):   

        
    def objective_linear(params):
        
        #Convert the normalized params to actual params
        params_actual = []
        for i in range(0,len(params)):
            param = params[i]
            bounds = bounds_lin_norm[i]
            param_actual1 = param * (bounds[1]-bounds[0]) + bounds[0]
            param_actual = 10**param_actual1
            params_actual.append(param_actual)

        G1pred = []
        G2pred = []
        for strain in strain_lin:

            #Generate Model Data for given strains in linear region
            G1pred_elem, G2pred_elem = fitting_func_SAOS_freq_and_strain(params_actual, strain=strain,\
                            freq=freq,cycles=cycles,\
                            ntimesteps_cycle=ntimesteps_cycle)
            G1pred.append(G1pred_elem)
            G2pred.append(G2pred_elem)

        

        G1exp = data_exp_frswp[0][1] #there is only one elemnt in the set
        G2exp = data_exp_frswp[0][2]

        if len(freqvals) > 1: #freq sweep data
            
            #get the index for the frequency
            l=0
            for freqval in freqvals:
                if freqval == freq:
                    index=l
                l+=1

            #difference between exp and model
            sumsqdiff = []
            for k in range(0,len(G1pred)):
                diffG1 = ( (G1exp[index]- G1pred[k]) / G1exp[index] )**2
                diffG2 = ( (G2exp[index]- G2pred[k]) / G2exp[index] )**2

        if len(freqvals) == 1 and len(strain_lin) > 1: #ampl sweep data
            #difference between exp and model
            sumsqdiff = []
            for k in range(0,len(G1pred)):
                diffG1 = ( (G1exp[k]- G1pred[k]) / G1exp[k] )**2
                diffG2 = ( (G2exp[k]- G2pred[k]) / G2exp[k] )**2
             

        sumsqdiff.append(diffG1)
        sumsqdiff.append(diffG2)

        meansumsqdiff = np.mean(sumsqdiff)
        
        return meansumsqdiff
    
    os.system('rm strain_sweep_lin.txt')
    #Minimize the objective function for linear region
    result = minimize(fun=objective_linear, x0=initial_guess_lin, method='Nelder-Mead', bounds=bounds_lin)          
    optimized_params_lin = result.x
    min_obj_lin = result.fun

    
    optim_actual = []
    for i in range(0,len(optimized_params_lin)):
            param = optimized_params_lin[i]
            bounds = bounds_lin_norm[i]
            param_actual1 = param * (bounds[1]-bounds[0]) + bounds[0]
            param_actual = 10**param_actual1
            optim_actual.append(param_actual)


    return optim_actual, min_obj_lin, result

def fitting_func_SAOS_perfreq(params, strain,\
                      freq_vals,cycles,ntimesteps_cycle, gammadotlist):
    #takes strain as unitless not percnetage!!
    generate_input=1
    run_the_cases=1
    read_the_cases=1

    # Part 1: generate input 
    if generate_input == 1:
    
        os.system('rm -r inputs_lin')
        os.system('rm -r results_lin')
        os.system('mkdir -p inputs_lin')
        os.system('mkdir -p results_lin')
        os.system('mkdir results_lin/terminal_output')
        os.system('mkdir results_lin/temporal_output')

        #print('calculating for params:{}'.format(params))
        modulus,t_rel = params

        gdlist =  ', '.join(map(str, np.array(gammadotlist)))
       
  
                
        ii=0
        for freq in freq_vals:

                f = open('inputs_lin/input{:04}.dat'.format(ii), 'w')       #write the input file
                f.write('&comppar\n\
        modulus = {}       \n\
        lambda = {}        \n\
        gamma_0 = {}       \n\
        freq = {}          \n\
        cycles = {}        \n\
        gammadotlist = {}        \n\
        numtimesteps_percycle = {}  \n\
        filename = "fitout_{:04}"   \n\
    /\n\n\n\n'.format(modulus,t_rel, strain, freq, cycles, gdlist, ntimesteps_cycle,ii))
                f.close
                ii+=1
            
                
        #generate one extra set of input for the last one since very last iput file has an issue with reading
                    
        f = open('inputs_lin/input{:04}.dat'.format(ii), 'w')       #write the input file
        f.write('&comppar\n\
    modulus = {}       \n\
    lambda = {}        \n\
    gamma_0 = {}       \n\
    freq = {}          \n\
    cycles = {}        \n\
    gammadotlist = {}        \n\
    numtimesteps_percycle = {}  \n\
    filename = "fitout_{:04}"   \n\
/\n\n\n\n'.format(modulus,t_rel, strain, freq, cycles, gammadotlist, ntimesteps_cycle,ii))
        f.close

    # Part 2: run the cases   
    if run_the_cases == 1:
        
        
        ii=0
        for freq in freq_vals:

                os.system('touch strain_sweep_lin.txt')
                #os.system('make')
                os.system("./startup_fitting_linear_datainput < inputs_lin/input{:04}.dat > results_lin/terminal_output/{:04}".format(ii,ii))
                os.system('mv fitout_* results_lin/temporal_output')
                os.system('mv strain_sweep_lin.txt results_lin/strain_sweep_{:04}.txt'.format(ii))

                ii+=1

    # Part 3: read the cases
    if read_the_cases == 1:

        path_list=[]
        col_names=['cycle', 'pertime', 'tauxx', 'tauxy', 'tauyy',\
                    'vonmis','gammadot', 'gammadot_eff']
        

        res_list=sorted(os.listdir('results_lin/temporal_output'))
        sample=[]
            
        for res in res_list:
            sample_interval=[]

            path = 'results_lin/temporal_output/{}'.format(res)
            var1 = pd.read_fwf(filepath_or_buffer=path , header=None, names=col_names,\
                                widths=[12,12,12,12,12,12,12,12])
                            
            time = var1['pertime'].tolist()
            xy = var1['tauxy'].tolist()
            
            element = [time,xy]
            sample.append(element)
            
    return sample

def optimization_lin_perfreq(initial_guess_lin, bounds_lin_norm, strain_val_lin, 
                freq_vals, wholefreq, cycles,ntimesteps_cycle, col_index_stress, all_data ):   

        
    def objective_linear(params):

        ## Experimental ##

        #find the right index forthe freq sweep data
        extracted_ydata = []
        max_val = []
        l=0
        for freqval in wholefreq:
            if freq_vals[0] == freqval:
                index=l
            l+=1

        #extract the data
        last2cycles = all_data[0][0][col_index_stress][index*257:(index+1)*257] #waveform data at given frequency
        gammadotlist = all_data[0][0][3][index*257:(index+1)*257] #gammadot for inputtin into simulation
        for elem in last2cycles:
            extracted_ydata.append(elem)
        maxelem = max(last2cycles) 
        arraymaxelem = [maxelem]*257 
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
        predicted_whole_data = fitting_func_SAOS_perfreq(params_actual,strain_val_lin,\
                        freq_vals=freq_vals,cycles=cycles,\
                        ntimesteps_cycle=ntimesteps_cycle, gammadotlist=gammadotlist)
        
        #read model data 
        pred_extracted_ydata = []
        for m in range(0,len(freq_vals)):
            last2cycles = predicted_whole_data[m][col_index_stress][-257:]
            for elem in last2cycles:
                pred_extracted_ydata.append(elem)


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


    os.system('rm strain_sweep_lin.txt')
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

def optimization_nonlin(initial_guess_nonlin, bounds_nonlin_norm,bounds_nonlin,\
                 strain_values_nonlin, freq,cycles,eta_inf,ntimesteps_cycle,\
                 col_index_stress, all_data, optimized_params_lin, objfunc, alg, options):

    modulus = optimized_params_lin[0]
    t_rel = optimized_params_lin[1]

    def objective_maxnorm(params):

        #calculate the params with units
        params_actual = []
        for i in range(0,len(params)):
            param = params[i]
            bounds = bounds_nonlin_norm[i]
            param_actual = param * (bounds[1]-bounds[0]) + bounds[0]

            if i==0 or i == 3 or i == 4:
                param_actual = 10**(param_actual)
            params_actual.append(param_actual)


        #Generate Model Data
        predicted_whole_data = fitting_func_LAOS(params_actual,strain_values=strain_values_nonlin,\
                        freq=freq,cycles=cycles,eta_inf=eta_inf,\
                        ntimesteps_cycle=ntimesteps_cycle, modulus=modulus,t_rel=t_rel)
        
        #Read Model Data 
        pred_extracted_ydata = []
        for strain in predicted_whole_data:
            last2cycles = strain[col_index_stress][-257*4:-257*3]
            for elem in last2cycles:
                pred_extracted_ydata.append(elem)

        #Read Experimental Data
        extracted_ydata = []
        max_val = []
        for i in range(0,len(strain_values_nonlin)):
                
                strain = all_data[i] 
                last2cycles = strain[col_index_stress][-257*4:-257*3]
                
                for elem in last2cycles:
                    extracted_ydata.append(elem)

                maxelem = max(last2cycles) 
                arraymaxelem = [maxelem]*257 #spesific
                    
                for elem in arraymaxelem:
                    max_val.append(elem)

                                
        #difference between exp and model
        diff = []
        i=0
        for elemdata in extracted_ydata:
            diff_elem =  abs( extracted_ydata[i]  - pred_extracted_ydata[i] )

            diff_elem_norm = diff_elem / max_val[i]
            diff.append(diff_elem_norm**2)
            i+=1
       
        #sum up the diff values
        diff_sqr_sum = 0 
        for elemdiffsqr in diff:
                diff_sqr_sum  += elemdiffsqr / len(extracted_ydata)
       
        diffsumsq_list_nonlin.append(diff_sqr_sum)

        return diff_sqr_sum

    if objfunc == 'MAX':
        function = objective_maxnorm

    os.system('rm strain_sweep.txt')
    
    #Minimize the objective function for non-linear region
    diffsumsq_list_nonlin = []
    options = options
    result = minimize(fun=function, x0=initial_guess_nonlin, method=alg,\
                      bounds=bounds_nonlin, options=options) #, jac='2-point') 
    num_iterations_nonlin = result.nit
    optimized_params_nonlin = result.x
    min_obj_nonlin = result.fun
    jac = approx_fprime(result.x, objective_maxnorm, options['eps'])
    
    return optimized_params_nonlin, min_obj_nonlin, result, jac

def optimization_nonlin_HBgiven(initial_guess_nonlin, HBparam, bounds_nonlin_norm,\
                 strain_values_nonlin, freq,cycles,eta_inf,ntimesteps_cycle,\
                 col_index_stress, all_data, optimized_params_lin):


    def objective_maxnorm(params): #it is only alpha, I1c
        
        penalty =0 
        for paramval in params:
            if paramval < 0 or paramval > 1: 
                print(params)
                print('penalty applied')
                penalty= 1e8

        #calculate the params with units
        params_actual = [0]*5

        #input the HB params (TAU, dotgammacr, n)
        params_actual[0] = HBparam[0]
        params_actual[1] = HBparam[1]
        params_actual[3] = HBparam[2]

        #calculate the params with units for alpha and I1c
        for i in range(0,len(params)):
            param = params[i]
            bounds = bounds_nonlin_norm[i]
            param_actual = param * (bounds[1]-bounds[0]) + bounds[0]

            if i==0:
                params_actual[2] = param_actual
            
            if i==1:
                param_actual = 10**(param_actual)
                params_actual[4] = param_actual


        #doing everything per frequency
        pred_extracted_ydata = [] #all data will be here
        for freqvals in freq:
            modulus = optimized_params_lin[0][0]*freqvals**optimized_params_lin[0][1]
            t_rel = optimized_params_lin[1][0]*freqvals**optimized_params_lin[1][1]

            if freqvals == freq[0] : freqindex = 0
            if freqvals == freq[1] : freqindex = 1
            if freqvals == freq[2] : freqindex = 2

            #Generate Model Data for all strains of a freq
            predicted_whole_data_freq = fitting_func_LAOS(params_actual,strain_values=strain_values_nonlin,\
                            freq=freqvals,cycles=cycles,eta_inf=eta_inf,\
                            ntimesteps_cycle=ntimesteps_cycle, modulus=modulus,t_rel=t_rel,\
                            data_exp=all_data, freqindex=freqindex)
            
            #Read Model Data 
            pred_extracted_ydata_freq = []
            for strain in predicted_whole_data_freq:
                last2cycles = strain[col_index_stress][-257:]

                i=0
                for elem in last2cycles:
                    pred_extracted_ydata_freq.append(elem)
            
            for elem in pred_extracted_ydata_freq:
                pred_extracted_ydata.append(elem)


        #Read Experimental Data
        extracted_ydata = []
        max_val = []
        for i_fq in range(0,len(freq)):
            for i_st in range(0,len(strain_values_nonlin)): 
                last2cycles = all_data[i_fq*7+i_st][1][-257:]
                for elem in last2cycles:
                    extracted_ydata.append(elem)

                maxelem = max(all_data[i_fq*7+i_st][1][-257:]) 
                arraymaxelem = [maxelem]*257 #spesific
                        
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
    result = minimize(fun=objective_maxnorm, x0=initial_guess_nonlin, method='Nelder-Mead') 
    optimized_params_nonlin = result.x
    min_obj_nonlin = result.fun
    
    return optimized_params_nonlin, min_obj_nonlin, result

def optimization_nonlin_together(initial_guess_nonlin, bounds_nonlin_norm,
                 strain_values_nonlin, freq,cycles,eta_inf,ntimesteps_cycle,\
                 col_index_stress, all_data, optimized_params_lin):

    
    def objective_maxnorm(params):

        #print(params)

        penalty =0 
        for paramval in params:
            if paramval < 0 or paramval > 1: 
                #print('penalty applied')
                penalty= 1e8

        params_actual = []
        for i in range(0,len(params)):
            param = params[i]
            bounds = bounds_nonlin_norm[i]
            param_actual = param * (bounds[1]-bounds[0]) + bounds[0]

            if i==0 or i == 3 or i == 4:
                param_actual = 10**(param_actual)
            params_actual.append(param_actual)
        
        #doing everything per frequency
        pred_extracted_ydata = [] #all data will be here
        for freqvals in freq:
            modulus = optimized_params_lin[0][0]*freqvals**optimized_params_lin[0][1]
            t_rel = optimized_params_lin[1][0]*freqvals**optimized_params_lin[1][1]

            #print(f'freq:{freqvals} modulus:{modulus} trel:{t_rel}')

            if freqvals == freq[0] : freqindex = 0
            if freqvals == freq[1] : freqindex = 1
            if freqvals == freq[2] : freqindex = 2

            #Generate Model Data for all strains of a freq
            predicted_whole_data_freq = fitting_func_LAOS(params_actual,strain_values=strain_values_nonlin,\
                            freq=freqvals,cycles=cycles,eta_inf=eta_inf,\
                            ntimesteps_cycle=ntimesteps_cycle, modulus=modulus,t_rel=t_rel,\
                            data_exp=all_data, freqindex=freqindex)
            
            #Read Model Data 
            pred_extracted_ydata_freq = []
            for strain in predicted_whole_data_freq:
                last2cycles = strain[col_index_stress][-257:]

                i=0
                for elem in last2cycles:
                    pred_extracted_ydata_freq.append(elem)
            
            for elem in pred_extracted_ydata_freq:
                pred_extracted_ydata.append(elem)


        #Read Experimental Data
        extracted_ydata = []
        max_val = []
        for i_fq in range(0,len(freq)):
            for i_st in range(0,len(strain_values_nonlin)): 
                last2cycles = all_data[i_fq*len(strain_values_nonlin)+i_st][1][-257:]
                for elem in last2cycles:
                    extracted_ydata.append(elem)

                maxelem = max(all_data[i_fq*len(strain_values_nonlin)+i_st][1][-257:]) 
                arraymaxelem = [maxelem]*257 #spesific
                        
                for elem in arraymaxelem:
                    max_val.append(elem)


        # print(f'pred:{len(pred_extracted_ydata)}')
        # print(f'exp:{len(extracted_ydata)}')
                                
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


    os.system('rm strain_sweep.txt')
    #Minimize the objective function for non-linear region
    result = minimize(fun=objective_maxnorm, x0=initial_guess_nonlin, method='Nelder-Mead') #, options={'xatol': 1e-4, 'fatol': 1e-4})  
    optimized_params_nonlin = result.x
    min_obj_nonlin = result.fun

    if result.success:
        print('Minimization is finished successfully for the non-linear region') 
    else: 
        print('Minimization is not successful')

    return optimized_params_nonlin, min_obj_nonlin, result

def optimization_nonlin_individualfreq(initial_guess_nonlin, bounds_nonlin_norm,
                 strain_values_nonlin, freq,cycles,eta_inf,ntimesteps_cycle,\
                 col_index_stress, all_data, optimized_params_lin, freqselect):

    
    def objective_maxnorm(params):

        penalty =0 
        for paramval in params:
            if paramval < 0 or paramval > 1: 
                print('penalty applied')
                penalty= 1e8

        params_actual = []
        for i in range(0,len(params)):
            param = params[i]
            bounds = bounds_nonlin_norm[i]
            param_actual = param * (bounds[1]-bounds[0]) + bounds[0]

            if i==0 or i == 3 or i == 4:
                param_actual = 10**(param_actual)
            params_actual.append(param_actual)
        
        #doing everything per frequency
        pred_extracted_ydata = [] #all data will be here
        for freqvals in freqselect:
            modulus = optimized_params_lin[0][0]*freqvals**optimized_params_lin[0][1]
            t_rel = optimized_params_lin[1][0]*freqvals**optimized_params_lin[1][1]

            if freqvals == freq[0] : freqindex = 0
            if freqvals == freq[1] : freqindex = 1
            if freqvals == freq[2] : freqindex = 2

            #Generate Model Data for all strains of a freq
            predicted_whole_data_freq = fitting_func_LAOS(params_actual,strain_values=strain_values_nonlin,\
                            freq=freqvals,cycles=cycles,eta_inf=eta_inf,\
                            ntimesteps_cycle=ntimesteps_cycle, modulus=modulus,t_rel=t_rel,\
                            data_exp=all_data, freqindex=freqindex)
            
            #Read Model Data 
            pred_extracted_ydata_freq = []
            for strain in predicted_whole_data_freq:
                last2cycles = strain[col_index_stress][-257:]

                i=0
                for elem in last2cycles:
                    pred_extracted_ydata_freq.append(elem)
            
            for elem in pred_extracted_ydata_freq:
                pred_extracted_ydata.append(elem)


        #Read Experimental Data
        extracted_ydata = []
        max_val = []
        for freqvals in freqselect:

            if freqvals == freq[0] : freqindex = 0
            if freqvals == freq[1] : freqindex = 1
            if freqvals == freq[2] : freqindex = 2

            for i_st in range(0,len(strain_values_nonlin)): 
                last2cycles = all_data[freqindex*len(strain_values_nonlin)+i_st][1][-257:]
                for elem in last2cycles:
                    extracted_ydata.append(elem)

                maxelem = max(all_data[freqindex*len(strain_values_nonlin)+i_st][1][-257:]) 
                arraymaxelem = [maxelem]*257 #spesific
                        
                for elem in arraymaxelem:
                    max_val.append(elem)


        # print(f'pred:{len(pred_extracted_ydata)}')
        # print(f'exp:{len(extracted_ydata)}')
                                
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


    os.system('rm strain_sweep.txt')
    #Minimize the objective function for non-linear region
    result = minimize(fun=objective_maxnorm, x0=initial_guess_nonlin, method='Nelder-Mead') #, options={'xatol': 1e-5, 'fatol': 1e-4})  
    optimized_params_nonlin = result.x
    min_obj_nonlin = result.fun

    if result.success:
        print('Minimization is finished successfully for the non-linear region') 
    else: 
        print('Minimization is not successful')

    return optimized_params_nonlin, min_obj_nonlin, result

def steady_shearing(params, shear_rate_vals, eta_inf):      

    generate_input=1
    run_the_cases=1
    read_the_cases=1

    #Part0: Set the time array
    ntimesteps = 60000
    mintime = 0
    maxtime = 1.5
    ndatapoints = ntimesteps+1


    times = np.logspace(mintime,maxtime,ndatapoints)
    timesteps = []
    for i in range(0,len(times)-1):
        timesteps.append(times[i+1]-times[i])

    #print(max(timesteps))
    

    # Part 1: generate input 
    if generate_input == 1:
    
        os.system('rm -r inputs_steady')
        os.system('rm -r results_steady')
        os.system('mkdir -p inputs_steady')
        os.system('mkdir -p results_steady')
        os.system('mkdir results_steady/terminal_output')
        os.system('mkdir results_steady/temporal_output')

        modulus,t_rel, tau_y,nexp,alpha,gammadot_cr,I1c  = params
        #print(params)
               
        ii=0
        for gammadot in shear_rate_vals:

            f = open('inputs_steady/input{:04}.dat'.format(ii), 'w')       #write the input file
            f.write('&comppar\n\
    modulus = {}       \n\
    lambda = {}        \n\
    tau_y = {}         \n\
    nexp = {}          \n\
    alpha = {}         \n\
    gammadot_cr = {}   \n\
    eta_inf = {}       \n\
    I1c = {}           \n\
    gamma_dot = {}       \n\
    ntimesteps = {}       \n\
    timestepval = {}  \n\
    filename = "fitout_{:04}"   \n\
/\n\n\n\n'.format(modulus,t_rel, tau_y,nexp, alpha, gammadot_cr, eta_inf, I1c, gammadot,\
                        ntimesteps,timesteps,ii))
            f.close
            ii+=1
                
        #generate one extra set of input for the last one since very last iput file has an issue with reading
        
        f = open('inputs_steady/inputextra{:04}.dat'.format(ii), 'w')       #write the input file
        f.write('&comppar\n\
    modulus = {}       \n\
    lambda = {}        \n\
    tau_y = {}         \n\
    nexp = {}          \n\
    alpha = {}         \n\
    gammadot_cr = {}   \n\
    eta_inf = {}       \n\
    I1c = {}           \n\
    gamma_dot = {}       \n\
    ntimesteps = {}       \n\
    timestepval = {}  \n\
    filename = "fitout_{:04}"   \n\
/\n\n\n\n'.format(modulus,t_rel, tau_y,nexp, alpha, gammadot_cr, eta_inf, I1c, gammadot,\
                        ntimesteps,timesteps,ii))
        f.close


        input_file_list = sorted(os.listdir('inputs_steady'))

        for input in input_file_list:
            
            with open(f"inputs_steady/{input}", "r") as fp:
                text = fp.readlines()

                if len(text) >= 12:  # Check if the file has at least 12 lines
                    text[11] = text[11].replace('[', "").replace(']', "")

                    with open(f"inputs_steady/{input}", "w") as fw:
                        for line in text:
                            fw.write(line)
            
            if input == input_file_list[-2]:
                os.system(f'cp  inputs_steady/{input_file_list[-2]} inputs_steady/input{ii:04}.dat')


        os.system(f'rm  inputs_steady/inputextra*')


    # Part 2: run the cases   
    if run_the_cases == 1:

        #os.system('make')
       
        ii=0
        for gammadot in shear_rate_vals:

            os.system("./startup_fitting_steady < inputs_steady/input{:04}.dat  > results_steady/terminal_output/{:04}".format(ii,ii))
            os.system('mv fitout_* results_steady/temporal_output')

            ii+=1 

    # Part 3: read the cases
    if read_the_cases == 1:

        path_list=[]
        col_names=['step', 'time', 'strain', 'tauxx', 'tauxy', 'tauyy',\
                    'vonmis','gammadot', 'gammadot_eff', 'eta_eff', 'G_mod']
        

        res_list=sorted(os.listdir('results_steady/temporal_output'))
        sample=[]
            
        for res in res_list:
            sample_interval=[]

            path = 'results_steady/temporal_output/{}'.format(res)
            var1 = pd.read_fwf(filepath_or_buffer=path , header=None, names=col_names,\
                                widths=[12,12,12,12,12,12,12,12,12,12,12])
            
            time = var1['time'].tolist()
            strain = var1['strain'].tolist()
            xy = var1['tauxy'].tolist()

                
            element = [time,strain,xy]
            sample.append(element)
            
    return sample

def extract_predicted_G1G2(strains):

    col_names=['freq', 'G1', 'G2']
    sampleG1=[]   
    sampleG2=[]

    #Linear region
    for i in range(0,len(strains)):
        path = f'results/strain_sweep_{i:04}.txt'
        var1 = pd.read_fwf(filepath_or_buffer=path , header=None, names=col_names,\
                                widths=[12,12,12])
                
        G1 = var1['G1'].tolist()
        G2 = var1['G2'].tolist()
            
        sampleG1.append(G1[-1])
        sampleG2.append(G2[-1])
            
    return sampleG1, sampleG2

# working on the HB fitting of LAOS and stress growth data
def shear_rate_sweep(target_strain, data_exp_stgrw, i_strain, i_stress):
    
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

def HB(x, tau0, K, n):
    return tau0 + K*np.power(x,n)

def sum_of_squared_residuals_lin_modulus(params_actual, strain_lin, freqvals, freq, cycles, ntimesteps_cycle, data_exp_frswp):
    
    G1pred = []
    G2pred = []
    for strain in strain_lin:

        #Generate Model Data for given strains in linear region
        G1pred_elem, G2pred_elem = fitting_func_SAOS_freq_and_strain(params_actual, strain=strain,\
                        freq=freq,cycles=cycles,\
                        ntimesteps_cycle=ntimesteps_cycle)
        G1pred.append(G1pred_elem)
        G2pred.append(G2pred_elem)

    

    G1exp = data_exp_frswp[0][1] #there is only one elemnt in the set
    G2exp = data_exp_frswp[0][2]

    if len(freqvals) > 1: #freq sweep data
        
        #get the index for the frequency
        l=0
        for freqval in freqvals:
            if freqval == freq:
                index=l
            l+=1

        #difference between exp and model
        sumsqdiff = []
        for k in range(0,len(G1pred)):
            diffG1 = ( (G1exp[index]- G1pred[k]) / G1exp[index] )**2
            diffG2 = ( (G2exp[index]- G2pred[k]) / G2exp[index] )**2

    if len(freqvals) == 1 and len(strain_lin) > 1: #ampl sweep data
        #difference between exp and model
        sumsqdiff = []
        for k in range(0,len(G1pred)):
            diffG1 = ( (G1exp[k]- G1pred[k]) / G1exp[k] )**2
            diffG2 = ( (G2exp[k]- G2pred[k]) / G2exp[k] )**2
            

    sumsqdiff.append(diffG1)
    sumsqdiff.append(diffG2)

    meansumsqdiff = np.mean(sumsqdiff)
    
    return meansumsqdiff

def sum_of_squared_residuals_perfreq(params, freq_vals, strain_val_lin,\
                                     wholefreq, cycles,ntimesteps_cycle, col_index_stress,\
                                     all_data):
    
    ## Experimental ##
    #find the right index forthe freq sweep data
    extracted_ydata = []
    max_val = []
    l=0
    for freqval in wholefreq:
        if freq_vals[0] == freqval:
            index=l
        l+=1

    #extract the data
    last2cycles = all_data[0][0][1][index*257:(index+1)*257] #waveform data at given frequency
    gammadotlist = all_data[0][0][3][index*257:(index+1)*257] #gammadot for inputtin into simulation
    for elem in last2cycles:
        extracted_ydata.append(elem)
    maxelem = max(last2cycles) 
    arraymaxelem = [maxelem]*257 
    for elem in arraymaxelem:
        max_val.append(elem)

    ## Predicted ##

    #generate model data for given strains in linear region
    predicted_whole_data = fitting_func_SAOS_perfreq(params,strain_val_lin,\
                    freq_vals=freq_vals,cycles=cycles,\
                    ntimesteps_cycle=ntimesteps_cycle, gammadotlist=gammadotlist)

    #read model data 
    pred_extracted_ydata = []
    for m in range(0,len(freq_vals)):
        last2cycles = predicted_whole_data[m][col_index_stress][-257:]
        for elem in last2cycles:
            pred_extracted_ydata.append(elem)


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

def sum_of_squared_residuals_laos(params_actual, strain_values_nonlin, freq, cycles, eta_inf, ntimesteps_cycle,\
                                   optimized_params_lin, col_index_stress, all_data ):
        
    #doing everything per frequency
    pred_extracted_ydata = [] #all data will be here
    for freqvals in freq:
        modulus = optimized_params_lin[0][0]*freqvals**optimized_params_lin[0][1]
        t_rel = optimized_params_lin[1][0]*freqvals**optimized_params_lin[1][1]

        if freqvals == freq[0] : freqindex = 0
        if freqvals == freq[1] : freqindex = 1
        if freqvals == freq[2] : freqindex = 2

        #Generate Model Data for all strains of a freq
        predicted_whole_data_freq = fitting_func_LAOS(params_actual,strain_values=strain_values_nonlin,\
                        freq=freqvals,cycles=cycles,eta_inf=eta_inf,\
                        ntimesteps_cycle=ntimesteps_cycle, modulus=modulus,t_rel=t_rel,\
                        data_exp=all_data, freqindex=freqindex)
        
        #Read Model Data 
        pred_extracted_ydata_freq = []
        for strain in predicted_whole_data_freq:
            last2cycles = strain[col_index_stress][-257:]

            i=0
            for elem in last2cycles:
                pred_extracted_ydata_freq.append(elem)
        
        for elem in pred_extracted_ydata_freq:
            pred_extracted_ydata.append(elem)


    #Read Experimental Data
    extracted_ydata = []
    max_val = []
    for i_fq in range(0,len(freq)):
        for i_st in range(0,len(strain_values_nonlin)): 
            last2cycles = all_data[i_fq*7+i_st][1][-257:]
            for elem in last2cycles:
                extracted_ydata.append(elem)

            maxelem = max(all_data[i_fq*7+i_st][1][-257:]) 
            arraymaxelem = [maxelem]*257 #spesific
                    
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

def extract_steady_from_laos(data_exp_avg):

    ultimate_shrate = []
    ultimate_stress = []

    for k in range(0,len(data_exp_avg)):

        x = data_exp_avg[k][3][-257:] #shear rate
        y = data_exp_avg[k][1][-257:] #shear stress

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

def sum_of_squared_residuals_laos_perfreq(params_actual, strain_values_nonlin, freq, cycles, eta_inf, ntimesteps_cycle,\
                                   optimized_params_lin, col_index_stress, all_data, freqselect ):
        
    #doing everything per frequency
    pred_extracted_ydata = [] #all data will be here
    for freqvals in freqselect:
        modulus = optimized_params_lin[0][0]*freqvals**optimized_params_lin[0][1]
        t_rel = optimized_params_lin[1][0]*freqvals**optimized_params_lin[1][1]

        if freqvals == freq[0] : freqindex = 0
        if freqvals == freq[1] : freqindex = 1
        if freqvals == freq[2] : freqindex = 2

        #Generate Model Data for all strains of a freq
        predicted_whole_data_freq = fitting_func_LAOS(params_actual,strain_values=strain_values_nonlin,\
                        freq=freqvals,cycles=cycles,eta_inf=eta_inf,\
                        ntimesteps_cycle=ntimesteps_cycle, modulus=modulus,t_rel=t_rel,\
                        data_exp=all_data, freqindex=freqindex)
        
        #Read Model Data 
        pred_extracted_ydata_freq = []
        for strain in predicted_whole_data_freq:
            last2cycles = strain[col_index_stress][-257:]

            i=0
            for elem in last2cycles:
                pred_extracted_ydata_freq.append(elem)
        
        for elem in pred_extracted_ydata_freq:
            pred_extracted_ydata.append(elem)


    #Read Experimental Data
    extracted_ydata = []
    max_val = []
    for freqvals in freqselect:

        if freqvals == freq[0] : freqindex = 0
        if freqvals == freq[1] : freqindex = 1
        if freqvals == freq[2] : freqindex = 2

        for i_st in range(0,len(strain_values_nonlin)): 
            last2cycles = all_data[freqindex*len(strain_values_nonlin)+i_st][1][-257:]
            for elem in last2cycles:
                extracted_ydata.append(elem)

            maxelem = max(all_data[freqindex*len(strain_values_nonlin)+i_st][1][-257:]) 
            arraymaxelem = [maxelem]*257 #spesific
                    
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

def power_law(x,k,n): 
    return np.multiply(k,np.power(x,n))

def linear(x,a,b):
    return np.multiply(x,a)+b

def replace_nan_with_avg(arr):
    for i in range(1, len(arr)-1):
        if np.isnan(arr[i]):
            # Replace NaN with the average of the previous and next elements
            arr[i] = (arr[i-1] + arr[i+1]) / 2
    return arr