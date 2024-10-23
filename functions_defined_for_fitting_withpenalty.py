
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

def optimization_lin_perfreq_python(initial_guess_lin, bounds_lin_norm, strain_val_lin, 
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

def Ruud_Single_Cycle(gamma_0,omega,G,sigma_Y,gammadot_cr,nxp,eta_w,tau_M,I1c,alpha, cycle, numtimesteps):

    #global A_11,A_22,A_33,A_12

    A_11 = 1 
    A_12 = 0
    A_22 = 1
    
    T_E = numtimesteps
    dt = cycle  / (T_E * omega)

    sigma = np.zeros(T_E)
    gammadot = np.zeros(T_E)
    time = np.zeros(T_E)
    G_val = np.zeros(T_E)
    eta_eff_val = np.zeros(T_E)

    
    sum_G1 = 0
    sum_G2 = 0

    
    for t in range (T_E):
    
        time[t] = t*dt    

        # define shearrate
        gammadot[t] = omega * gamma_0 * np.cos(omega*2*math.pi*t*dt) 

        #define the shearrate dependent viscosity
        eta_eff = ( sigma_Y/ np.abs(gammadot[t]) ) * ( 1 + (np.abs(gammadot[t])/gammadot_cr)**nxp )
        
        #define the damping function to use in strain dependent elastic modulus
        I1 = A_11 + A_22 -2 
        h_gamma = 1.0 / ( 1 + abs(I1/I1c)**alpha)

        G_val[t] = G*h_gamma
        eta_eff_val[t] = eta_eff+eta_w

        #Shear rate dependent viscosity + dtrain dependent elstic modulus = relaxation time
        tau_eff = (eta_eff+eta_w)/(G*h_gamma)

        #Finalize relaxation time by adding the solvent relaxation time
        itau = 1.0/tau_M + 1.0/tau_eff;  # Cf. Miyazaki ExtraMaxwell mode
        tau = 1.0/itau   
        
        #Solve for the structure tensor
        A_11 += dt*(2*gammadot[t]*A_12 - (A_11 - 1)/tau )
        A_12 += dt*(  gammadot[t]*A_22 - A_12/tau )
        #structure tensor is updated
        
       # det = A_11*A_22 - A_12*A_12       
        #if (t%50==0):
            #print(det)

        I1 = (A_11 + A_22 -2)
        h_gamma = 1.0/(1+(abs(I1/I1c))**alpha)

        #calcuate sigma by eqn 7 from Ruud's paper    
        
        sigma[t] = G * h_gamma * A_12 # eta * gammadot
        
        sum_G1 = sum_G1 + sigma[t] * np.sin(omega*2*math.pi*t* dt)
        sum_G2 = sum_G2 + sigma[t] * np.cos(omega*2*math.pi*t* dt)
    
    t_e = time[T_E-1]
    
    #Calculate first harmoincs G1,G2 by stress and strain, eqn 8 in Ruud's article
    G1 = np.abs(np.pi * sum_G1 / t_e / gamma_0)
    G2 = np.abs(np.pi * sum_G2 / t_e / gamma_0)

    return G1, G2, sigma, time, gammadot, G_val, eta_eff_val

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
            sum_G1 = sum_G1 + sigma_elem * np.sin(freq*t*dt*2*3.14)*dt
            sum_G2 = sum_G2 + sigma_elem * np.cos(freq*t*dt*2*3.14)*dt
        
        
        G1.append(np.abs(2*np.pi * sum_G1 /  strain_lin*np.pi))
        G2.append(np.abs(2*np.pi * sum_G2 /  strain_lin*np.pi))

    return G1, G2, sigma, time

def optimization_lin_perfreq_ruud(initial_guess_lin, bounds_lin_norm, strain_val_lin, 
                current_freq, wholefreq, cycles,ntimesteps_cycle, col_index_stress, all_data ):   

        
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
        G1, G2, sigma, time = Ruud_Single_Cycle_Oldroydb(strain_val_lin,current_freq,\
                                                           params_actual, cycles, ntimesteps_cycle, gammadotlist)
        
        #read model data 
        pred_extracted_ydata = sigma[-257:]


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

def sum_of_squared_residuals_ruud(params_actual, current_freq, strain_val_lin,\
                                     wholefreq, cycles,ntimesteps_cycle, col_index_stress,\
                                     all_data):
    
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
        last2cycles = all_data[0][0][col_index_stress][index*257:(index+1)*257] #waveform data at given frequency
        gammadotlist = all_data[0][0][3][index*257:(index+1)*257] #gammadot for inputtin into simulation
        for elem in last2cycles:
            extracted_ydata.append(elem)
        maxelem = max(last2cycles) 
        arraymaxelem = [maxelem]*257 
        for elem in arraymaxelem:
            max_val.append(elem)

        ## Predicted ##

        #generate model data for given strains in linear region
        G1, G2, sigma, time = Ruud_Single_Cycle_Oldroydb(strain_val_lin,current_freq,\
                                                           params_actual, cycles, ntimesteps_cycle, gammadotlist)
        
        #read model data 
        pred_extracted_ydata = sigma[-257:]


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

def optimization_nonlin_together_ruud(initial_guess_nonlin, bounds_nonlin_norm,
                 strain_values_nonlin, freq,cycles,eta_inf,ntimesteps_cycle,\
                 col_index_stress, all_data, optimized_params_lin):

    
    def objective_maxnorm(params):

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
        
        #Extract predicted data
        #doing everything per frequency
        pred_extracted_ydata = [] #all data will be here
        freqindex=0
        for freqvals in freq:
            modulus = optimized_params_lin[0][0]*freqvals**optimized_params_lin[0][1]
            t_rel = optimized_params_lin[1][0]*freqvals**optimized_params_lin[1][1]
            
            strain_index=0
            for strainval in strain_values_nonlin:

                gammadot_data =all_data[strain_index+len(strain_values_nonlin)*freqindex][3][-257:]
                G1, G2, sigma, time = Ruud_Single_Cycle_Full(modulus, t_rel, freqvals, strainval,\
                                                              params_actual, cycles, eta_inf,\
                                                              ntimesteps_cycle, gammadot_data)
                
                for elem in sigma[-257:]:
                    pred_extracted_ydata.append(elem)
                
                strain_index+=1

            freqindex  += 1


        #Read Experimental Data
        extracted_ydata = []
        max_val = []
        for i_fq in range(0,len(freq)):
            for i_st in range(0,len(strain_values_nonlin)): 
                last2cycles = all_data[i_fq*len(strain_values_nonlin)+i_st][col_index_stress][-257:]
                for elem in last2cycles:
                    extracted_ydata.append(elem)

                maxelem = max(all_data[i_fq*len(strain_values_nonlin)+i_st][col_index_stress][-257:]) 
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


    #Minimize the objective function for non-linear region
    result = minimize(fun=objective_maxnorm, x0=initial_guess_nonlin, method='Nelder-Mead') #, options={'xatol': 1e-4, 'fatol': 1e-4})  
    optimized_params_nonlin = result.x
    min_obj_nonlin = result.fun

    if result.success:
        print('Minimization is finished successfully for the non-linear region') 
    else: 
        print('Minimization is not successful')

    return optimized_params_nonlin, min_obj_nonlin, result

def Ruud_Single_Cycle_Full( G, tau_M, freq, strain, params, ncycles, eta_w,\
                            ntimesteps_percycle, gammadot_data):

# arams, strain_values,\
# freq,cycles,eta_inf,ntimesteps_cycle, modulus, t_rel, data_exp, freqindex
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
        for t in range (0,ntimesteps_percycle+1):
     
            gammadot = gammadot_data[t]

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
            sum_G1 = sum_G1 + sigma_elem * np.sin(freq*t*dt*2*3.14)*dt
            sum_G2 = sum_G2 + sigma_elem * np.cos(freq*t*dt*2*3.14)*dt
        
        
        G1.append(np.abs(2*np.pi * sum_G1 /  strain*np.pi))
        G2.append(np.abs(2*np.pi * sum_G2 /  strain*np.pi))

    return G1, G2, sigma, time

def sum_of_squared_residuals_laos_ruud(params_actual, strain_values_nonlin, freq, cycles, eta_inf, ntimesteps_cycle,\
                                   optimized_params_lin, col_index_stress, all_data ):
        
         #Extract predicted data
        #doing everything per frequency
        pred_extracted_ydata = [] #all data will be here
        freqindex=0
        for freqvals in freq:
            modulus = optimized_params_lin[0][0]*freqvals**optimized_params_lin[0][1]
            t_rel = optimized_params_lin[1][0]*freqvals**optimized_params_lin[1][1]
            
            strain_index=0
            for strainval in strain_values_nonlin:

                gammadot_data =all_data[strain_index+len(strain_values_nonlin)*freqindex][3][-257:]
                G1, G2, sigma, time = Ruud_Single_Cycle_Full(modulus, t_rel, freqvals, strainval,\
                                                              params_actual, cycles, eta_inf,\
                                                              ntimesteps_cycle, gammadot_data)
                
                for elem in sigma[-257:]:
                    pred_extracted_ydata.append(elem)
                
                strain_index+=1

            freqindex  += 1


        #Read Experimental Data
        extracted_ydata = []
        max_val = []
        for i_fq in range(0,len(freq)):
            for i_st in range(0,len(strain_values_nonlin)): 
                last2cycles = all_data[i_fq*len(strain_values_nonlin)+i_st][col_index_stress][-257:]
                for elem in last2cycles:
                    extracted_ydata.append(elem)

                maxelem = max(all_data[i_fq*len(strain_values_nonlin)+i_st][col_index_stress][-257:]) 
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
    
        return diff_sqr_sum 