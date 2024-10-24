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

