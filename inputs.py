#running from scratch?
run_all = 'yes' #yes or no

#if not running from scratch
checkpoint_file_date  = '2024_10_21_18_05' #yyyy_mm_dd_hh_mm

explaianation_on_the_run = ''

#Waveform universal
number_of_data_per_cycle = 256

#SAOS data inputs
saos_amplitude = 2
saos_frequencies = [] 
saos_number_of_cycles= 2


#LAOS data inputs
laos_amplitudes = []
laos_frequencies = []


#SAOS fitting inputs
modulus_boundary = []
trel_boundary = []
saos_ncycles = 2
initial_guess = [0.5,0.5]


#LAOS fitting inputs
target_cycle_in_fitting = 1
bounds_sigy = [1e2,1e44] 
bounds_nexp = [1e-1, 8e-1] 
bounds_alph =  [1e-1, 8e-1] 
bounds_gdcr = [1e-2,1e3] 
bounds_I1cr = [1e-4,1e3] 

##do you also want to perform flow informed fitting?
flow_informed = 'yes' #yes or no

###Do you want to extract HB parameters from LAOS?
extract_HB_param_from_laos = 'yes ' #yes or no

###Do you want to provide external HB parameters?
HB_external = 'yes'
HB_external_params = ['sigy', 'n', 'dgcr']



