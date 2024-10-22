color_palette = plt.cm.tab20(np.linspace(0, 1, 20))
markersize=10

#Running the HB fitting for HB informed fitting
#Option 1: Read a previous run case
run_HB = 'no'
readHB = 'HBfit_2024_10_22_13_50'

#Option2: Run the HB fitting from scratch
run_HB = 'yes'
extract_HB_param_from_laos = 'yes' #use LAOS data?
data_index_for_laos_built_flow_curves = [0,1,2,7,8,9,14,15,16]
skip_first = 3 #skip the first * of laos extracted flow data

extract_HB_param_from_stress_growth = 'yes' #use stress growth data?
current_time = datetime.datetime.now()
formatted_time = current_time.strftime('%Y_%m_%d_%H_%M')

extract_HB_param_from_external = 'yes' #provide an external result?
HB_external_params = [1000, 10, 0.2] #sigma_y, dot_gamma_cr, n
error_HBext = [0,0,0] #sigma_y, dot_gamma_cr, n




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



#Stress Growth Data if exsits
str_grw_shear_rates = [0.063, 0.1, 0.17, 0.28, 0.45, 0.73, 1.2,\
                1.96, 3.21, 5.25, 8.58, 14,23,37,61]
target_strain =120 # to extract the stress at



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






