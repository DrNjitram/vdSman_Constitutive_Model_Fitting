import numpy as np
from matplotlib import pyplot as plt

from functions_fitting import read_data_CCR

color_palette = plt.cm.hsv(np.linspace(0, 1, 16))
n_header= 18 #header rows just from the begining to the first data
cols_btw_data = 12 #number20 of rows between two measurements
header_list=['Strain %', "G' [kPa]", "G'' [kPa]", 'G* [kPa]', 'Tan Delta', 'Temperature [°C]', 'Pressure [bar]', "Viscosity n' [Pa*s]", 'Viscosity n" [Pa*s]', 'Viscosity n* [Pa*s]', 'Shear Rate [1/sec]', 'Shear Stress [Pa]']
what_to_collect = ['Strain %', "G' [kPa]", "G'' [kPa]"]#[ 'Strain', 'Storage Modulus','Loss Modulus'  ]
path_read=r"D:\Data\Rheology\CCR\202411 DoE pH T M Alph 8\Overview 1.xls"

sample_index = [0, 8]
sample_name = ['M55 T130 0M', 'M55 T130 0M']

#   Function
infos, data_exp_strswp = read_data_CCR(n_header, header_list, what_to_collect, path_read, cols_btw_data)
data_exp_strswp = [data_exp_strswp[i] for i in sample_index]
print(data_exp_strswp)

G1rep1 = data_exp_strswp[0][1]
G1rep2 = data_exp_strswp[1][1]

G2rep1 = data_exp_strswp[0][2]
G2rep2 = data_exp_strswp[1][2]

avg = [[],[]]
std_str_swp = [[],[]]
for i in range(0,len(G1rep1)):
    avg[0].append(np.average([G1rep1[i], G1rep2[i]]))
    avg[1].append(np.average([G2rep1[i], G2rep2[i]]))
    std_str_swp[0].append(np.std([G1rep1[i], G1rep2[i]]))
    std_str_swp[1].append(np.std([G2rep1[i], G2rep2[i]]))

data_exp_strswp_avg = [[]]
for i in range(0,3):
    if i == 0: data_exp_strswp_avg[0].append(data_exp_strswp[0][0])
    if i == 1: data_exp_strswp_avg[0].append(avg[0])
    if i == 2: data_exp_strswp_avg[0].append(avg[1])

# Create a figure and an array of axes with 1 row and 3 columns
fig, axs = plt.subplots(figsize=(8, 6))
axs.set_xscale('log')
axs.set_yscale('log')

axs.set_xlabel('Strain Amplitude (%)')
axs.set_ylabel("G', G'' (Pa)")

#axs.set_ylim(1e2, 1e5)

axs.errorbar(data_exp_strswp_avg[0][0], data_exp_strswp_avg[0][1], yerr=std_str_swp[0], marker='o', color=color_palette[0],\
              label='Data 1', markersize=5, linewidth =0, ecolor = 'black', elinewidth=3)
axs.errorbar(data_exp_strswp_avg[0][0], data_exp_strswp_avg[0][2], yerr=std_str_swp[1], marker='o', color=color_palette[0],\
              label='Data 1', markersize=5, linewidth =0, ecolor = 'black', elinewidth=3, alpha =0.5)

plt.show()



