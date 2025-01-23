import numpy as np
from matplotlib import pyplot as plt

from data_reading_inputs import ndata_cycle

path_read=r"W:\Data\Rheology\CCR\202411 DoE pH T M Alph 8\Overview 1.xls"
ndata_cycle = 35
color_palette = plt.cm.hsv(np.linspace(0, 1, 16))

font = {'family' : 'arial',
        'weight' : 'bold',
        'size'   : 24}

plt.rcParams.update({'font.size': 24})

