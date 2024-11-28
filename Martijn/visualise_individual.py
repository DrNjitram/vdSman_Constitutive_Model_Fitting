import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

def cyclical_mean(data, signal_length):
    filtered_data = data[:len(data)//signal_length * signal_length]
    periodic = data.reshape(len(filtered_data)//signal_length, signal_length)
    return periodic.mean(axis=0)

path = r"D:\Data\Rheology\CCR\202409 DOE pH T M Arcon F\Block4"
filepaths = []
for path, subdirs, files in os.walk(path):
    for file in files:
        filepaths.append(os.path.join(path, file))

fig, ax = plt.subplots()

colors = ["b", "r", "k"]
for i, extension in enumerate([".asc"]):
    for filepath in filepaths[::-1]:
        if filepath.endswith(extension):
            if "8_0" in filepath or "0_5" in filepath:
                continue
            data = np.loadtxt(filepath, delimiter=";")
            x_data = cyclical_mean(data[:, 0], 90)
            plt.plot(x_data/np.max(x_data), cyclical_mean(data[:, 1], 90), label=os.path.basename(filepath)[:-4].replace("_", "."))

#ax.legend(handles=[Patch(facecolor='blue',label='T=130'),Patch(facecolor='red', label='T=170')])
plt.xlabel("Shear Strain [-]")
plt.ylabel("Normalised Shear Stress, σ [Pa]")
plt.title("130°C")
plt.legend()
plt.show()