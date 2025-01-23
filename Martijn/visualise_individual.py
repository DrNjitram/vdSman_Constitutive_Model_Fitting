import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

def cyclical_mean(data, signal_length):
    filtered_data = data[:len(data)//signal_length * signal_length]
    periodic = data.reshape(len(filtered_data)//signal_length, signal_length)
    return periodic.mean(axis=0)



def plot_lissajouis(path, extension=(".asc"), title=""):
    filepaths = []
    for path, subdirs, files in os.walk(path):
        for file in files:
            filepaths.append(os.path.join(path, file))

    fig, ax = plt.subplots()

    colors = ["b", "r", "k"]
    for i, extension in enumerate(extension):
        for filepath in filepaths[::-1]:
            if filepath.endswith(extension):

                if "8_0" in filepath or "0_5" in filepath:
                    continue
                data = np.loadtxt(filepath, delimiter=";")
                x_data = cyclical_mean(data[:, 0], 90)
                plt.plot(x_data/np.max(x_data), cyclical_mean(data[:, 1], 90), label=os.path.basename(filepath).split(".")[0].replace("_", "."))

    #ax.legend(handles=[Patch(facecolor='blue',label='T=130'),Patch(facecolor='red', label='T=170')])
    plt.xlabel("Shear Strain [-]")
    plt.ylabel("Normalised Shear Stress, σ [Pa]")
    plt.title(title if title else ", ".join(extension))
    plt.legend(title="Max Amplitude")
    plt.show()

path = r"D:\Data\Rheology\CCR\202409 DOE pH T M Arcon F\Block4"
plot_lissajouis(path, title="Arcon F 130°C 55%")
path = r"D:\Data\Rheology\CCR\202411 DoE pH T M Alph 8\Block 1"
plot_lissajouis(path, ["3"], title="Alpha 8 100°C 55%")