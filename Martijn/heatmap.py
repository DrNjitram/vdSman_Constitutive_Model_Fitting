import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

data = np.loadtxt("heat_map.tsv", skiprows=1, usecols=[1, 2, 3, 4], delimiter="\t")
print(data)

fits = ["Tau" ,"n", "γ dot cr", "log G0"]
parameters = [
   # "Intercept",
    "T",
    "Moisture",
    "T*Moisture",
    "pH",
    "T^2",
    "pH^2"
]

estimate = data[1:6, :]
t_ratio = data[8:, :]

cmap = mpl.colors.ListedColormap(['#FF4F4F', '#EA4FFF', '#4F4FFF'])
fig, ax = plt.subplots()
t_ratio_filtered = t_ratio.copy()
t_ratio_filtered[np.abs(t_ratio)>=2] = 1
t_ratio_filtered[np.abs(t_ratio)>=3] = 2
t_ratio_filtered[np.abs(t_ratio)<2] = 0
im = ax.imshow(np.abs(t_ratio_filtered.T), cmap=cmap)

#Show all ticks and label them with the respective list entries
ax.set_xticks(range(len(parameters)), labels=parameters,
              rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(fits)), labels=fits)

#Loop over data dimensions and create text annotations.
for i in range(len(fits)):
    for j in range(len(parameters)):
        text = ax.text(j, i, t_ratio[j, i],
                       ha="center", va="center", color="k")

#ax.set_title("Harvest of local farmers (in tons/year)")
fig.tight_layout()
plt.show()