import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np



PCA = pd.read_csv("C://Users//agarm//Desktop//TFM 2.0//TFM-2.0//Python//PCA10_Covariables.csv", index_col = 0)

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

x = PCA['EV1']
y = PCA['EV2']
z = PCA['EV3']

num_labels = len(np.unique(PCA['pop']))
colors =["#000000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF",
            "#00FFFF", "#FF8000", "#80FF00", "#8000FF", "#FFFF80", "#FF80FF",
            "#80FFFF", "#008000", "#0080FF", "#800080", "#808080", "#C0C0C0",
            "#A0A0A0", "#606060", "#404040", "#202020", "#AA12AA", "#158752", 
            "#CCCC56", "#949513"]
colors_dict = {etiq: color for etiq, color in zip(np.unique(PCA['pop']), colors)}
 
colores=[]

for label, color in colors_dict.items():
    label_indices = PCA['pop'] == label
    ax.scatter(x[label_indices], y[label_indices], z[label_indices], c=color, label=label)


ax.set_xlabel('Componente Principal 1')
ax.set_ylabel('Componente Principal 2')
ax.set_zlabel('Componente Principal 3')
ax.legend()

plt.show()
