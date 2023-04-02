## Plots the pairwise RMSD of the full backbone

import numpy as np
import matplotlib
import seaborn as sb
import matplotlib.pyplot as plt

data = np.fromfile('rmsd.dat', np.float32)
rmsd = data.reshape(int(np.sqrt(len(data))), int(np.sqrt(len(data))))

cmap = matplotlib.cm.get_cmap("jet", 6)
plt.figure(figsize=(12,8))
sb.heatmap(rmsd*10, cmap = cmap, vmax=15)
plt.xlabel('Mobile')
plt.ylabel('Reference')
plt.title('Frame-wise full Backbone RMSD  ($\AA$)')
plt.savefig('bbrmsd_map.jpg')