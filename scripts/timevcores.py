import numpy as np
import os
import matplotlib.pyplot as plt
cores = [1, 2, 4, 8, 16, 32]
realtimediv = [3.391, 3.391/2, 3.391/4, 3.391/8, 3.391/16, 3.391/32]
realtime = [3.391, 2.812, 2.374, 1.950, 1.669, 1.462]
plt.errorbar(cores, realtimediv,  linestyle = 'None', marker = 'o', color = 'green')
plt.errorbar(cores, realtime,  linestyle = 'None', marker = 'o', color = 'blue')
plt.xlabel('Number of Cores')
# plt.legend(loc = 'upper right', fontsize = 12)
plt.savefig('timevcoreplot.pdf', bbox_inches = 'tight')
