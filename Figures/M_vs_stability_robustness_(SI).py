# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 16:13:25 2026

@author: jamil
"""

import numpy as np
import pandas as pd
import os
import sys

from matplotlib import pyplot as plt

abspath = os.path.abspath(__file__)
file_directory_name = os.path.dirname(abspath)
os.chdir(file_directory_name)

sys.path.insert(0, file_directory_name.removesuffix("\\figures") + \
                "\\Modules")
from simulation_functions import le_pivot_r, generic_heatmaps

# %%

df_mu_c_M_si = pd.read_csv(file_directory_name.removesuffix("\\figures") + \
                           "\\Data\\simulation_dataframes\\M_vs_mu_c_consumer_inhibition.csv")
    
resource_pool_sizes = np.unique(df_mu_c_M_si['M'])
    
fig, axs = generic_heatmaps(df_mu_c_M_si,
                            'M', 'mu_c', 
                            'resource pool size, ' + r'$M$',
                            'average total consumption  rate, ' + r'$\mu_c$',
                            ['Max. lyapunov exponent'], 'Purples_r',
                            'Large resource pools still stabilise communities with\ndirect consumer self-inhibition',
                            (1, 1), (5.5, 5.5),
                            pivot_functions = {'Max. lyapunov exponent' : le_pivot_r},
                            specify_min_max={'Max. lyapunov exponent' : [0,1]})

axs.set_xticks(np.arange(0.5, len(resource_pool_sizes) + 0.5, 2),
                          labels = resource_pool_sizes[::2], fontsize = 14)

cbar = axs.collections[0].colorbar
cbar.set_label(label = 'Proportion of simulations with stable dynamics',
               size = '14')
cbar.ax.tick_params(labelsize = 12)

fig.suptitle("Large resource pools still stability communities with\n" + \
             "direct consumer self-inhibition",
             fontsize = 14, weight = "bold")

plt.show()

# %%

df_mu_c_M_fs = pd.read_csv(file_directory_name.removesuffix("\\figures") + \
                           "\\Data\\simulation_dataframes\\M_vs_mu_c_fixed_supply.csv")
    
fig, axs = generic_heatmaps(df_mu_c_M_fs,
                            'M', 'mu_c', 
                            'resource pool size, ' + r'$M$',
                            'average total consumption  rate, ' + r'$\mu_c$',
                            ['Max. lyapunov exponent'], 'Purples_r',
                            'Large resource pools still stabilise communities with\ndirect consumer self-inhibition',
                            (1, 1), (5.5, 5.5),
                            pivot_functions = {'Max. lyapunov exponent' : le_pivot_r},
                            specify_min_max={'Max. lyapunov exponent' : [0,1]})

axs.set_xticks(np.arange(0.5, len(resource_pool_sizes) + 0.5, 2),
                          labels = resource_pool_sizes[::2], fontsize = 14)

cbar = axs.collections[0].colorbar
cbar.set_label(label = 'Proportion of simulations with stable dynamics',
               size = '14')
cbar.ax.tick_params(labelsize = 12)

fig.suptitle("Large resource pools still stability communities with\n" + \
             "under constant (total) resource flux",
             fontsize = 14, weight = "bold")


plt.show()
    
