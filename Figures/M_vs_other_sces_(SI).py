# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 16:17:16 2026

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
from simulation_functions import generic_heatmaps

# %%

def plot_sces(sces):
    
    """
    
    Plot self-consistency equations
    
    """
    
    quantities = ['phi_N', 'N_mean', 'q_N',
                  'phi_R', 'R_mean', 'q_R',
                  'v_N', 'chi_R']
    
    titles = ['species survival fraction, ' + r'$\phi_N$',
              'average species abundance, ' + r'$\langle N \rangle$', 
              'fluctuations in species abundances, ' + r'$\langle N^2 \rangle$',
              'resource survival fraction, ' + r'$\phi_R$',
              'average resource abundance, ' + r'$\langle R \rangle$', 
              'fluctuations in resource abundances, ' + r'$\langle R^2 \rangle$',
              'Avg. susceptibility of a species to\nperturbations by the cavity ' + \
                  'variables, ' + r'$v^{(N)}$',
              'Avg. susceptibility of a resource to\nperturbations by the' + \
                  ' cavity variables, ' + r'$\chi^{(R)}$']
    
    fig, axs = generic_heatmaps(sces,
                                'M', 'mu_c',
                                'resource pool size, $M$',
                                'avaerage total consumption coefficient, $\mu_c$',
                                quantities,
                                ['Blues', 'Greens', 'Oranges',
                                 'Blues', 'Greens', 'Oranges',
                                 'magma', 'inferno'],
                                titles,
                                (3, 3),
                                (8.5, 7),
                                mosaic = [['phi_N', 'N_mean', 'N_squared'],
                                          ['phi_R', 'R_mean', 'R_squared'],
                                          ['v_N', 'chi_R', '.']],
                                gridspec_kw = {'hspace' : 0.25,
                                               'wspace' : 0.1})
    
    plt.show()

# %%

# read in data
sces = pd.read_pickle(file_directory_name.removesuffix("\\figures") + \
                      "\\Data\\self_consistency_equations\\M_vs_mu_c.pkl")
 
# plot data
plot_sces(sces)