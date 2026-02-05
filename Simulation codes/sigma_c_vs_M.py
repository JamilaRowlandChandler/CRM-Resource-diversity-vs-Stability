# -*- coding: utf-8 -*-
"""
Created on Wed May  7 19:27:31 2025

@author: jamil
"""

import numpy as np
import pandas as pd
import sys
import os
from matplotlib import pyplot as plt

# %%

abspath = os.path.abspath(__file__)
file_directory_name = os.path.dirname(abspath)
os.chdir(file_directory_name)

sys.path.insert(0, file_directory_name.removesuffix("\\simulation codes") + \
                "\\Modules")
from simulation_functions import CRM_across_parameter_space

sys.path.insert(0,  file_directory_name.removesuffix("\\simulation codes") + \
                "\\Modules\\Cavity method functions")
import self_consistency_equation_functions as sce

# %%

def M_effect_sigma_c(M_range, sigma_C_range, mu_C, n, fixed_parameters):
    
    parameters = generate_parameters(M_range, sigma_C_range, mu_C, n,
                                     fixed_parameters)
    
    CRM_across_parameter_space(parameters,
                               file_directory_name.removesuffix("\\simulation codes") + \
                                   "\\Data\\simulation_data\\M_vs_sigma_c",
                               ['M', 'sigma_c'])
                    
# %%

def generate_parameters(M_range, sigma_C_range, mu_C, n,
                        fixed_parameters):
    
    M_sigma_C_combinations = np.unique(sce.parameter_combinations([M_range,
                                                                sigma_C_range],
                                                               n),
                                    axis = 1)
    
    variable_parameters = np.vstack([M_sigma_C_combinations[0, :]/fixed_parameters['gamma'],
                                     M_sigma_C_combinations[0, :],
                                     np.repeat(mu_C, M_sigma_C_combinations.shape[1])/M_sigma_C_combinations[0, :],
                                     M_sigma_C_combinations[1, :]/np.sqrt(M_sigma_C_combinations[0, :])])

    # array of all parameter combinations
    parameters = sce.variable_fixed_parameters(variable_parameters,
                                               fixed_parameters,
                                               ['S', 'M', 'mu_c', 'sigma_c'])
    
    for parms in parameters:
        
        parms['S'] = np.int32(parms['S'])
        parms['M'] = np.int32(parms['M'])
    
    return parameters
    
# %%

resource_pool_sizes = np.arange(50, 275, 25)

# %%

# mu_c = 160

M_effect_sigma_c(resource_pool_sizes, (0.5, 2.5), 160,
                 11, {'mu_y': 1, 'sigma_y' : 1.6/np.sqrt(150), 'b' : 1,
                      'd' : 1, 'gamma' : 1})
