# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 23:32:42 2025

@author: jamil
"""

# -*- coding: utf-8 -*-
"""
Created on Fri May  9 15:32:30 2025

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

def M_effect_sigma_y(M_range, sigma_y_range, mu_C, sigma_C, n, fixed_parameters):
    
    parameters = generate_parameters(M_range, sigma_y_range, mu_C, sigma_C,
                                     n, fixed_parameters)
    
    CRM_across_parameter_space(parameters,
                               file_directory_name.removesuffix("\\simulation codes") + \
                                   "\\Data\\simulation_data\\M_vs_sigma_y"
                               ['M', 'sigma_y'])
                
# %%

def generate_parameters(M_range, sigma_y_range, mu_C, sigma_C, n,
                        fixed_parameters):
    
    M_sigma_y_combinations = np.unique(sce.parameter_combinations([M_range,
                                                                   sigma_y_range],
                                                                  n),
                                    axis = 1)
    
    variable_parameters = np.vstack([M_sigma_y_combinations[0, :]/fixed_parameters['gamma'],
                                     M_sigma_y_combinations,
                                     np.repeat(mu_C, M_sigma_y_combinations.shape[1])/M_sigma_y_combinations[0, :],
                                     np.repeat(sigma_C, M_sigma_y_combinations.shape[1])/np.sqrt(M_sigma_y_combinations[0, :])])

    # array of all parameter combinations
    parameters = sce.variable_fixed_parameters(variable_parameters,
                                               fixed_parameters,
                                               ['S', 'M', 'sigma_y', 'mu_c', 'sigma_c'])
    
    for parms in parameters:
        
        parms['S'] = np.int32(parms['S'])
        parms['M'] = np.int32(parms['M'])
    
    return parameters
    
# %%

resource_pool_sizes = np.arange(50, 275, 25)

# %%

M_effect_sigma_y(resource_pool_sizes, (0.05, 0.25), 160, 1.6, 9,
                 {'mu_y' : 1, 'b' : 1, 'd' : 1, 'gamma' : 1})
