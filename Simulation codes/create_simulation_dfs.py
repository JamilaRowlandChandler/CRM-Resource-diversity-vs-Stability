# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 17:59:43 2025

@author: jamil

This script is not for users. I used it to generate dataframes from my simulation
data stored in my local files.
I haven't uploaded all the simulations to GitHub as they use too much storage.

But, you can simulate some data yourself using the scripts in "Simulation codes",
then save them like so

"""

import numpy as np
import pandas as pd
import sys
import os

# %%

abspath = os.path.abspath(__file__)
file_directory_name = os.path.dirname(abspath)
os.chdir(file_directory_name)

sys.path.insert(0, file_directory_name.removesuffix("\\simulation codes") + \
                "\\Modules")
from simulation_functions import generate_simulation_df

# %%

########################## M vs mu_c ##########################

df_mu_c_M = generate_simulation_df("C:/Users/jamil/Documents/PhD/Data files and figures/Ecological-Dynamics-and-Community-Selection/Ecological Dynamics/Data/" \
                                   + 'resource_diversity_stability/simulations/M_vs_mu_c') 
    # or replace with file_directory_name.removesuffix("\\simulation codes") + "\\Data\\simulation_data\\M_vs_mu_c"

df_mu_c_M.to_csv(file_directory_name.removesuffix("\\simulation codes") + \
                "\\Data\\simulation_dataframes\\M_vs_mu_c.csv")
    
del df_mu_c_M

# %%

########################## M vs mu_c (with direct consumer self-inhibition) ##########################

df_mu_c_M = generate_simulation_df("C:/Users/jamil/Documents/PhD/Data files and figures/Ecological-Dynamics-and-Community-Selection/Ecological Dynamics/Data/" \
                                   + 'resource_diversity_stability/simulations/M_vs_mu_c_consumer_inhibition') 
    # or replace with file_directory_name.removesuffix("\\simulation codes") + "\\Data\\simulation_data\\M_vs_mu_c_consumer_inhibition"

df_mu_c_M.to_csv(file_directory_name.removesuffix("\\simulation codes") + \
                "\\Data\\simulation_dataframes\\M_vs_mu_c_consumer_inhibition.csv")
    
del df_mu_c_M

# %%

########################## M vs mu_c (with fixed total resource supply) ##########################

df_mu_c_M = generate_simulation_df("C:/Users/jamil/Documents/PhD/Data files and figures/Ecological-Dynamics-and-Community-Selection/Ecological Dynamics/Data/" \
                                   + 'resource_diversity_stability/simulations/M_vs_mu_c_fixed_supply') 
    # or replace with file_directory_name.removesuffix("\\simulation codes") + "\\Data\\simulation_data\\M_vs_mu_c_fixed_supply"

df_mu_c_M.to_csv(file_directory_name.removesuffix("\\simulation codes") + \
                "\\Data\\simulation_dataframes\\M_vs_mu_c_fixed_supply.csv")
    
del df_mu_c_M

# %%

########################## M vs sigma_c ##########################

df_sigma_c_M = generate_simulation_df("C:/Users/jamil/Documents/PhD/Data files and figures/Ecological-Dynamics-and-Community-Selection/Ecological Dynamics/Data/" \
                                   + 'resource_diversity_stability/simulations/M_vs_sigma_c')
    
df_sigma_c_M.to_csv(file_directory_name.removesuffix("\\simulation codes") + \
                "\\Data\\simulation_dataframes\\M_vs_sigma_c.csv")
    
del df_sigma_c_M

# %%

########################## M vs sigma_y ##########################

df_sigma_y_M = generate_simulation_df("C:/Users/jamil/Documents/PhD/Data files and figures/Ecological-Dynamics-and-Community-Selection/Ecological Dynamics/Data/" \
                                   + 'resource_diversity_stability/simulations/M_vs_sigma_y')

df_sigma_y_M.to_csv(file_directory_name.removesuffix("\\simulation codes") + \
                "\\Data\\simulation_dataframes\\M_vs_sigma_y.csv")
    
del df_sigma_y_M