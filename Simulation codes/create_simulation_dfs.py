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

def generate_and_save_simulation_df(subdirectory):
    
    df = generate_simulation_df(file_directory_name.removesuffix("\\simulation codes") + \
                                "\\Data\\simulation_data\\" + subdirectory)
        
    df.to_csv(file_directory_name.removesuffix("\\simulation codes") + \
              "\\Data\\simulation_dataframes\\" + subdirectory + ".csv")

# %%

########################## M vs mu_c ##########################

generate_and_save_simulation_df("M_vs_mu_c")

# %%

########################## M vs mu_c (with direct consumer self-inhibition) ##########################

generate_and_save_simulation_df("M_vs_mu_c_consumer_inhibition")

# %%

########################## M vs mu_c (with fixed total resource supply) ##########################

generate_and_save_simulation_df("M_vs_mu_c_fixed_supply")

# %%

########################## M vs sigma_c ##########################

generate_and_save_simulation_df("M_vs_sigma_c")

# %%

########################## M vs sigma_y ##########################

generate_and_save_simulation_df("M_vs_sigma_y")
