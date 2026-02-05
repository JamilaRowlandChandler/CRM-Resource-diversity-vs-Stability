# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 11:50:49 2025

@author: jamil
"""


import numpy as np
import pandas as pd
import seaborn as sns
import os
import sys
from scipy.optimize import curve_fit

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patheffects as patheffects

abspath = os.path.abspath(__file__)
file_directory_name = os.path.dirname(abspath)
os.chdir(file_directory_name)

sys.path.insert(0, file_directory_name.removesuffix("\\figures") + \
                "\\Modules")
from simulation_functions import le_pivot_r

# %%

################################ sigma_c ################################

df_simulation_c = pd.read_csv(file_directory_name.removesuffix("\\figures") + \
                            "\\Data\\simulation_dataframes\\M_vs_sigma_c.csv")
    
globally_solved_sces_c = pd.read_pickle(file_directory_name.removesuffix("\\figures") + \
                                      "\\Data\\self_consistency_equations\\M_vs_sigma_c.pkl") 

solved_boundary_c = pd.read_pickle(file_directory_name.removesuffix("\\figures") + \
                                 "\\Data\\self_consistency_equations\\stability_bound\\M_vs_sigma_c.pkl")

# %%

################################ sigma_y ################################

df_simulation_y = pd.read_csv(file_directory_name.removesuffix("\\figures") + \
                            "\\Data\\simulation_dataframes\\M_vs_sigma_y.csv")
    
globally_solved_sces_y = pd.read_pickle(file_directory_name.removesuffix("\\figures") + \
                                      "\\Data\\self_consistency_equations\\M_vs_sigma_y.pkl") 

solved_boundary_y = pd.read_pickle(file_directory_name.removesuffix("\\figures") + \
                                 "\\Data\\self_consistency_equations\\stability_bound\\M_vs_sigma_y.pkl")

# %%

############# Phase diagram + analytically-derived boundary #####################

def Stability_Plot(df_simulation_c, globally_solved_sces_c, solved_boundary_c,
                   df_simulation_y, globally_solved_sces_y, solved_boundary_y):
    
    def hyperbolic_fit(variable, boundary):
        
        fit_p, _ = curve_fit(hyperbolic, boundary['M'], boundary[variable],
                             bounds = [0, [1e6, 1e6]])
        
        return fit_p
    
    def hyperbolic(x, a, b):
        
        return a + b/x
    
    ###################
    
    def quadratic_fit(variable, boundary):
        
        smoother = np.poly1d(np.polyfit(boundary['M'], boundary[variable], 2))
        
        return smoother
        
    def sigma_vs_stability(ax, variable, pivot, sigmas, xlabel,
                           example_M = 225):
        
        stability = pivot.loc[:, example_M].to_frame()
        stability.reset_index(inplace = True)
        stability.rename(columns = {example_M : 'P(stability)'}, inplace = True)
        
        sns.lineplot(data = stability, x = variable, y = 'P(stability)',
                     ax = ax, linewidth = 3, color = 'black',
                     err_style = "bars", errorbar = ("pi", 100),
                     marker = "o", markersize = 9)
        
        se_95 = 1.96*np.sqrt((stability['P(stability)'] * \
                              (1 - stability['P(stability)']))/20)
            
        # add standard error bars on proportions (95 % confidence interval)
        # 20 = number of samples/communities 
        error_bars = ax.errorbar(x = stability[variable], 
                                 y =  stability['P(stability)'],
                                 yerr = se_95,
                                 fmt = 'none',
                                 ecolor = 'black',
                                 linewidth = 2) 

        ax.set_xticks(sigmas[::2], labels = sigmas[::2], fontsize = 10, 
                      rotation = 0)
        
        ax.yaxis.set_tick_params(labelsize = 10)
        
        ax.set_xlabel(xlabel, fontsize = 10, weight = 'bold')
        ax.set_ylabel('Probability(stability)', fontsize = 10, weight = 'bold')
        
        sns.despine(ax = ax)
        
    ############################################
    
    def stability_condition(ax, variable, sigmas, sces, boundary, boundary_method,
                            xlabel, example_M = 225):

        df_plot = sces.loc[sces['M'] == example_M, :]
        
        dfl = pd.melt(df_plot[[variable, 'rho', 'Species packing']], [variable])
        
        dfl.loc[dfl['variable'] == 'Species packing', 'value'] = \
            np.sqrt(dfl.loc[dfl['variable'] == 'Species packing', 'value'])
            
        good_solves = boundary.loc[boundary['loss'] <= -28, :]
            
        match boundary_method:
            
            case "hyperbolic":
                
                fit_p = hyperbolic_fit(variable, good_solves)
                fitted_boundary = hyperbolic(example_M, *fit_p)
            
            case "quadratic":
                
                smoother = quadratic_fit(variable, good_solves)
                fitted_boundary = smoother(example_M)
                
        if variable == "sigma_c":
    
            ax.add_patch(Rectangle((np.min(sigmas) - 0.1, np.min(dfl['value']) - 0.1),
                                   fitted_boundary - np.min(sigmas) + 0.1,
                                   np.max(dfl['value']) - np.min(dfl['value']) + 0.12,
                                   fill = True, color = '#8f8cc0ff', zorder = 0))
            
            ax.set_xlim([np.min(sigmas) - 0.1, np.max(sigmas) + 0.1])
            
        elif variable == "sigma_y":
            
            ax.add_patch(Rectangle((fitted_boundary, np.min(dfl['value']) - 0.1),
                                   np.max(sigmas) + 0.1 - fitted_boundary,
                                   np.max(dfl['value']) + 0.12 - np.min(dfl['value']),
                                   fill = True, color = '#8f8cc0ff', zorder = 0))
            
            ax.set_xlim([np.min(sigmas) - 0.01, np.max(sigmas) + 0.01])
        
        ax.vlines(fitted_boundary, np.min(dfl['value']) - 0.02,
                  np.max(dfl['value']) + 0.02,
                  color = 'black', linewidth = 2.5, zorder = 0)
        
        
        subfig1 = sns.lineplot(dfl, x = variable, y = 'value', hue = 'variable',
                               ax = ax, linewidth = 3,
                               palette = sns.color_palette(['black', 'black'], 2),
                               zorder = 10)
    
        subfig2 = sns.lineplot(dfl, x = variable, y = 'value', hue = 'variable',
                               ax = ax, linewidth = 2, marker = 'o', markersize = 8,
                               palette = sns.color_palette(['#00557aff', '#3dc27aff'], 2),
                               zorder = 10, markeredgewidth = 0.4, markeredgecolor = 'black')
        
        ax.set_ylim([np.min(dfl['value']) - 0.02, np.max(dfl['value']) + 0.02])
        
        ax.set_xlabel(xlabel, fontsize = 10, weight = 'bold')
        ax.set_ylabel('')
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.set_xticks(sigmas[::2], labels = [str(sigma) #+ '^2$' 
                                             for sigma in sigmas[::2]])
        
        ax.legend_.remove()
        

    #############################################################################################
    
    df_simulation_c = df_simulation_c.loc[df_simulation_c['sigma_c'] > 0.75, :]
    globally_solved_sces_c = globally_solved_sces_c.loc[globally_solved_sces_c['sigma_c'] > 0.75, :]
    solved_boundary_c = solved_boundary_c.loc[solved_boundary_c['sigma_c'] > 0.75, :]
    
    df_simulation_y = df_simulation_y.loc[df_simulation_y['sigma_y'] < 0.225, :]
    globally_solved_sces_y = globally_solved_sces_y.loc[globally_solved_sces_y['sigma_y'] < 0.225, :]
    solved_boundary_y = solved_boundary_y.loc[solved_boundary_y['sigma_y'] < 0.225, :]
    
    resource_pool_sizes = np.unique(df_simulation_c['M'])
    sigma_cs = np.unique(df_simulation_c['sigma_c'])
    sigma_ys = np.unique(df_simulation_y['sigma_y'])
    
    sns.set_style('ticks')
  
    mosaic = [["P_c", ".", "I_C_c"],
              ["P_y", ".", "I_C_y"]]
    
    fig, axs = plt.subplot_mosaic(mosaic, figsize = (6.5, 5.4),
                                  width_ratios = [5.2, 1.8, 5.4],
                                  height_ratios = [1, 1],
                                  gridspec_kw = {'hspace' : 0.5, 'wspace' : 0.1})
    
    ######################## Stability diagram ######################################
    
    # Simulation data
    
    stability_pivot_c = le_pivot_r(df_simulation_c, columns = 'M',
                                   index = 'sigma_c')[0]
    
    stability_pivot_y = le_pivot_r(df_simulation_y, columns = 'M',
                                     index = 'sigma_y')[0]
    
    sigma_vs_stability(axs["P_c"], 'sigma_c', stability_pivot_c, sigma_cs,
                       'std. deviation in total\nconsumption rate, ' + r'$\sigma_c$')
    
    sigma_vs_stability(axs["P_y"], 'sigma_y', stability_pivot_y, sigma_ys,
                       'std. deviation in yield\nconversion, ' + r'$\sigma_y$')
    
    axs["P_c"].set_title("Different sources of interaction heterogeneity" + \
                          "\ninduce opposing stability transitions",
                          fontsize = 11, weight = "bold", y = 1.1)
        
    #axs["P_c"].set_xlabel("")
        
    #################### Instability condition vs M #####################
    
    stability_condition(axs["I_C_c"], 'sigma_c', sigma_cs, globally_solved_sces_c,
                        solved_boundary_c, "hyperbolic",
                        'std. deviation in total\nconsumption rate, ' + r'$\sigma_c$')
    
    stability_condition(axs["I_C_y"], 'sigma_y', sigma_ys, globally_solved_sces_y,
                        solved_boundary_y, "quadratic",
                        'std. deviation in yield\nconversion, ' + r'$\sigma_y$')
    
    axs["I_C_c"].set_title("by having opposing " + \
                           "effects on\nreciprocity and the packing ratio",
                           fontsize = 11, weight = "bold", y = 1.1)
    
    ###############################################
    
    plt.show()

Stability_Plot(df_simulation_c, globally_solved_sces_c, solved_boundary_c,
               df_simulation_y, globally_solved_sces_y, solved_boundary_y)
