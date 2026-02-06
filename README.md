# Code repository for Rowland-Chandler et al., "Resource diversity begets stability in complex ecosystems".

To access the corresponding data files, please go to...
The code is written assuming the user has the data saved in a folder called "Data" in the same directory as these codes. 

## Repository contents:

### Example notebooks: Jupyter notebooks using the user how to do the following. (These codes are not used in the paper itself.)
- `Running simulations.ipynb`: Randomly generate a consumer-resource model community, run a simulation, and determine community stability.
- `Solving self consistency equations.ipynb`: Solve the self-consistency equations for a given parameter set (parameters determining the distribution of model coefficients), and solve for the stability boundary.

### Figures: Code used for generating figures in the Main text and Supplementary information.
- `eLV_vs_CRM.py`: Main text, figure 4; Supplementary information, figures S1, S2, and S5
- `heterogeneity_vs_stability.py`: Main text, figure 5
- `heterogeneity_vs_stabilit_eLV_(SI).py`: Supplementary information, figures S3, S4, and S6.
- `M_vs_stability.py`: Main text, figure 3
- `M_vs_other_sces.py`: Suppplementary information, figure S9
- `M_vs_stability_robustness_(SI)`: Supplementary information, figures S10 and S11
- `other_stability_transitions.py`: Supplementary information, figure S8

### Modules:
- #### Cavity method functions: Subdirectory of function files for
	- generating the self-consistency equations (`self_limiting_gc_v_finite_equations.py`)
 	- and numerically solving them (`self_consistency_equations_functions.py`).
- #### Consumer-resource models:
	- `models.py`: Contains the consumer-resource model classes. These classes inherit methods from the other files in the subdirectory excluding `effective_LV_models.py`.
	- `parameters.py`: Interface for randomly generating model coefficients from their respective distributions
	- `initial_abundances.py`: Interface for randomly generating initial species and resource abundances,
	- `differential_equations.py`: Interace for simulating dynamics. (Although the method describing each model's ODE is actually in `models.py` - sorry for the confusion)
	- `community_level_properties.py`: Interface for calculating emergent community properties like stability.
 	- `effective_LV_models.py`: Contains the effective Lotka-Volterra model class. eLVs are generated from already-initialised consumer-resource model objects. 
- #### `simulation_functions.py`: Collection of functions for simulating many communities with different parameter distributions, generating dataframes of community properties from simulations, and generating quick plots.

### SCEs codes: Codes for running the solver routine for the self-consistency equations (in "Cavity method functions")
- `all_solves.py` contains codes for solving the self-consistency equations for all varying parameters (in Main text fig. 3 and 5., and supplementary information fig. S8). 

### Simulation codes: Codes for running simulations (calling methods from "Consumer-resource models")
- `create_simulation_dfs.py`: Code for extracting parameter distributions and emergent properties from "communities" of the consumer-resource model classes
- `mu_c_vs_M.py`: Creates consumer-resource models and runs simulations for different resource pool sizes ($M$) and average total consumption coefficients ($\mu_c$).
- `mu_c_vs_M_consumer_si.py`: Same as `mu_c_vs_M.py`, except it generates models where consumers directly inhibit each other (as well as compete for resources).
- `mu_c_vs_M_egLV.py`: Creates effective Lotka-Volterra models from the consumer-resource models generated in `mu_c_vs_M.py`, then runs simulations.
- `mu_c_vs_M_scaled_supply.py`: Same as `mu_c_vs_M.py`, except it generates models where the intrinsic resource growth coefficient ($b_\alpha$) is also scaled with the resource pool size ($M$).
- `sigma_c_vs_M.py`: Creates consumer-resource models and runs simulations for different resource pool sizes ($M$) and standard deviations in consumption coefficients ($\sigma_c$).
- `sigma_vs_M_egLV`: Creates effective Lotka-Volterra models from the consumer-resource models generated in `sigma_c_vs_M.py` and `sigma_y_vs_M.py`, then runs simulations.
- `sigma_y_vs_M.py`: Creates consumer-resource models and runs simulations for different resource pool sizes ($M$) and standard deviations in the yield conversion factor ($\sigma_y$).
