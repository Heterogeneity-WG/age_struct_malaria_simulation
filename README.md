# Project Team

This repository is comprised of codes related to the project "Modeling Immunity to Malaria with an Age-Structured PDE Framework" arising from the American Mathematical Society's Mathematical Research Community on infectious diseases (2020-2022). The team members on this project are:

Lauren Childs (Virginia Tech), Christina Edholm (Scripps College), Denis Patterson (Princeton University), Joan Ponce (Arizona State University), Olivia Propser (University of Tennessee, Knoxville), Zhuolin Qu (University of Texas at San Antonio) and Lihong Zhao (University of California Merced)

# Codebase functionality and user guide

Start by running _set_workspace.m_

**Parameters for the model (including the vacciation scheme/type are set in the following files:**

- _Malaria_parameters_baseline.m_ sets all demographic parameter (fertility/mortality rates), disease progression rates, immunity parameters and paramters for mosquito fertility/mortality (including seasonality),
  - the user can switch between sterilizing immunity and "blood statge-type vaccination" by adjusting the parameter P.z (0 for sterilizing, 1 for blood stage),
  - the user specifies the vaccination rate with variable P.v0.
- _Malaria_parameters_transform_vac.m_ sets the vaccination scheme by specifying the age range for which vaccination is applied.

**I want to:**

1. Run a single simulation of the model -> _run_age_structured_Malaria_vac.m_
    - this will generate one run of the model with the vaccination scheme and parameters set as described above.
2. Compare different vaccination schemes or vary other parameters systematicaly -> _run_compare.m_
3. Run simulations and save results (quantities of interest and solutions) to text files -> _run_save_data_and_figure.m_
4. Calculate the efficacy of vaccination -> _run_efficacy.m_
5. Run sensitivity analysis on the model -> _run_SA.m_
6. Calibrate sigmoids (immunity response functions) -> _run_parameter_search.m_ -> _Data_Fitting/run_parameter_fit.m_

**NB** Some results are saved to named directories so the codebase folder should be downloaded in full.

# Summary of key files

Malaria_parameters_baseline.m (script) -> sets the basic model parameters as global variables. This file calls Malaria_parameters_transform.m to find the stable age distribution via calls to balance_fertility.m to find a balanced birth rate function, if this option is selected, and then a call to find_stable_age.m

_age_structured_Malaria_vac.m_ (function) -> time-stepping algorithm for solving the ODE-PDE system

_run_age_structured_Malaria_vac.m_ (script) -> Calls _age_structured_Malaria.m_ for time-stepping algorithm to solve the ODE-PDE system. Calls Malaria_parameters_baseline.m to set parameters and find stable age distribution. Calls Malaria_IC and Immunity IC for initial/boundary conditions

_run_SA.m_ (script) -> runs sensitivity analysis of the model against each of a list of selected parameters. Calls _Malaria_parameters_baseline.m_ to set the parameters

_run_parameter_search.m_ (script) -> version of run_age_structured_Malaria.m to search the parameter space and assess impact of parameter changes in the sigmoids psi, phi and rho on the immunity distributions and endemic equilibrium













