# Project Team

This repository is comprised of codes related to the project "Modeling Immunity to Malaria with an Age-Structured PDE Framework" arising from the American Mathematical Society's Mathematical Research Community on infectious diseases (2020-2022). The team members on this project are:

Lauren Childs (Virginia Tech), Christina Edholm (Scripps College), Denis Patterson (Princeton University), Joan Ponce (Arizona State University), Olivia Prosper (University of Tennessee, Knoxville), Zhuolin Qu (University of Texas at San Antonio) and Lihong Zhao (University of California Merced)

# Codebase functionality and user guide

Start by running _set_workspace.m_

**Parameters for the model (including the vaccination scheme/type are set in the following files:**

- _Malaria_parameters_baseline.m_ sets all demographic parameters (fertility/mortality rates), disease progression rates, immunity parameters, and parameters for mosquito fertility/mortality (including seasonality),
  - the user can switch between sterilizing immunity and "blood stage-type vaccination" by adjusting the parameter P.z (0 for sterilizing, 1 for blood stage),
  - the user specifies the default vaccination rate with variable P.v0.
  - the user can turn off the seasonality by uncommenting the line "P.ss_c = 1; P.ss_S0 = 1;".
- _Malaria_parameters_transform.m_ sets the functional form for the time-dependent seasonality, age-dependent demographics, and immunity functions.
- _Malaria_parameters_transform_vac.m_ sets the vaccination scheme by specifying the age range for which vaccination is applied.

**I want to:**

1. Run a single simulation of the model -> _run_age_structured_Malaria_vac.m_
    - this will generate one run of the model with the vaccination scheme and parameters set as described above.
2. Compare different vaccination schemes or vary other parameters systematically -> _run_compare.m_
3. Run simulations and save results (quantities of interest and solutions) to text files -> _run_save_data_and_figure.m_
4. Calculate the efficacy of vaccination -> _run_efficacy.m_
5. Run sensitivity analysis on the model 
  -> local and extended SA: _run_SA.m_  
  -> global SA: _run_SA_PRCC.m_ _run_SA_eFAST.m_
6. Calibrate sigmoids (immunity response functions) -> _run_parameter_search.m_ -> _Data_Fitting/run_parameter_fit.m_

**NB** Some results are saved to named directories so the codebase folder should be downloaded in full.

# Summary of key files and dependencies

_run_age_structured_Malaria_vac.m_ (script) 
 - Calls _Malaria_parameters_baseline.m_ _Malaria_parameters_transform.m_ and _Malaria_parameters_transform_vac.m_ to set parameters.
 - Calls _age_structured_Malaria_IC_vac.m_ and _steady_state_vac.m_ for initial condition
 - Calls _age_structured_Malaria_vac.m_ for time evolution of the ODE-PDE system.

_Malaria_parameters_baseline.m_ (script) -> sets the basic model parameters as global variables. This file is often run together with
 - _Malaria_parameters_transform.m_  
 - _Malaria_parameters_transform_vac.m_ 

_age_structured_Malaria_vac.m_ (function) -> time-stepping algorithm for solving the ODE-PDE system
 - Calls _mosquito_ODE.m_ for mosquito ODE subsystem
 - Calls _biting_rate.m_ for compromised biting rates between humans and mosquitoes; 
 - Calls _FOI_H.m_ and _FOI_M.m_ for force of infections
 - Calls _sigmoid_prob.m_ for immunity linking functions

_age_structured_Malaria_eff.m_ (function) -> time-stepping algorithm for solving the three-patch model for vaccine efficacy calculations

_run_SA.m_ (script) -> runs sensitivity analysis of the model against each of a list of selected parameters. Calls _Malaria_parameters_baseline.m_ to set the parameters

_run_SA_eFAST.m_ (script) -> runs the global eFAST sensitivity analysis. Calls the functions under the eFAST/ directory.

_run_parameter_search.m_ (script) -> version of run_age_structured_Malaria.m to search the parameter space and assess the impact of parameter changes in the sigmoids psi, phi, and rho on the immunity distributions and endemic equilibrium













