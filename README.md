# R PROJECT README
# Stochastic Compartmental Model for SARS-COV-2 (OCTA)

This R project is conceived by members of SIMID group during the COVID-19 pandemic to estimate and explore the effect of social distancing and vaccine uptake in Belgium. OCTA stands for st**O**chastic **C**ompartment model for infec**T**ious dise**A**ses.

Current Version: 1.0.0

**Scientific output** related to this model:

* [Modeling the early phase of the Belgian COVID-19 epidemic (Abrams et al. Epidemics 2021)](https://doi.org/10.1016/j.epidem.2021.100449)
* [Technical Notes](https://www.simid.be/news/technical-note-covid19/)
* [RESTORE reports](https://covid-en-wetenschap.github.io/restore)
* The impact of quality-adjusted life years on evaluating COVID-19 mitigation strategies: lessons from age-specific vaccination roll-out and variants of concern in Belgium (2020-2022) (Willem et al. BMC Public Health, 2024)

**GitHub repository**

We have developed and continuously managed this R project using GitHub for robust version control, complemented by unit-testing to ensure the accuracy of our model throughout its development. Upon acceptance of our manuscript on vaccination roll-out and variants of concern, we made all the code and input data necessary to replicate the published scenarios available in a public archive on ZENODO.

**Project directory structure**


| Folder          | Content                                                                                     |
|-----------------|---------------------------------------------------------------------------------------------|
| data            | Houses all input data and configuration files. Temporary files (e.g., reported cases or hospital admissions from Sciensano) are stored in a "download" subfolder and are excluded from the git repository. |
| doc             | Contains a user manual and documentation on the incorporation of social contact data                       |
| R               | Holds the main R files for running the stochastic model, estimating parameters, and plotting results. |
| unit_testing    | Contains reference model output used in the main_workbench for comparing against the latest runs. |


**Main project files**

| File                                       | Content                                                                                     |
|--------------------------------------------|---------------------------------------------------------------------------------------------|
| projections_wave1.R                      | This script runs the stochastic model for the first COVID-19 wave in 2020 in Belgium. |
| projections_contacts.R                   | This script runs the stochastic model for Belgium, varying social contact patterns. |
| projections_uptake.R                     | This script runs the stochastic model for Belgium, varying vaccine uptake files. |
| plot_manuscript.R                        | This script (re)creates simulation output figures and/or combines different simulation output.   |
| calibration_mcmc.R                       | Use this script to initiate parameter estimations using MCMC (Markov Chain Monte Carlo).    |
| main_vaccination.R                       | This script loads all principal data and functions essential for the stochastic model.    |
| lib_model_core.R                         | Contains the core functions of the stochastic transmission model, including the main function and the log-likelihood function. |
| lib_xxx.R                                | Contains help functions for data manipulation, to visualize simulation output, case studies, or parameter estimation. |
| workbench_testing.R                      | This script runs different versions of the stochastic transmission model and verifies the latest model output against previously stored results. |


**Where to start?**

To get started with the OCTA model, we recommend beginning with the `R/projections_wave1.R` and `R/projections_contacts.R` file, which provide a basic introduction to running the transmission model. For comprehensive guidance, refer to the `doc/user_manual.pdf`.

To explore further, the `R/projections_uptake.R` and `R/plot_manuscript.R` files contain the simulation and plotting code used for our publication in BMC Public Health (2024). These scripts demonstrate more advanced uses of the model and how to generate outputs for academic dissemination.


**Contributors** (in alphabetical order):

* Steven Abrams 
* Philippe Beutels
* Christel Faes
* Nicolas Franco
* Niel Hens
* Lennert Verboven
* James Wambua
* Lander Willem (lander.willem@uantwerpen.be)
* The SIMID COVID-19 modelling team.
