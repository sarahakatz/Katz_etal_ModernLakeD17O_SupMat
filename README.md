# Katz_etal_ModernLakeD17O_SI
# Supplemental code for “Using lake carbonate D'17O and D47 for reconstructing precipitation d18O in humid systems” by Sarah A. Katz, Naomi E. Levin, Donald T. Rodbell, David P. Gillikin, Phoebe G. Aron, and Benjamin H. Passey


### Title: Using lake carbonate D'17O and D47 for reconstructing precipitation d18O in humid systems
### Supplementary script: lake.bal.supplement.R
### Code Author: Sarah Katz
### Contact: Sarah Katz (skatzees@umich.edu) or Naomi Levin (nelevin@umich.edu) with questions or comments

#########################################################
##########  0) Description of script  ###################
#########################################################

###     This script models the isotopic composition (d18O, d17O, d2H, D'17O, d-excess) of evaporated lake waters
###     based on steady state mass balance equations. Here, we use both the traditional mass balance equations of 
###     Criss (1999) for flow-through lakes experiencing evaporation (Eq. 1), and equations that account for the
###     influence of evaporated waters on the isotopic composition of atmospheric moisture (Eq. 2; e.g., Benson 
###     and White, 1994; Passey and Ji, 2019).
###
###     This script models the isotopic composition of lake waters using a Monte Carlo approach, in which the user 
###     prescribes a numerical range for each model parameter. For each parameter, values are randomly generated and used 
###     to solve the water budget equations. This work is based upon the Monte Carlo approach used by Passey
###     and Ji (2019). Our work expands upon that work by including code to calculate lake water isotopes using
###     either of the balance equations and includes calculations for d2H and d-excess. Though  modeled d-excess is not
###     discussed in this paper, our code may be useful for studies which combine lake water d-excess and D'17O data.
###
###     This script also includes the code to generate Figures 6 and 7.


###     Contents:

###     0) Description of script
###     1) Required packages & datasets from this study and Passey & Ji (2019)
###     2) Lake budget steady state equations
###     3) Constants
###     4) User defined model parameters and random selection of parameter values
###     5) Steady state lake water calculations - Eq. 1
###     6) Steady state lake water calculations - Eq. 2
###     7) Plots of modeled lake waters and lambda_lake (used to create Fig. 6)
###     8) Calculate modeled lambda_lake - D'17Ow relationship
###     9) 3D plots (Used to create Fig. 7)
