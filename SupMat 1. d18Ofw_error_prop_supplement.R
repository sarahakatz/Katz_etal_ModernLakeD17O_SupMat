### Title: Detecting hydrologic distinctions among Andean lakes using clumped and triple oxygen isotopes
### Supplementary script: Error propagation for d18O values of carbonate formation waters (d18Orfw) from carbonate d18Oc and D47 values
### Author: Sarah Katz
### Contact: Sarah Katz (skatzees@umich.edu) or Naomi Levin (nelevin@umich.edu) with questions or comments

################################
##  0) Description of script  ##
################################

###     This script uses Monte Carlo resampling to propagate uncertainty in reconstructed carbonate formation water
###     d18O calculations (d18Orfw) from carbonate d18Oc and D47 values. 

###     Contents: 
###     0) Description of script
###     1) Description of file format and data import
###     2) Equations for reference
###     3) For loop
###       3.1) For loop set-up
###       3.2) Open for loop
###       3.3) Gaussian re-sampling of d18Oc and D47 values 
###       3.4) Gaussian re-sampling of slope and intercept values for D47-temp transfer function
###       3.5) Calculate new T_D47 values from re-sampled D47 values, slopes, and intercepts
###       3.6) Calculates new vector of alpha values (for carb-water fractionation) from re-sampled T_D47 values
###       3.7) Calculate new d18Ofw values using re-sampled alpha and d18Oc values
###       3.8) Calculate re-sampled d18Ofw and T_D47 average and sd values
###     4) Compile means and standard deviations; create new data frame; export .csv file

####################
## 1) Import Data ##
####################

### File format and description of columns. File type should be .csv
###     Column A      Sample_ID     <- sample name
###     Column B      N             <- number of replicate analyses
###     Column C      d18Oc_pdb     <- carbonate d18O value, reported on the VPBD scale in units of per mil
###     Column D      d18Oc_stdev   <- 1 sigma standard deviation on 'd18Oc_pdb' among the replicate analyses
###     Column E      D47           <- carbonate D47 value, reported in units of per mil
###     Column F      D47_stdev     <- 1 sigma standard deviation on 'D47' among the replicate analyses


### Define path.to.data on your local drive. User should update based on their own file organizational structure
    path.to.data <- "~/Dropbox (University of Michigan)/Katz - Junin modern paper/R_scripts/"
### Define data file. Update file name, if different than the suggested file name.
    df <- list.files(path=path.to.data, pattern="d18Ofw_error_prop_R") 
### Import data
    df <- do.call("rbind", lapply(df, function(x) 
      read.csv(paste(path.to.data, x, sep=''), stringsAsFactors = FALSE)))

################################
## 2) Equations for reference ##
################################

### rnorm()                                                    ## Produce a normal distribution of values around a mean (inputs: N, mean, sd)

### d18Oc_smow = (1.03091*d18Oc_pdb) + 30.91                   ##convert d18O from pdb to smow reference frame
### D47 = 0.0422 (+- 0.0019) *10^6/T^2 + 0.1262 (+- 0.0207)    ## D47-temp transfer function (Bonifacie et al., 2017). Note 95% confidence intervals are given. 1 sigma SD on slope and intercept are +- 0.0003 and 0.0035 per mil, respectively 
### 1000ln(a) = 18.03*(10^3/T) - 32.42                         ## temperature-dependent fractionation factor between water and calcite (Kim and O'Neil, 1997)

    
## Option for manual data input. Comment-out if using for loop
## lab = c("sample_name", d18Oc_pdb", "d18Oc_err", "D47", "D47_err")
##dat = c(test, -13.351,	0.104,	0.6437,	0.0125) 


####################
## 3) For loop #####
####################   

### 3.1) For loop set-up
    
z = 10000                                                     ##desired number of re-samplings

mean_d18Orfw_rs <- vector()                                    ## Create empty vectors to hold products from for loop
sd_d18Orfw_rs <- vector()
mean_T_D47_rs <- vector()
sd_T_D47_rs <- vector()

### 3.2) Open for loop


for (i in 1:nrow(df)){


## 3.3) Gaussian re-sampling of d18Oc and D47 values

d18Oc_rs = rnorm(z, df[i,3], df[i,4])                         ## d18Oc_rs gives re-sampled d18Oc values in pdb.
d18Oc_rs_smow = (1.03091*d18Oc_rs) + 30.91                    ## d18Oc_rs-smow gives re-sampled d18Oc values in smow.
D47_rs = rnorm(z, df[i,5], df[i,6])                           ## D47_rs gives re-sampled D47 values

## 3.4) Gaussian re-sampling of slope and intercept values for D47-temp transfer function
## D47 = 0.0422 (+- 0.0019) *10^6/T^2 + 0.1262 (+- 0.0207)    ## D47_temp transfer function (Bonifacie et al., 2017)
D47_T_tf = c(0.0422, 0.0003, 0.1262, 0.0035)                  ## Update these values if using a different temperature calibration
                                                              ## Note that the 1 sigma stdev is used here.

slope_rs = rnorm(z, D47_T_tf[1], D47_T_tf[2])
int_rs = rnorm(z, D47_T_tf[3], D47_T_tf[4])

## 3.5) Calculate new T_D47 values from re-sampled D47 values, slopes, and intercepts 
## T_D47_rs = ((slope_rs * 10^6) / (D47 - int_rs))^0.5 - 273.15

T_D47_rs = ((slope_rs * 10^6) / (D47_rs - int_rs))^0.5 - 273.15


## 3.6) Calculates new vector of alpha values (for carb-water fractionation) from re-sampled T_D47 values
## 1000ln(a) = 18.03*(10^3/T) - 32.42                         ## temperature-dependent fractionation factor between water and calcite
##                                                                      (Kim and O'Neil, 1997)

alpha_rs = exp((18.03*(10^3/(T_D47_rs+273.15)) - 32.42)/1000)


## 3.7) Calculate new d18Orfw values using re-sampled alpha and d18Oc values

## alpha_rs = (1000 + d18Oc_rs)/(1000 + d18Ofw_rs)  <- solve to isolate d18Orfw_rs
d18Orfw_rs = ((1000 + d18Oc_rs_smow)/alpha_rs)-1000


## 3.8) Calculate re-sampled d18Orfw and T(D47) average and sd values
mean_d18Orfw_rs[i] = round(mean(d18Orfw_rs), digits = 2)
sd_d18Orfw_rs[i] = round(sd(d18Orfw_rs), digits = 2)

mean_T_D47_rs[i] = round(mean(T_D47_rs), digits = 0)
sd_T_D47_rs[i] = round(sd(T_D47_rs), digits = 0)

}

## 4) Compile means and standard deviations; create new data frame; export .csv file

print(mean_d18Orfw_rs)      
print(sd_d18Orfw_rs) 
print(mean_T_D47_rs)
print(sd_T_D47_rs)

df_output = cbind.data.frame(df, mean_d18Orfw_rs, sd_d18Orfw_rs, mean_T_D47_rs, sd_T_D47_rs)    ## Creates a new data frame containing the calculated products from the script

write.csv(df_output, paste0("~/Dropbox (University of Michigan)/Katz - Junin modern paper/R_scripts/df_output ", format(Sys.time(), "%Y-%b-%d %H.%M"), ".csv"))
