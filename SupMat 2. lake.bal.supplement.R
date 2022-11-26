### Title: Detecting hydrologic distinctions among Andean lakes using clumped and triple oxygen isotopes
### Supplementary script: lake.bal.supplement.R
### Author: Sarah Katz
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
###     This script also includes the code to generate Figure 6.


###     Contents:

###     0) Description of script
###     1) Required packages 
###     2) Lake budget steady state equations
###     3) Constants
###     4) User defined model parameters and random selection of parameter values
###     5) Steady state lake water calculations - Eq. 1
###     6) Steady state lake water calculations - Eq. 2
###     7) Plots of modeled lake waters and lambda_lake (used to create Fig. 6)
###     8) Calculate modeled lambda_lake - D'17Ow relationship



#########################################################
#############  1) Required Packages  ####################
#########################################################

## Required packages
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("rgl")


library(ggplot2)
library(ggpubr)
library(rgl)


plot.path <- "~/Desktop/"             ## user may update plot path

#########################################################
##############  2) EQUATIONS  ###########################
#########################################################

## Eq. 1  
##  Rw = (adiff*(1-h)*Ri + aeq*h*Xe*Rv)/(Xe+adiff*(1-h)*(1-Xe))
##  Throughflow lake with evaporation (Criss, 1999, Eq 4.46c)

## Eq. 2
## Rw = (aeq*Ri*(adiff(a-h) + h*(1-F)) + aeq*h*Xe*Rv*F) / (Xe + aeq*(1-Xe)*(adiff*(1-h) + h(1-F))
## Throughflow lake where evaporated water contributes to atmospheric vapor (Benson and White, 1994; Passey and Ji, 2019, Eq. 6)


#########################################################
##############  3) CONSTANTS  ###########################
#########################################################

R18smow = 0.0020052   ## Baertschi, 1976; IAEA Reference sheet
R17smow = 0.0003799   ## Li et al., 1988; IAEA Reference sheet
R2smow = 0.0001558    ## Hagemann et al., 1970; IAEA Reference sheet

theta.eq = 0.529      ## Barkan and Luz, 2005
theta.diff = 0.5185   ## Barkan and Luz, 2007
theta.ref = 0.528

gmwl= c(8,10)         ## Craig, 1961

diffratio18 = 1/0.9723  ## Merlivat, 1978
                ## adiff18 = 1.028489
diffratio2H = 1/0.9755  ## Merlivat, 1978
                ## adiff18 = 1.025115

#########################################################
##############  4) USER DEFINED  ########################
#########################################################

z = 1000        ## number of model simulations

temp = 14       ## Temperature in degrees Celsius


########### RANDOMLY  DISTRIBUTED VARIABLES #############
                 ### Min & Max values ###

## Proportion of turbulent fractionation. May range from 0-1.
Phimin = 0.2     ## MIN 
Phimax = 0.7     ## MAX

## Proportion of Evaporative losses to Inputs. May range from 0 (no evap loss) to 1 (all loss from evap)
Xemin = 0.01    ## MIN
Xemax = 1.0     ## MAX

## Relative humidity normalized to lake surface temperature. May range from 0 (dry atmosphere) to 1 (saturated atmosphere)
hmin = 0.3      ## MIN
hmax = 0.9      ## MAX

## Fraction of local vapor derived from non-lake water sources. May range from 0 (the lake provides 100% of local moisture) to 1 (the lake provides 0% of local moisture)
Fmin = 0.7      ## MIN
Fmax = 1        ## MAX
          

########### Random selection of parameter values between provided min-max values #############

Phi = runif(z, min = Phimin, max = Phimax)
Xe = runif(z, min = Xemin, max = Xemax)
h = runif(z, min = hmin, max = hmax)
Fr = runif(z, min = Fmin, max = Fmax)


##################  User-assigned values for input water  #######################
## Based on amount-weighted mean annual precipitation at Junin, Peru (Katz et al., Section 5.2.1). All in units of per mil.
D17Oi = 0.031
dp18Oi = -14.1
dp17Oi = D17Oi + (dp18Oi * theta.ref)  
d2Hi = -100
dxsi = d2Hi-(8*dp18Oi)

################## Calculate R values for input waters  ##################

Ri18 = exp(dp18Oi/1000)*R18smow
Ri17 = exp(dp17Oi/1000)*R17smow
Ri2H = ((d2Hi/1000)+1)*R2smow

################## Isotopic ratio of water vapor in equilibrium with regional water  ##################

d18Owv = runif(z, min = dp18Oi-2, max = dp18Oi+2)  ## random selection between min and max values for liquid water that atmospheric vapor is in equilibrium with. Water isotopic composition based on amount-weighted mean annual precipitation. 
D17Owv = runif(z, min = 0.015, max = 0.040)        ## random selection between min and max values for liquid water that atmospheric vapor is in equilibrium with. Water isotopic composition based on amount-weighted mean annual precipitation.
dxswv = runif(z, min = 0, max = 30)                ## random selection between min and max values for liquid water that atmospheric vapor is in equilibrium with. Water isotopic composition based on amount-weighted mean annual precipitation.

dp18Owv = log(d18Owv/1000+1)*1000                   ## convert liquid water that atmospheric vapor is in equilibrium with d18O to d'18O
dp17Owv = D17Ov + theta.ref*dp18Owv                 ## calculate liquid water that atmospheric vapor is in equilibrium with d'17O
d2Hwv = dxswv + 8*d18Owv                            ## calculate liquid water that atmospheric vapor is in equilibrium with d2H

## Temperature dependent equilibrium fractionation factor between vapor and liquid water
aeq18 = exp((1137/((temp + 273.15)^2)) - (0.4156/(temp+273.15)) - 0.0020667)    ## Majoube 1971
aeq17 = exp(theta.eq * log(aeq18))                                              
aeq2H = exp((24844/((temp + 273.15)^2)) - (76.248/(temp+273.15)) + 0.052612)    ## Majoube 1971

## Calculate R values for vapor
Rv18 = (exp(dp18Owv/1000)* R18smow)/aeq18
Rv17 = (exp(dp17Owv/1000)*R17smow)/aeq17
Rv2H = (R2smow*((d2Hwv/1000)+1))/aeq2H    

## Isotopic composition of regional vapor
dp18Ovapor = log(Rv18/R18smow)*1000                 ## atmospheric water vapor d'18O
dp17Ovapor = log(Rv17/R17smow)*1000                 ## atmospheric water vapor d'17O
Dp17Ovapor = (dp17Ovapor - (0.528*dp18Ovapor))      ## atmospheric water vapor D'17O. in per mil

d18Ovapor = (Rv18/R18smow - 1)*1000                 ## atmospheric water vapor d18O
d2Hvapor = (Rv2H/R2smow - 1)*1000                   ## atmospheric water vapor d2H
dxsvapor = d2Hvapor - 8*d18Ovapor                   ## atmospheric water vapor d-excess

## Diffusion vs. pure turbulence (i.e. no fractionation). When Phi = 1, all diffusive fractionation; when Phi = 0, no diffusive fractionation (all turbulent)
adiff18 = Phi*diffratio18 + (1-Phi)
adiff17 = exp(theta.diff*log(adiff18))
adiff2H = Phi*diffratio2H + (1-Phi)  


#############################################################################################
#######################  5) LAKE WATER CALCULATIONS - EQ. 1  ################################
##  Assumes evaporated water DOES NOT influence isotopic composition of atmospheric vapor ###
#############################################################################################

## Calculate R values for lake waters
Rw18.eq1 = ((adiff18*(1-h)*Ri18) + (aeq18*h*Xe*Rv18))/(Xe + (adiff18*(1-h)*(1-Xe)))
Rw17.eq1 = ((adiff17*(1-h)*Ri17) + (aeq17*h*Xe*Rv17))/(Xe + (adiff17*(1-h)*(1-Xe)))
Rw2H.eq1 = ((adiff2H*(1-h)*Ri2H) + (aeq2H*h*Xe*Rv2H))/(Xe + (adiff2H*(1-h)*(1-Xe)))

## Calculate delta (d), delta prime (dp), D'17O (Dp), d-excess (dsx) values for lake waters in units of per mil and D'17O (Dp) in units of per meg.
dp18Ow.eq1 = (log(Rw18.eq1/R18smow))*1000
d18Ow.eq1 = ((Rw18.eq1/R18smow)-1)*1000
dp17Ow.eq1 = (log(Rw17.eq1/R17smow))*1000
d2Hw.eq1 = ((Rw2H.eq1/R2smow)-1)*1000
dp2Hw.eq1 = (log(Rw2H.eq1/R2smow))*1000

Dp17Ow.eq1 = (dp17Ow.eq1 - (theta.ref*dp18Ow.eq1))*1000
dxsw.eq1 = d2Hw.eq1 - gmwl[1]*d18Ow.eq1
dpxs.eq1= dp2Hw.eq1 - gmwl[1]*dp18Ow.eq1

lam.lake.eq1 = (dp17Ow.eq1-dp17Oi)/(dp18Ow.eq1-dp18Oi)
EL.eq1 = (d2Hw.eq1- (log(d2Hi/1000+1)*1000) )/(dp18Ow.eq1-dp18Oi)

#############################################################################################
#######################  6) LAKE WATER CALCULATIONS - EQ. 2  ################################
#####  Assumes evaporated water DO influence isotopic composition of atmospheric vapor ######
#############################################################################################

## Calculate R values for lake waters
Rw18.eq2 = ((aeq18*Ri18*(adiff18*(1-h)+h*(1-Fr)))+(aeq18*h*Xe*Rv18*Fr))/
  (Xe+aeq18*(1-Xe)*(adiff18*(1-h)+h*(1-Fr)))
  
Rw17.eq2 = ((aeq17*Ri17*(adiff17*(1-h)+h*(1-Fr)))+(aeq17*h*Xe*Rv17*Fr))/
  (Xe+aeq17*(1-Xe)*(adiff17*(1-h)+h*(1-Fr)))

Rw2H.eq2 = ((aeq2H*Ri2H*(adiff2H*(1-h)+h*(1-Fr)))+(aeq2H*h*Xe*Rv2H*Fr))/
  (Xe+aeq2H*(1-Xe)*(adiff2H*(1-h)+h*(1-Fr)))

## Calculate delta (d), delta prime (dp), d-excess (dsx) values for lake waters in units of per mil and D'17O (Dp) in units of per meg.
dp18Ow.eq2 = (log(Rw18.eq2/R18smow))*1000
d18Ow.eq2 = ((Rw18.eq2/R18smow)-1)*1000
dp17Ow.eq2 = (log(Rw17.eq2/R17smow))*1000
d2Hw.eq2 = ((Rw2H.eq2/R2smow)-1)*1000

Dp17Ow.eq2 = (dp17Ow.eq2 - (theta.ref*dp18Ow.eq2))*1000
dxsw.eq2 = d2Hw.eq2 - gmwl[1]*d18Ow.eq2

## Calculate lambda lake values
lam.lake.eq2 = (dp17Ow.eq2-dp17Oi)/(dp18Ow.eq2-dp18Oi)

#########################################################
#####################  7) PLOTS  ########################
#########################################################
## Requires packages 'ggplot2' and 'ggpubr'

##  Assumes evaporated water does not influence isotopic composition of atmospheric vapor (Eq. 1) ##
##  p1 = D'17O vs d'18O of modeled lake waters using EQ. 1
p1 <- ggplot()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  geom_point(aes(x=dp18Ow.eq1, y = Dp17Ow.eq1), size =2, shape = 23, color="grey")+  
  geom_point(aes(y=D17Oi*1000, x=dp18Oi), size =2, shape = 23, color = "black", stroke = 0.75, fill="white")+ 
  scale_x_continuous(limits = c(-20, 20), expand = c(0, 0))+
  scale_y_continuous(limits = c(-150, 50),expand = c(0, 0), labels = scales::number_format(accuracy = 1))+
  theme(text = element_text(size=16))+
  labs(x=expression(delta*"\u02B9"^"18"*"O (\u2030)"), y=expression(Delta*"\u02B9"^"17"*"O (per meg; "*lambda[ref]*" = 0.528)"))
p1

##  p2 = d-excess vs d'18O of modeled lake waters using EQ. 1
p2 <- ggplot()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  geom_point(aes(x=d18Ow.eq1, y = dxsw.eq1), size =1, shape = 24, color="grey")+  
  geom_point(aes(y=dxsi, x=dp18Oi), size =2, shape = 24, color = "black", stroke = 1, fill="white")+ 
  scale_x_continuous(limits = c(-20, 10), expand = c(0, 0))+
  scale_y_continuous(limits = c(-100, 75),expand = c(0, 0), labels = scales::number_format(accuracy = .1))+
  theme(text = element_text(size=16))+
  labs(x=expression(delta*""^"18"*"O (\u2030)"), y=expression("d-excess (\u2030)"))
p2


#####  Assumes evaporated water do influence influence isotopic composition of atmopspheric vapor (Eq. 2) #####
##  p3 = D'17O vs d'18O of modeled lake waters using EQ. 2
## ** Creates plots for top row of Figure 6
p3 <- ggplot()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  geom_point(aes(x=dp18Ow.eq2, y = Dp17Ow.eq2), size =1.5, shape = 23, color="grey")+  
  geom_point(aes(y=D17Oi*1000, x=dp18Oi), size =2, shape = 23, color = "black", stroke = 0.75, fill="white")+ 
  scale_x_continuous(limits = c(-20, 20), expand = c(0, 0))+
  scale_y_continuous(limits = c(-150, 50),expand = c(0, 0), labels = scales::number_format(accuracy = 1))+
  theme(text = element_text(size=16))+
  labs(x=expression(delta*"\u02B9"^"18"*"O (\u2030)"), y=expression(Delta*"\u02B9"^"17"*"O (per meg; "*lambda[ref]*" = 0.528)"))
p3


##  p4 = d-excess vs d'18O of modeled lake waters using EQ. 2
p4 <- ggplot()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  geom_point(aes(x=d18Ow.eq2, y = dxsw.eq2), size =2, shape = 24, color="grey")+  
  geom_point(aes(y=dxsi, x=dp18Oi), size =2, shape = 24, color = "black", stroke = 1, fill="white")+ 
  scale_x_continuous(limits = c(-20, 10), expand = c(0, 0))+
  scale_y_continuous(limits = c(-100, 75),expand = c(0, 0), labels = scales::number_format(accuracy = .1))+
  theme(text = element_text(size=16))+
  labs(x=expression(delta*""^"18"*"O (\u2030)"), y=expression("d-excess (\u2030)"))
p4

##  p5 = modeled lake water D'17O (Eq 2) vs lambda_lake
## ** Creates plots for bottom row of Figure 6
p5 <- ggplot()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  geom_point(aes(y=lam.lake.eq2, x = Dp17Ow.eq2-D17Oi*1000), size =2, shape = 23, color="grey")+
  geom_smooth(method = 'lm', se=TRUE, aes(y=lam.lake.eq2, x = Dp17Ow.eq2-D17Oi*1000),  formula = y ~ poly(x,3), fill="blue", color="black", linetype=2)+     ## Third-order polynomial fit. Blue envelop is 95% CI
  #geom_smooth(method = 'lm', se=TRUE, aes(y=lam.lake.eq2, x = Dp17Ow.eq2-D17Oi*1000),  formula = y ~ poly(x,2), color="red")+                               ## Second-order polynomial fit. Red envelop is 95% CI
  scale_x_continuous(limits = c(-150, 10), expand = c(0, 0))+
  scale_y_continuous(limits = c(0.51, 0.53),expand = c(0, 0), labels = scales::number_format(accuracy = .001))+
  theme(text = element_text(size=16))+
  labs(y=expression(lambda[lake]), x=expression(Delta*"\u02B9"^"17"*"O"[lake]*" - "*Delta*"\u02B9"^"17"*"O"[input]*" (per meg; "*lambda[ref]*" = 0.528)"))

p5


## Compile p3 and p5
compiled = ggarrange(p3, p5, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
compiled
#ggsave(filename="compiled plots", plot = compiled, path=plot.path, device=cairo_pdf, height=10, width=5 )

#################################################################
#####  8) Calculate modeled lam.lake - D'17Ow relationship  #####
#################################################################

## Using a third order polynomial <-- preferred
Dp17Ow.eq2.permil=Dp17Ow.eq2/1000
model.poly3 = lm(lam.lake.eq2 ~ (Dp17Ow.eq2.permil + I(Dp17Ow.eq2.permil^2) + I(Dp17Ow.eq2.permil^3)))
#summary(model.poly3)

## Using a second order polynomial <-- poorer fit than a third order polynomial
model.poly2 = lm(lam.lake.eq2 ~ (Dp17Ow.eq2.permil + I(Dp17Ow.eq2.permil^2)))



