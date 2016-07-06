###########################################################################
############## RMARK CJS Survival Analysis of Charadrius Chicks ###########
###########################################################################
########################## Luke Eberhart-Phillips #########################
###########################################################################

# Set the working directory
setwd("/home/luke/comparative_ASR/Chick_survival_analysis/") 

# Reference RMark functions and load ggplot functions
library(RMark) 
library(ggplot2)
library(stringr)
library(extrafont)

# Find fonts from computer that are candara or Candara
font_import(pattern="[C/c]andara") 
loadfonts(device = "win") # load these into R

# Adds the location of MARK the RProfile.site file
MarkPath <- "/usr/local/bin/mark"
MarkViewer<-"nano"

# Import prepared capture history file (10 columns: "ch": concatenated daily
# presence/absence, "Ring": unique individual ID, "Year": year that the chick
# hatched ranging from 2009 to 2012, "Hatch_Date_simp": hatch date in the 3 
# character format (i.e., m%d), "Sex": male or female, "Hatch_ID": the nest ID
# from which the chick hatched, "Care_ID" the nest ID from which parental care
# was given (i.e., for cross-fostered chicks, the hatch ID and care IDs will
# be different), "CF": whether or not the chick was cross-fostered (i.e., 1 = 
# yes, 0 = no), "Hatch_brood_size": the number of chicks that hatched from the
# nest where the chick came from, "Care_site": the site at which the chick was
# brooded.
KiP_chicks <- 
  read.table("Data_files/Andava_KiP_chick_capture_history_2009-2015_MARK.txt",
             header=T,colClasses=c("character","factor","factor","integer","numeric",
                                   "factor","factor","integer"))

MP_chicks <- 
  read.table("Data_files/Andava_MP_chick_capture_history_2009-2015_MARK.txt",
             header=T,colClasses=c("character","factor","factor","integer","numeric",
                                   "factor","factor","integer"))

WfP_chicks <- 
  read.table("Data_files/Andava_WfP_chick_capture_history_2009-2015_MARK.txt",
             header=T,colClasses=c("character","factor","factor","integer","numeric",
                                   "factor","factor","integer"))

Maio_chicks <- 
  read.table("Data_files/Maio_chick_capture_history_2007-2015_MARK.txt",
             header=T,colClasses=c("character","factor","integer","numeric",
                                   "factor","factor","factor"))

Tuzla_chicks <- 
  read.table("Data_files/Tuzla_chick_capture_history_1996-2000_MARK.txt",
             header=T,colClasses=c("character","factor","integer","numeric",
                                   "factor","factor"))

# Subset the data to assess survival without chicks from 2008 because encounter
# rates were very low.  Each chick was only encountered once in 2008 because
# siblings were not uniquely ringed and chicks were not re-captured after the
# initial ringing.
#Maio_chicks_no_2008 <- Maio_chicks[which(Maio_chicks$Year != "2008"),]

# Create processed RMARK data format as CJS with 3 groups (sex, and 
# cross-fostering treatment).
KiP_chicks.proc=process.data(KiP_chicks,model="CJS",
                             groups=c("Sex","Year"))
MP_chicks.proc=process.data(MP_chicks,model="CJS",
                            groups=c("Sex","Year"))
WfP_chicks.proc=process.data(WfP_chicks,model="CJS",
                             groups=c("Sex","Year"))
Maio_chicks.proc=process.data(Maio_chicks,model="CJS",
                              groups=c("Sex","Site","Year"))
Tuzla_chicks.proc=process.data(Tuzla_chicks,model="CJS",
                               groups=c("Sex","Year"))
# Create the design data
KiP_chicks.ddl=make.design.data(KiP_chicks.proc)
MP_chicks.ddl=make.design.data(MP_chicks.proc)
WfP_chicks.ddl=make.design.data(WfP_chicks.proc)
Maio_chicks.ddl=make.design.data(Maio_chicks.proc)
Tuzla_chicks.ddl=make.design.data(Tuzla_chicks.proc)

# # check parameter matrices to see if groups were binned correctly
# PIMS(mark(KiP_chicks.proc,KiP_chicks.ddl,model.parameters=list(Phi=list(formula=~Sex+Year)),output=F),"Phi")
# PIMS(mark(MP_chicks.proc,MP_chicks.ddl,model.parameters=list(Phi=list(formula=~Sex+Year)),output=F),"Phi")
# PIMS(mark(WfP_chicks.proc,WfP_chicks.ddl,model.parameters=list(Phi=list(formula=~Sex+Year)),output=F),"Phi")
# PIMS(mark(Maio_chicks.proc,Maio_chicks.ddl,model.parameters=list(Phi=list(formula=~Sex+Site+Year)),output=F),"Phi")
# PIMS(mark(Tuzla_chicks.proc,Tuzla_chicks.ddl,model.parameters=list(Phi=list(formula=~Sex+Year)),output=F),"Phi")

# create a quadratic time variabel so that it can be tested along side the time
# model
time <- c(0:(KiP_chicks.proc$nocc[1]-1))
Quadratic <- time^2
Cubic <- time^3
quad_time <- data.frame(time, Quadratic, Cubic)
KiP_chicks.ddl$p=merge_design.covariates(KiP_chicks.ddl$Phi,quad_time,bygroup=F,bytime=T)
KiP_chicks.ddl$Phi=merge_design.covariates(KiP_chicks.ddl$Phi,quad_time,bygroup=F,bytime=T)

time <- c(0:(MP_chicks.proc$nocc[1]-1))
Quadratic <- time^2
Cubic <- time^3
quad_time <- data.frame(time, Quadratic, Cubic)
MP_chicks.ddl$p=merge_design.covariates(MP_chicks.ddl$Phi,quad_time,bygroup=F,bytime=T)
MP_chicks.ddl$Phi=merge_design.covariates(MP_chicks.ddl$Phi,quad_time,bygroup=F,bytime=T)

time <- c(0:(WfP_chicks.proc$nocc[1]-1))
Quadratic <- time^2
Cubic <- time^3
quad_time <- data.frame(time, Quadratic, Cubic)
WfP_chicks.ddl$p=merge_design.covariates(WfP_chicks.ddl$Phi,quad_time,bygroup=F,bytime=T)
WfP_chicks.ddl$Phi=merge_design.covariates(WfP_chicks.ddl$Phi,quad_time,bygroup=F,bytime=T)

time <- c(0:(Maio_chicks.proc$nocc[1]-1))
Quadratic <- time^2
Cubic <- time^3
quad_time <- data.frame(time, Quadratic, Cubic)
Maio_chicks.ddl$p=merge_design.covariates(Maio_chicks.ddl$Phi,quad_time,bygroup=F,bytime=T)
Maio_chicks.ddl$Phi=merge_design.covariates(Maio_chicks.ddl$Phi,quad_time,bygroup=F,bytime=T)

time <- c(0:(Tuzla_chicks.proc$nocc[1]-1))
Quadratic <- time^2
Cubic <- time^3
quad_time <- data.frame(time, Quadratic, Cubic)
Tuzla_chicks.ddl$p=merge_design.covariates(Tuzla_chicks.ddl$Phi,quad_time,bygroup=F,bytime=T)
Tuzla_chicks.ddl$Phi=merge_design.covariates(Tuzla_chicks.ddl$Phi,quad_time,bygroup=F,bytime=T)
# # Goodness of fit testing on the most general model
# # Define the variable to test and their interactions
# Phi.general.model=list(formula=~Time*Sex*Year*Site) # general Phi model
# p.general.model=list(formula=~Time*Sex*Year*Site) # general p model
# # set wd so that results go to the correct folder
# setwd("H:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/RMARK_model_output/Maio/Chick/General_Model")
# # run model
# general.model=mark(Maio_chicks.proc,Maio_chicks.ddl, 
#                    model.parameters=list(Phi=Phi.general.model,
#                                          p=p.general.model))
# # export results to folder in wd
# export.MARK(Maio_chicks.proc, 
#             project.name="Maio_General_Model",
#             general.model,replace=T) 
# # Open MARK gui, File > "RMark Import", open ".Rinp" file and run median c-hat
# # diagnostics on the model.
# # cleanup old MARK files
# cleanup() 
# 
# # run the program RELEASE GOF test in addition to median c-hat
# release.gof(Maio_chicks_sex.proc)

# Phi.general.model=list(formula=~Time*Sex*Year*Site) # general Phi model
# p.general.model=list(formula=~Time*Sex*Year*Site) # general p model
# # set wd so that results go to the correct folder
# setwd("H:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/RMARK_model_output/Maio/Chick/General_Model")
# # run model
# general.model=mark(Maio_chicks.proc,Maio_chicks.ddl, 
#                    model.parameters=list(Phi=Phi.general.model,
#                                          p=p.general.model))
# # export results to folder in wd
# export.MARK(Maio_chicks.proc, 
#             project.name="Maio_General_Model",
#             general.model,replace=T) 
# # Open MARK gui, File > "RMark Import", open ".Rinp" file and run median c-hat
# # diagnostics on the model.
# # cleanup old MARK files
# cleanup() 
# 
# # run the program RELEASE GOF test in addition to median c-hat
# release.gof(Maio_chicks_sex.proc)

# Setup models to test.
# First assess the top model describing variation in the encounter rate (p)
KiP_chicks_full_analysis=function() 
{
  setwd("/home/luke/comparative_ASR/Chick_survival_analysis/RMark_output_files/KiP") 
  Phi.dot=list(formula=~1) #p=Phi(.)
  Phi.Time=list(formula=~Time) #Phi(T)
  Phi.Quadratic=list(formula=~Quadratic) #Phi(Quadratic)
  Phi.Cubic=list(formula=~Cubic) #Phi(Quadratic)
  Phi.sex=list(formula=~Sex) #p=Phi(sex(.))
  Phi.year=list(formula=~Year)
  Phi.hatchdate=list(formula=~Day_of_Season) #Phi(hatchdate)
  Phi.hatchdate.x.Time=list(formula=~Day_of_Season*Time)
  Phi.hatchdate.x.Quadratic=list(formula=~Day_of_Season*Quadratic)
  Phi.hatchdate.x.Cubic=list(formula=~Day_of_Season*Cubic)
  Phi.hatchdate.x.sex=list(formula=~Day_of_Season*Sex)
  Phi.hatchdate.x.sex.Time=list(formula=~Day_of_Season*Sex*Time)
  Phi.hatchdate.x.sex.Quadratic=list(formula=~Day_of_Season*Sex*Quadratic)
  Phi.hatchdate.x.sex.Cubic=list(formula=~Day_of_Season*Sex*Cubic)
  Phi.hatchdate.Time=list(formula=~Day_of_Season+Time)
  Phi.hatchdate.Quadratic=list(formula=~Day_of_Season+Quadratic)
  Phi.hatchdate.Cubic=list(formula=~Day_of_Season+Cubic)
  Phi.hatchdate.sex=list(formula=~Day_of_Season+Sex)
  Phi.hatchdate.sex.Time=list(formula=~Day_of_Season+Sex+Time)
  Phi.hatchdate.sex.Quadratic=list(formula=~Day_of_Season+Sex+Quadratic)
  Phi.hatchdate.sex.Cubic=list(formula=~Day_of_Season+Sex+Cubic)
  Phi.Time.x.sex=list(formula=~Sex*Time) #Phi(sex(T))
  Phi.Year.x.sex=list(formula=~Sex*Year) #Phi(sex(T))
  Phi.Quadratic.x.sex=list(formula=~Sex*Quadratic) #Phi(sex(Quadratic))
  Phi.Cubic.x.sex=list(formula=~Sex*Cubic) #Phi(sex(Quadratic))
  Phi.Time.sex=list(formula=~Sex+Time) #Phi(sex(T))
  Phi.Year.sex=list(formula=~Sex+Year) #Phi(sex(T))
  Phi.Quadratic.sex=list(formula=~Sex+Quadratic) #Phi(sex(Quadratic))
  Phi.Cubic.sex=list(formula=~Sex+Cubic) #Phi(sex(Cubic))
  Phi.Time.x.sex.x.year=list(formula=~Sex*Time*Year) #Phi(sex(T))
  Phi.Quadratic.x.sex.x.year=list(formula=~Sex*Quadratic*Year) #Phi(sex(Quadratic))
  Phi.Time.sex.year=list(formula=~Sex+Time+Year) #Phi(sex(T))
  Phi.Quadratic.sex.year=list(formula=~Sex+Quadratic+Year) #Phi(sex(Quadratic))
  p.sex=list(formula=~Sex)
  p.sex.Time=list(formula=~Sex+Time)
  p.sex.x.Time=list(formula=~Sex*Time)
  p.sex.Quadratic=list(formula=~Sex+Quadratic)
  p.sex.x.Quadratic=list(formula=~Sex*Quadratic)
  p.sex.x.Year=list(formula=~Sex*Year)
  p.sex.Year=list(formula=~Sex+Year)
  p.sex.Year.Time=list(formula=~Sex+Year+Time)
  p.sex.x.Year.x.Time=list(formula=~Sex*Year*Time)
  p.sex.Year.Quadratic=list(formula=~Sex+Year+Quadratic)
  p.sex.x.Year.x.Quadratic=list(formula=~Sex*Year*Quadratic)
  p.sex.Year.Cubic=list(formula=~Sex+Year+Cubic)
  p.sex.x.Year.x.Cubic=list(formula=~Sex*Year*Cubic)
  p.dot=list(formula=~1) #p(.)
  p.Time=list(formula=~Time)
  p.Quadratic=list(formula=~Quadratic)
  p.Cubic=list(formula=~Cubic)
  p.year=list(formula=~Year)
  p.year.x.Time=list(formula=~Year*Time)
  p.year.x.Quadratic=list(formula=~Year*Quadratic)
  p.year.Time=list(formula=~Year+Time)
  p.year.Quadratic=list(formula=~Year+Quadratic)
  cml=create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
  model.list=mark.wrapper(cml,data=KiP_chicks.proc,ddl=KiP_chicks.ddl,
                          threads = 8, output = FALSE) #runs model list in MARK
  export.MARK(KiP_chicks.proc, project.name="KiP_chick_full_analysis",model.list,replace=T)
  return(model.list) #stores completed model list
}

MP_chicks_full_analysis=function() 
{
  setwd("/home/luke/comparative_ASR/Chick_survival_analysis/RMark_output_files/MP") 
  Phi.dot=list(formula=~1) #p=Phi(.)
  Phi.Time=list(formula=~Time) #Phi(T)
  Phi.Quadratic=list(formula=~Quadratic) #Phi(Quadratic)
  Phi.Cubic=list(formula=~Cubic) #Phi(Quadratic)
  Phi.sex=list(formula=~Sex) #p=Phi(sex(.))
  Phi.year=list(formula=~Year)
  Phi.hatchdate=list(formula=~Day_of_Season) #Phi(hatchdate)
  Phi.hatchdate.x.Time=list(formula=~Day_of_Season*Time)
  Phi.hatchdate.x.Quadratic=list(formula=~Day_of_Season*Quadratic)
  Phi.hatchdate.x.Cubic=list(formula=~Day_of_Season*Cubic)
  Phi.hatchdate.x.sex=list(formula=~Day_of_Season*Sex)
  Phi.hatchdate.x.sex.Time=list(formula=~Day_of_Season*Sex*Time)
  Phi.hatchdate.x.sex.Quadratic=list(formula=~Day_of_Season*Sex*Quadratic)
  Phi.hatchdate.x.sex.Cubic=list(formula=~Day_of_Season*Sex*Cubic)
  Phi.hatchdate.Time=list(formula=~Day_of_Season+Time)
  Phi.hatchdate.Quadratic=list(formula=~Day_of_Season+Quadratic)
  Phi.hatchdate.Cubic=list(formula=~Day_of_Season+Cubic)
  Phi.hatchdate.sex=list(formula=~Day_of_Season+Sex)
  Phi.hatchdate.sex.Time=list(formula=~Day_of_Season+Sex+Time)
  Phi.hatchdate.sex.Quadratic=list(formula=~Day_of_Season+Sex+Quadratic)
  Phi.hatchdate.sex.Cubic=list(formula=~Day_of_Season+Sex+Cubic)
  Phi.Time.x.sex=list(formula=~Sex*Time) #Phi(sex(T))
  Phi.Year.x.sex=list(formula=~Sex*Year) #Phi(sex(T))
  Phi.Quadratic.x.sex=list(formula=~Sex*Quadratic) #Phi(sex(Quadratic))
  Phi.Cubic.x.sex=list(formula=~Sex*Cubic) #Phi(sex(Quadratic))
  Phi.Time.sex=list(formula=~Sex+Time) #Phi(sex(T))
  Phi.Year.sex=list(formula=~Sex+Year) #Phi(sex(T))
  Phi.Quadratic.sex=list(formula=~Sex+Quadratic) #Phi(sex(Quadratic))
  Phi.Cubic.sex=list(formula=~Sex+Cubic) #Phi(sex(Cubic))
  Phi.Time.x.sex.x.year=list(formula=~Sex*Time*Year) #Phi(sex(T))
  Phi.Quadratic.x.sex.x.year=list(formula=~Sex*Quadratic*Year) #Phi(sex(Quadratic))
  Phi.Time.sex.year=list(formula=~Sex+Time+Year) #Phi(sex(T))
  Phi.Quadratic.sex.year=list(formula=~Sex+Quadratic+Year) #Phi(sex(Quadratic))
  p.sex=list(formula=~Sex)
  p.sex.Time=list(formula=~Sex+Time)
  p.sex.x.Time=list(formula=~Sex*Time)
  p.sex.Quadratic=list(formula=~Sex+Quadratic)
  p.sex.x.Quadratic=list(formula=~Sex*Quadratic)
  p.sex.x.Year=list(formula=~Sex*Year)
  p.sex.Year=list(formula=~Sex+Year)
  p.sex.Year.Time=list(formula=~Sex+Year+Time)
  p.sex.x.Year.x.Time=list(formula=~Sex*Year*Time)
  p.sex.Year.Quadratic=list(formula=~Sex+Year+Quadratic)
  p.sex.x.Year.x.Quadratic=list(formula=~Sex*Year*Quadratic)
  p.sex.Year.Cubic=list(formula=~Sex+Year+Cubic)
  p.sex.x.Year.x.Cubic=list(formula=~Sex*Year*Cubic)
  p.dot=list(formula=~1) #p(.)
  p.Time=list(formula=~Time)
  p.Quadratic=list(formula=~Quadratic)
  p.Cubic=list(formula=~Cubic)
  p.year=list(formula=~Year)
  p.year.x.Time=list(formula=~Year*Time)
  p.year.x.Quadratic=list(formula=~Year*Quadratic)
  p.year.Time=list(formula=~Year+Time)
  p.year.Quadratic=list(formula=~Year+Quadratic)
  cml=create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
  model.list=mark.wrapper(cml,data=MP_chicks.proc,ddl=MP_chicks.ddl,
                          threads = 8, output = FALSE) #runs model list in MARK
  export.MARK(MP_chicks.proc, project.name="MP_chick_full_analysis",model.list,replace=T)
  return(model.list) #stores completed model list
}

WfP_chicks_full_analysis=function() 
{
  setwd("/home/luke/comparative_ASR/Chick_survival_analysis/RMark_output_files/WfP") 
  Phi.dot=list(formula=~1) #p=Phi(.)
  Phi.Time=list(formula=~Time) #Phi(T)
  Phi.Quadratic=list(formula=~Quadratic) #Phi(Quadratic)
  Phi.Cubic=list(formula=~Cubic) #Phi(Quadratic)
  Phi.sex=list(formula=~Sex) #p=Phi(sex(.))
  Phi.year=list(formula=~Year)
  Phi.hatchdate=list(formula=~Day_of_Season) #Phi(hatchdate)
  Phi.hatchdate.x.Time=list(formula=~Day_of_Season*Time)
  Phi.hatchdate.x.Quadratic=list(formula=~Day_of_Season*Quadratic)
  Phi.hatchdate.x.Cubic=list(formula=~Day_of_Season*Cubic)
  Phi.hatchdate.x.sex=list(formula=~Day_of_Season*Sex)
  Phi.hatchdate.x.sex.Time=list(formula=~Day_of_Season*Sex*Time)
  Phi.hatchdate.x.sex.Quadratic=list(formula=~Day_of_Season*Sex*Quadratic)
  Phi.hatchdate.x.sex.Cubic=list(formula=~Day_of_Season*Sex*Cubic)
  Phi.hatchdate.Time=list(formula=~Day_of_Season+Time)
  Phi.hatchdate.Quadratic=list(formula=~Day_of_Season+Quadratic)
  Phi.hatchdate.Cubic=list(formula=~Day_of_Season+Cubic)
  Phi.hatchdate.sex=list(formula=~Day_of_Season+Sex)
  Phi.hatchdate.sex.Time=list(formula=~Day_of_Season+Sex+Time)
  Phi.hatchdate.sex.Quadratic=list(formula=~Day_of_Season+Sex+Quadratic)
  Phi.hatchdate.sex.Cubic=list(formula=~Day_of_Season+Sex+Cubic)
  Phi.Time.x.sex=list(formula=~Sex*Time) #Phi(sex(T))
  Phi.Year.x.sex=list(formula=~Sex*Year) #Phi(sex(T))
  Phi.Quadratic.x.sex=list(formula=~Sex*Quadratic) #Phi(sex(Quadratic))
  Phi.Cubic.x.sex=list(formula=~Sex*Cubic) #Phi(sex(Quadratic))
  Phi.Time.sex=list(formula=~Sex+Time) #Phi(sex(T))
  Phi.Year.sex=list(formula=~Sex+Year) #Phi(sex(T))
  Phi.Quadratic.sex=list(formula=~Sex+Quadratic) #Phi(sex(Quadratic))
  Phi.Cubic.sex=list(formula=~Sex+Cubic) #Phi(sex(Cubic))
  Phi.Time.x.sex.x.year=list(formula=~Sex*Time*Year) #Phi(sex(T))
  Phi.Quadratic.x.sex.x.year=list(formula=~Sex*Quadratic*Year) #Phi(sex(Quadratic))
  Phi.Time.sex.year=list(formula=~Sex+Time+Year) #Phi(sex(T))
  Phi.Quadratic.sex.year=list(formula=~Sex+Quadratic+Year) #Phi(sex(Quadratic))
  p.sex=list(formula=~Sex)
  p.sex.Time=list(formula=~Sex+Time)
  p.sex.x.Time=list(formula=~Sex*Time)
  p.sex.Quadratic=list(formula=~Sex+Quadratic)
  p.sex.x.Quadratic=list(formula=~Sex*Quadratic)
  p.sex.x.Year=list(formula=~Sex*Year)
  p.sex.Year=list(formula=~Sex+Year)
  p.sex.Year.Time=list(formula=~Sex+Year+Time)
  p.sex.x.Year.x.Time=list(formula=~Sex*Year*Time)
  p.sex.Year.Quadratic=list(formula=~Sex+Year+Quadratic)
  p.sex.x.Year.x.Quadratic=list(formula=~Sex*Year*Quadratic)
  p.sex.Year.Cubic=list(formula=~Sex+Year+Cubic)
  p.sex.x.Year.x.Cubic=list(formula=~Sex*Year*Cubic)
  p.dot=list(formula=~1) #p(.)
  p.Time=list(formula=~Time)
  p.Quadratic=list(formula=~Quadratic)
  p.Cubic=list(formula=~Cubic)
  p.year=list(formula=~Year)
  p.year.x.Time=list(formula=~Year*Time)
  p.year.x.Quadratic=list(formula=~Year*Quadratic)
  p.year.Time=list(formula=~Year+Time)
  p.year.Quadratic=list(formula=~Year+Quadratic)
  cml=create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
  model.list=mark.wrapper(cml,data=WfP_chicks.proc,ddl=WfP_chicks.ddl,
                          threads = 8, output = FALSE) #runs model list in MARK
  export.MARK(WfP_chicks.proc, project.name="WfP_chick_full_analysis",model.list,replace=T)
  return(model.list) #stores completed model list
}

Maio_chicks_full_analysis=function() 
{
  setwd("/home/luke/comparative_ASR/Chick_survival_analysis/RMark_output_files/Maio") 
  Phi.dot=list(formula=~1) #p=Phi(.)
  Phi.Time=list(formula=~Time) #Phi(T)
  Phi.Quadratic=list(formula=~Quadratic) #Phi(Quadratic)
  Phi.Cubic=list(formula=~Cubic) #Phi(Quadratic)
  Phi.sex=list(formula=~Sex) #p=Phi(sex(.))
  Phi.year=list(formula=~Year)
  Phi.hatchdate=list(formula=~Day_of_Season) #Phi(hatchdate)
  Phi.hatchdate.x.Time=list(formula=~Day_of_Season*Time)
  Phi.hatchdate.x.Quadratic=list(formula=~Day_of_Season*Quadratic)
  Phi.hatchdate.x.Cubic=list(formula=~Day_of_Season*Cubic)
  Phi.hatchdate.x.sex=list(formula=~Day_of_Season*Sex)
  Phi.hatchdate.x.sex.Time=list(formula=~Day_of_Season*Sex*Time)
  Phi.hatchdate.x.sex.Quadratic=list(formula=~Day_of_Season*Sex*Quadratic)
  Phi.hatchdate.x.sex.Cubic=list(formula=~Day_of_Season*Sex*Cubic)
  Phi.hatchdate.Time=list(formula=~Day_of_Season+Time)
  Phi.hatchdate.Quadratic=list(formula=~Day_of_Season+Quadratic)
  Phi.hatchdate.Cubic=list(formula=~Day_of_Season+Cubic)
  Phi.hatchdate.sex=list(formula=~Day_of_Season+Sex)
  Phi.hatchdate.sex.Time=list(formula=~Day_of_Season+Sex+Time)
  Phi.hatchdate.sex.Quadratic=list(formula=~Day_of_Season+Sex+Quadratic)
  Phi.hatchdate.sex.Cubic=list(formula=~Day_of_Season+Sex+Cubic)
  Phi.Time.x.sex=list(formula=~Sex*Time) #Phi(sex(T))
  Phi.Year.x.sex=list(formula=~Sex*Year) #Phi(sex(T))
  Phi.Quadratic.x.sex=list(formula=~Sex*Quadratic) #Phi(sex(Quadratic))
  Phi.Cubic.x.sex=list(formula=~Sex*Cubic) #Phi(sex(Quadratic))
  Phi.Time.sex=list(formula=~Sex+Time) #Phi(sex(T))
  Phi.Year.sex=list(formula=~Sex+Year) #Phi(sex(T))
  Phi.Quadratic.sex=list(formula=~Sex+Quadratic) #Phi(sex(Quadratic))
  Phi.Cubic.sex=list(formula=~Sex+Cubic) #Phi(sex(Cubic))
  Phi.Time.x.sex.x.year=list(formula=~Sex*Time*Year) #Phi(sex(T))
  Phi.Quadratic.x.sex.x.year=list(formula=~Sex*Quadratic*Year) #Phi(sex(Quadratic))
  Phi.Time.sex.year=list(formula=~Sex+Time+Year) #Phi(sex(T))
  Phi.Quadratic.sex.year=list(formula=~Sex+Quadratic+Year) #Phi(sex(Quadratic))
  p.sex=list(formula=~Sex)
  p.sex.Time=list(formula=~Sex+Time)
  p.sex.x.Time=list(formula=~Sex*Time)
  p.sex.Quadratic=list(formula=~Sex+Quadratic)
  p.sex.x.Quadratic=list(formula=~Sex*Quadratic)
  p.sex.x.Year=list(formula=~Sex*Year)
  p.sex.Year=list(formula=~Sex+Year)
  p.sex.Year.Time=list(formula=~Sex+Year+Time)
  p.sex.x.Year.x.Time=list(formula=~Sex*Year*Time)
  p.sex.Year.Quadratic=list(formula=~Sex+Year+Quadratic)
  p.sex.x.Year.x.Quadratic=list(formula=~Sex*Year*Quadratic)
  p.sex.Year.Cubic=list(formula=~Sex+Year+Cubic)
  p.sex.x.Year.x.Cubic=list(formula=~Sex*Year*Cubic)
  p.dot=list(formula=~1) #p(.)
  p.Time=list(formula=~Time)
  p.Quadratic=list(formula=~Quadratic)
  p.Cubic=list(formula=~Cubic)
  p.year=list(formula=~Year)
  p.year.x.Time=list(formula=~Year*Time)
  p.year.x.Quadratic=list(formula=~Year*Quadratic)
  p.year.Time=list(formula=~Year+Time)
  p.year.Quadratic=list(formula=~Year+Quadratic)
  cml=create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
  model.list=mark.wrapper(cml,data=Maio_chicks.proc,ddl=Maio_chicks.ddl,
                          threads = 8, output = FALSE) #runs model list in MARK
  export.MARK(Maio_chicks.proc, project.name="Maio_chicks_full_analysis",model.list,replace=T)
  return(model.list) #stores completed model list
}

Tuzla_chicks_full_analysis=function() 
{
  setwd("/home/luke/comparative_ASR/Chick_survival_analysis/RMark_output_files/Tuzla/") 
  Phi.dot=list(formula=~1) #p=Phi(.)
  Phi.dot=list(formula=~1) #p=Phi(.)
  Phi.Time=list(formula=~Time) #Phi(T)
  Phi.Quadratic=list(formula=~Quadratic) #Phi(Quadratic)
  Phi.Cubic=list(formula=~Cubic) #Phi(Quadratic)
  Phi.sex=list(formula=~Sex) #p=Phi(sex(.))
  Phi.year=list(formula=~Year)
  Phi.hatchdate=list(formula=~Day_of_Season) #Phi(hatchdate)
  Phi.hatchdate.x.Time=list(formula=~Day_of_Season*Time)
  Phi.hatchdate.x.Quadratic=list(formula=~Day_of_Season*Quadratic)
  Phi.hatchdate.x.Cubic=list(formula=~Day_of_Season*Cubic)
  Phi.hatchdate.x.sex=list(formula=~Day_of_Season*Sex)
  Phi.hatchdate.x.sex.Time=list(formula=~Day_of_Season*Sex*Time)
  Phi.hatchdate.x.sex.Quadratic=list(formula=~Day_of_Season*Sex*Quadratic)
  Phi.hatchdate.x.sex.Cubic=list(formula=~Day_of_Season*Sex*Cubic)
  Phi.hatchdate.Time=list(formula=~Day_of_Season+Time)
  Phi.hatchdate.Quadratic=list(formula=~Day_of_Season+Quadratic)
  Phi.hatchdate.Cubic=list(formula=~Day_of_Season+Cubic)
  Phi.hatchdate.sex=list(formula=~Day_of_Season+Sex)
  Phi.hatchdate.sex.Time=list(formula=~Day_of_Season+Sex+Time)
  Phi.hatchdate.sex.Quadratic=list(formula=~Day_of_Season+Sex+Quadratic)
  Phi.hatchdate.sex.Cubic=list(formula=~Day_of_Season+Sex+Cubic)
  Phi.Time.x.sex=list(formula=~Sex*Time) #Phi(sex(T))
  Phi.Year.x.sex=list(formula=~Sex*Year) #Phi(sex(T))
  Phi.Quadratic.x.sex=list(formula=~Sex*Quadratic) #Phi(sex(Quadratic))
  Phi.Cubic.x.sex=list(formula=~Sex*Cubic) #Phi(sex(Quadratic))
  Phi.Time.sex=list(formula=~Sex+Time) #Phi(sex(T))
  Phi.Year.sex=list(formula=~Sex+Year) #Phi(sex(T))
  Phi.Quadratic.sex=list(formula=~Sex+Quadratic) #Phi(sex(Quadratic))
  Phi.Cubic.sex=list(formula=~Sex+Cubic) #Phi(sex(Cubic))
  Phi.Time.x.sex.x.year=list(formula=~Sex*Time*Year) #Phi(sex(T))
  Phi.Quadratic.x.sex.x.year=list(formula=~Sex*Quadratic*Year) #Phi(sex(Quadratic))
  Phi.Time.sex.year=list(formula=~Sex+Time+Year) #Phi(sex(T))
  Phi.Quadratic.sex.year=list(formula=~Sex+Quadratic+Year) #Phi(sex(Quadratic))
  p.sex=list(formula=~Sex)
  p.sex.Time=list(formula=~Sex+Time)
  p.sex.x.Time=list(formula=~Sex*Time)
  p.sex.Quadratic=list(formula=~Sex+Quadratic)
  p.sex.x.Quadratic=list(formula=~Sex*Quadratic)
  p.sex.x.Year=list(formula=~Sex*Year)
  p.sex.Year=list(formula=~Sex+Year)
  p.sex.Year.Time=list(formula=~Sex+Year+Time)
  p.sex.x.Year.x.Time=list(formula=~Sex*Year*Time)
  p.sex.Year.Quadratic=list(formula=~Sex+Year+Quadratic)
  p.sex.x.Year.x.Quadratic=list(formula=~Sex*Year*Quadratic)
  p.sex.Year.Cubic=list(formula=~Sex+Year+Cubic)
  p.sex.x.Year.x.Cubic=list(formula=~Sex*Year*Cubic)
  p.dot=list(formula=~1) #p(.)
  p.Time=list(formula=~Time)
  p.Quadratic=list(formula=~Quadratic)
  p.Cubic=list(formula=~Cubic)
  p.year=list(formula=~Year)
  p.year.x.Time=list(formula=~Year*Time)
  p.year.x.Quadratic=list(formula=~Year*Quadratic)
  p.year.Time=list(formula=~Year+Time)
  p.year.Quadratic=list(formula=~Year+Quadratic)
  cml=create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
  model.list=mark.wrapper(cml,data=Tuzla_chicks.proc,ddl=Tuzla_chicks.ddl,
                          threads = 8, output = FALSE) #runs model list in MARK
  export.MARK(Tuzla_chicks.proc, project.name="Tuzla_chicks_full_analysis",model.list,replace=T)
  return(model.list) #stores completed model list
}

# Run models -- started at 1130 March 9, 2016
KiP_chicks_full_analysis_run <- KiP_chicks_full_analysis()
MP_chicks_full_analysis_run <- MP_chicks_full_analysis()
WfP_chicks_full_analysis_run <- WfP_chicks_full_analysis()
Maio_chicks_full_analysis_run <- Maio_chicks_full_analysis()
Tuzla_chicks_full_analysis_run <- Tuzla_chicks_full_analysis()

write.table(KiP_chicks_full_analysis_run$model.table,
            file = "/home/luke/comparative_ASR/Chick_survival_analysis/Results/KiP_chicks_full_analysis_run_model_table.txt", sep = "\t")
write.table(WfP_chicks_full_analysis_run$model.table,
            file = "/home/luke/comparative_ASR/Chick_survival_analysis/Results/WfP_chicks_full_analysis_run_model_table.txt", sep = "\t")
write.table(MP_chicks_full_analysis_run$model.table,
            file = "/home/luke/comparative_ASR/Chick_survival_analysis/Results/MP_chicks_full_analysis_run_model_table.txt", sep = "\t")
write.table(Maio_chicks_full_analysis_run$model.table,
            file = "/home/luke/comparative_ASR/Chick_survival_analysis/Results/Maio_chicks_full_analysis_run_model_table.txt", sep = "\t")
write.table(Tuzla_chicks_full_analysis_run$model.table,
            file = "/home/luke/comparative_ASR/Chick_survival_analysis/Results/Tuzla_chicks_full_analysis_run_model_table.txt", sep = "\t")

# find models only with the encounter structure "~1" and subset the
# AIC model list. This will reveal the best structure for p.
# KiP, best p structure is 1) p(~Year * Time) or 2) p(~Year * Quadratic)
KiP_not_constant_phi_models <- rownames(KiP_chicks_full_analysis_run$model.table[which(KiP_chicks_full_analysis_run$model.table$Phi != "~1"), ])
KiP_chicks_constant_phi_mod_list <- remove.mark(KiP_chicks_full_analysis_run, as.numeric(KiP_not_constant_phi_models))

# find models only with the encounter structure "~Year * Quadratic" and subset the
# AIC model list. This wil reveal the best structure for phi, while using the
# best structure for p.
# KiP, best for Phi is Phi(~Sex * Time), therefore model 637 is best
KiP_not_Year_x_Quad_models <- rownames(KiP_chicks_full_analysis_run$model.table[which(KiP_chicks_full_analysis_run$model.table$p != "~Year * Quadratic"), ])
KiP_chicks_p_Year_x_Quad_mod_list <- remove.mark(KiP_chicks_full_analysis_run, as.numeric(KiP_not_Year_x_Quad_models))
rownames(KiP_chicks_full_analysis_run$model.table[which(KiP_chicks_full_analysis_run$model.table$Phi == "~Sex * Time" & 
                                                        KiP_chicks_full_analysis_run$model.table$p == "~Year * Quadratic"), ])

write.table(KiP_chicks_p_Year_x_Quad_mod_list$model.table,
            file = "/home/luke/comparative_ASR/Chick_survival_analysis/Results/KiP_chicks_p_Year_x_Quad_mod_list_model_table.txt", sep = "\t")

# WfP, best p structure is 1) p(~Year * Time) or 2) p(~Sex * Time)
WfP_not_constant_phi_models <- rownames(WfP_chicks_full_analysis_run$model.table[which(WfP_chicks_full_analysis_run$model.table$Phi != "~1"), ])
WfP_chicks_constant_phi_mod_list <- remove.mark(WfP_chicks_full_analysis_run, as.numeric(WfP_not_constant_phi_models))

# find models only with the encounter structure "~Year * Time" and subset the
# AIC model list. This wil reveal the best structure for phi, while using the
# best structure for p.
# WfP, best for Phi is Phi(~Sex * Cubic), therefore model 66 is best
WfP_not_Year_x_Time_models <- rownames(WfP_chicks_full_analysis_run$model.table[which(WfP_chicks_full_analysis_run$model.table$p != "~Year * Time"), ])
WfP_chicks_p_Year_x_Time_mod_list <- remove.mark(WfP_chicks_full_analysis_run, as.numeric(WfP_not_Year_x_Time_models))
rownames(WfP_chicks_full_analysis_run$model.table[which(WfP_chicks_full_analysis_run$model.table$Phi == "~Sex * Cubic" & 
                                                          WfP_chicks_full_analysis_run$model.table$p == "~Year * Time"), ])
write.table(WfP_chicks_p_Year_x_Time_mod_list$model.table,
            file = "/home/luke/comparative_ASR/Chick_survival_analysis/Results/WfP_chicks_p_Year_x_Time_mod_list_model_table.txt", sep = "\t")

# MP, best p structure is 1) p(~Year) or 2) p(~Year + Quadratic)
MP_not_constant_phi_models <- rownames(MP_chicks_full_analysis_run$model.table[which(MP_chicks_full_analysis_run$model.table$Phi != "~1"), ])
MP_chicks_constant_phi_mod_list <- remove.mark(MP_chicks_full_analysis_run, as.numeric(MP_not_constant_phi_models))

# find models only with the encounter structure "~Year + Quadratic" and subset the
# AIC model list. This wil reveal the best structure for phi, while using the
# best structure for p.
# MP, best for Phi is Phi(~Sex), therefore model 547 is best
MP_not_Year_Quad_models <- rownames(MP_chicks_full_analysis_run$model.table[which(MP_chicks_full_analysis_run$model.table$p != "~Year + Quadratic"), ])
MP_chicks_p_Year_Quad_mod_list <- remove.mark(MP_chicks_full_analysis_run, as.numeric(MP_not_Year_Quad_models))
rownames(MP_chicks_full_analysis_run$model.table[which(MP_chicks_full_analysis_run$model.table$Phi == "~Sex" & 
                                                          MP_chicks_full_analysis_run$model.table$p == "~Year + Quadratic"), ])
write.table(MP_chicks_p_Year_Quad_mod_list$model.table,
            file = "/home/luke/comparative_ASR/Chick_survival_analysis/Results/WfP_chicks_p_Year_x_Time_mod_list_model_table.txt", sep = "\t")

# Maio, best p structure is 1) p(~Sex * Year * Cubic)
Maio_not_constant_phi_models <- rownames(Maio_chicks_full_analysis_run$model.table[which(Maio_chicks_full_analysis_run$model.table$Phi != "~1"), ])
Maio_chicks_constant_phi_mod_list <- remove.mark(Maio_chicks_full_analysis_run, as.numeric(Maio_not_constant_phi_models))

# find models only with the encounter structure "~Sex * Year * Cubic" and subset the
# AIC model list. This wil reveal the best structure for phi, while using the
# best structure for p.
# Maio, best for Phi is Phi(~Sex * Cubic), therefore model 54 is best
Maio_not_Sex_x_Year_x_Cubic_models <- rownames(Maio_chicks_full_analysis_run$model.table[which(Maio_chicks_full_analysis_run$model.table$p != "~Sex * Year * Cubic"), ])
Maio_chicks_p_Sex_x_Year_x_Cubic_mod_list <- remove.mark(Maio_chicks_full_analysis_run, as.numeric(Maio_not_Sex_x_Year_x_Cubic_models))
rownames(Maio_chicks_full_analysis_run$model.table[which(Maio_chicks_full_analysis_run$model.table$Phi == "~Sex * Cubic" & 
                                                         Maio_chicks_full_analysis_run$model.table$p == "~Sex * Year * Cubic"), ])
write.table(Maio_chicks_p_Sex_x_Year_x_Cubic_mod_list$model.table,
            file = "/home/luke/comparative_ASR/Chick_survival_analysis/Results/Maio_chicks_p_Sex_x_Year_x_Cubic_mod_list_model_table.txt", sep = "\t")

# Tuzla, best p structure is 1) p(~Year + Time) or 2) p(~Year * Time) 
Tuzla_not_constant_phi_models <- rownames(Tuzla_chicks_full_analysis_run$model.table[which(Tuzla_chicks_full_analysis_run$model.table$Phi != "~1"), ])
Tuzla_chicks_constant_phi_mod_list <- remove.mark(Tuzla_chicks_full_analysis_run, as.numeric(Tuzla_not_constant_phi_models))

# find models only with the encounter structure "~Year + Time" and subset the
# AIC model list. This wil reveal the best structure for phi, while using the
# best structure for p.
# Tuzla, best for Phi is Phi(~Sex * Quadratic), therefore model 504 is best
Tuzla_not_Year_Time_models <- rownames(Tuzla_chicks_full_analysis_run$model.table[which(Tuzla_chicks_full_analysis_run$model.table$p != "~Year + Time"), ])
Tuzla_chicks_p_Year_Time_mod_list <- remove.mark(Tuzla_chicks_full_analysis_run, as.numeric(Tuzla_not_Year_Time_models))
rownames(Tuzla_chicks_full_analysis_run$model.table[which(Tuzla_chicks_full_analysis_run$model.table$Phi == "~Sex * Quadratic" & 
                                                           Tuzla_chicks_full_analysis_run$model.table$p == "~Year + Time"), ])
write.table(Tuzla_chicks_p_Year_Time_mod_list$model.table,
            file = "/home/luke/comparative_ASR/Chick_survival_analysis/Results/Tuzla_chicks_p_Year_Time_mod_list_model_table.txt", sep = "\t")

###############################################################################
# KiP
# Plot of Phi(~Sex * Quadratic)p(~Year * Quadratic)
# Assess how chick survival varies by sex and age
# ----------------------------------------------------------------
rownames(KiP_chicks_full_analysis_run$model.table[which(
  KiP_chicks_full_analysis_run$model.table$Phi == "~Sex * Quadratic" & 
    KiP_chicks_full_analysis_run$model.table$p == "~Year * Quadratic"), ])

Sex_Time_reals <- KiP_chicks_full_analysis_run[[637]]$results$real
Groups <- data.frame(str_split_fixed(rownames(Sex_Time_reals), " ", n = 5))
Sex_Time_reals <- cbind(Groups, Sex_Time_reals)
Sex_Time_reals <- Sex_Time_reals[which(Sex_Time_reals$X1 == "Phi"),]
Sex_Time_reals$age <- as.integer(unlist(str_extract_all(Sex_Time_reals$X4,"[0-9]+")))
Sex_Time_reals$Sex <- unlist(str_extract_all(Sex_Time_reals$X2,"[FM]"))
Sex_Time_reals$Sex <- as.factor(ifelse(Sex_Time_reals$Sex == "F","Female","Male"))
Sex_Time_reals$Year <- unlist(substr(Sex_Time_reals$X2, 5, 8))
Sex_Time_plot <- ggplot(Sex_Time_reals, aes(x = age, y = estimate, group = Sex)) + 
  theme_bw() +
  geom_line(size = 1.5, aes(colour = Sex)) +                                      
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) + #+
  theme(text=element_text(size=16, family="Candara"), # set the font as Candara
        legend.position = c(0.25, 0.75), 
        legend.justification = c(1, 0), panel.grid = element_blank(), 
        
        legend.text=element_text(size=30),
        legend.title=element_text(size=30),
        axis.title.x = element_text(size=35, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=35, vjust=1.2),
        axis.text.y  = element_text(size=30), 
        plot.title = element_text(vjust=1.2, size=20, face="bold"),
        panel.grid.major = element_blank()) +
  
  #legend.title = element_text(size=17),
  #legend.text  = element_text(size=17),
  #axis.title.x = element_text(size=20, vjust=-0.3),
  #axis.text.x  = element_text(size=17), 
  #axis.title.y = element_text(size=20, vjust=1.2),
  #axis.text.y  = element_text(size=17), 
  #plot.title = element_text(vjust=1.2, size=20, face="bold"),
  #panel.background = element_rect(colour = "black", size=1.2)) +
  scale_colour_brewer(palette = "Set1") + 
  xlab("Age (days)") + 
  ylab("Estimated daily survival rate (± 95% CI)")# +
#ggtitle("Figure 3")
Sex_Time_plot

Survival_to_Fledge_F <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("estimate")])
Survival_to_Fledge_F_lcl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("lcl")])
Survival_to_Fledge_F_ucl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("ucl")])
Survival_to_Fledge_M <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("estimate")])
Survival_to_Fledge_M_lcl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("lcl")])
Survival_to_Fledge_M_ucl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("ucl")])
estimate <- c(Survival_to_Fledge_F, Survival_to_Fledge_M)
lcl <- c(Survival_to_Fledge_F_lcl, Survival_to_Fledge_M_lcl)
ucl <- c(Survival_to_Fledge_F_ucl, Survival_to_Fledge_M_ucl)
Sex <- c("Female", "Male")
age <- c("Chick", "Chick")
KiP_Sex_chick_survival <- data.frame(Sex, age, estimate, lcl, ucl)
KiP_Sex_chick_survival$species <- "Kittlitz's"
###############################################################################
# WfP
# Plot of Phi(~Sex * Quadratic)p(~Year * Quadratic)
# Assess how chick survival varies by sex and age
# ----------------------------------------------------------------
Sex_Time_reals <- WfP_chicks_full_analysis_run[[66]]$results$real
Groups <- data.frame(str_split_fixed(rownames(Sex_Time_reals), " ", n = 5))
Sex_Time_reals <- cbind(Groups, Sex_Time_reals)
Sex_Time_reals <- Sex_Time_reals[which(Sex_Time_reals$X1 == "Phi"),]
Sex_Time_reals$age <- as.integer(unlist(str_extract_all(Sex_Time_reals$X4,"[0-9]+")))
Sex_Time_reals$Sex <- unlist(str_extract_all(Sex_Time_reals$X2,"[FM]"))
Sex_Time_reals$Sex <- as.factor(ifelse(Sex_Time_reals$Sex == "F","Female","Male"))
Sex_Time_reals$Year <- unlist(substr(Sex_Time_reals$X2, 5, 8))
Sex_Time_plot <- ggplot(Sex_Time_reals, aes(x = age, y = estimate, group = Sex)) + 
  theme_bw() +
  geom_line(size = 1.5, aes(colour = Sex)) +                                      
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) + #+
  theme(text=element_text(size=16, family="Candara"), # set the font as Candara
        legend.position = c(0.25, 0.75), 
        legend.justification = c(1, 0), panel.grid = element_blank(), 
        
        legend.text=element_text(size=30),
        legend.title=element_text(size=30),
        axis.title.x = element_text(size=35, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=35, vjust=1.2),
        axis.text.y  = element_text(size=30), 
        plot.title = element_text(vjust=1.2, size=20, face="bold"),
        panel.grid.major = element_blank()) +
  
  #legend.title = element_text(size=17),
  #legend.text  = element_text(size=17),
  #axis.title.x = element_text(size=20, vjust=-0.3),
  #axis.text.x  = element_text(size=17), 
  #axis.title.y = element_text(size=20, vjust=1.2),
  #axis.text.y  = element_text(size=17), 
  #plot.title = element_text(vjust=1.2, size=20, face="bold"),
  #panel.background = element_rect(colour = "black", size=1.2)) +
  scale_colour_brewer(palette = "Set1") + 
  xlab("Age (days)") + 
  ylab("Estimated daily survival rate (± 95% CI)")# +
#ggtitle("Figure 3")
Sex_Time_plot

Survival_to_Fledge_F <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("estimate")])
Survival_to_Fledge_F_lcl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("lcl")])
Survival_to_Fledge_F_ucl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("ucl")])
Survival_to_Fledge_M <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("estimate")])
Survival_to_Fledge_M_lcl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("lcl")])
Survival_to_Fledge_M_ucl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("ucl")])
estimate <- c(Survival_to_Fledge_F, Survival_to_Fledge_M)
lcl <- c(Survival_to_Fledge_F_lcl, Survival_to_Fledge_M_lcl)
ucl <- c(Survival_to_Fledge_F_ucl, Survival_to_Fledge_M_ucl)
Sex <- c("Female", "Male")
age <- c("Chick", "Chick")
WfP_Sex_chick_survival <- data.frame(Sex, age, estimate, lcl, ucl)
WfP_Sex_chick_survival$species <- "White-fronted"
###############################################################################
# MP
# Plot of Phi(~Sex * Quadratic)p(~Year * Quadratic)
# Assess how chick survival varies by sex and age
# ----------------------------------------------------------------
Sex_Time_reals <- MP_chicks_full_analysis_run[[547]]$results$real
Groups <- data.frame(str_split_fixed(rownames(Sex_Time_reals), " ", n = 5))
Sex_Time_reals <- cbind(Groups, Sex_Time_reals)
Sex_Time_reals <- Sex_Time_reals[which(Sex_Time_reals$X1 == "Phi"),]
Sex_Time_reals$age <- as.integer(unlist(str_extract_all(Sex_Time_reals$X4,"[0-9]+")))
Sex_Time_reals$Sex <- unlist(str_extract_all(Sex_Time_reals$X2,"[FM]"))
Sex_Time_reals$Sex <- as.factor(ifelse(Sex_Time_reals$Sex == "F","Female","Male"))
Sex_Time_reals$Year <- unlist(substr(Sex_Time_reals$X2, 5, 8))
Sex_Time_plot <- ggplot(Sex_Time_reals, aes(x = age, y = estimate, group = Sex)) + 
  theme_bw() +
  geom_line(size = 1.5, aes(colour = Sex)) +                                      
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) + #+
  theme(text=element_text(size=16, family="Candara"), # set the font as Candara
        legend.position = c(0.25, 0.75), 
        legend.justification = c(1, 0), panel.grid = element_blank(), 
        
        legend.text=element_text(size=30),
        legend.title=element_text(size=30),
        axis.title.x = element_text(size=35, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=35, vjust=1.2),
        axis.text.y  = element_text(size=30), 
        plot.title = element_text(vjust=1.2, size=20, face="bold"),
        panel.grid.major = element_blank()) +
  
  #legend.title = element_text(size=17),
  #legend.text  = element_text(size=17),
  #axis.title.x = element_text(size=20, vjust=-0.3),
  #axis.text.x  = element_text(size=17), 
  #axis.title.y = element_text(size=20, vjust=1.2),
  #axis.text.y  = element_text(size=17), 
  #plot.title = element_text(vjust=1.2, size=20, face="bold"),
  #panel.background = element_rect(colour = "black", size=1.2)) +
  scale_colour_brewer(palette = "Set1") + 
  xlab("Age (days)") + 
  ylab("Estimated daily survival rate (± 95% CI)")# +
#ggtitle("Figure 3")
Sex_Time_plot

Survival_to_Fledge_F <- 
  Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("estimate")]^30
Survival_to_Fledge_F_lcl <- 
  Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("lcl")]^30
Survival_to_Fledge_F_ucl <- 
  Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("ucl")]^30
Survival_to_Fledge_M <- 
  Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("estimate")]^30
Survival_to_Fledge_M_lcl <- 
  Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("lcl")]^30
Survival_to_Fledge_M_ucl <- 
  Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("ucl")]^30
estimate <- c(Survival_to_Fledge_F, Survival_to_Fledge_M)
lcl <- c(Survival_to_Fledge_F_lcl, Survival_to_Fledge_M_lcl)
ucl <- c(Survival_to_Fledge_F_ucl, Survival_to_Fledge_M_ucl)
Sex <- c("Female", "Male")
age <- c("Chick", "Chick")
MP_Sex_chick_survival <- data.frame(Sex, age, estimate, lcl, ucl)
MP_Sex_chick_survival$species <- "Madagascar"
###############################################################################
# Maio
# Plot of Phi(~Sex * Quadratic)p(~Year * Quadratic)
# Assess how chick survival varies by sex and age
# ----------------------------------------------------------------
Sex_Time_reals <- Maio_chicks_full_analysis_run[[54]]$results$real
Groups <- data.frame(str_split_fixed(rownames(Sex_Time_reals), " ", n = 5))
Sex_Time_reals <- cbind(Groups, Sex_Time_reals)
Sex_Time_reals <- Sex_Time_reals[which(Sex_Time_reals$X1 == "Phi"),]
Sex_Time_reals$age <- as.integer(unlist(str_extract_all(Sex_Time_reals$X4,"[0-9]+")))
Sex_Time_reals$Sex <- unlist(str_extract_all(Sex_Time_reals$X2,"[FM]"))
Sex_Time_reals$Sex <- as.factor(ifelse(Sex_Time_reals$Sex == "F","Female","Male"))
Sex_Time_reals$Year <- unlist(substr(Sex_Time_reals$X2, 5, 8))
Sex_Time_plot <- ggplot(Sex_Time_reals, aes(x = age, y = estimate, group = Sex)) + 
  theme_bw() +
  geom_line(size = 1.5, aes(colour = Sex)) +                                      
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) + #+
  theme(text=element_text(size=16, family="Candara"), # set the font as Candara
        legend.position = c(0.25, 0.75), 
        legend.justification = c(1, 0), panel.grid = element_blank(), 
        
        legend.text=element_text(size=30),
        legend.title=element_text(size=30),
        axis.title.x = element_text(size=35, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=35, vjust=1.2),
        axis.text.y  = element_text(size=30), 
        plot.title = element_text(vjust=1.2, size=20, face="bold"),
        panel.grid.major = element_blank()) +
  
  #legend.title = element_text(size=17),
  #legend.text  = element_text(size=17),
  #axis.title.x = element_text(size=20, vjust=-0.3),
  #axis.text.x  = element_text(size=17), 
  #axis.title.y = element_text(size=20, vjust=1.2),
  #axis.text.y  = element_text(size=17), 
  #plot.title = element_text(vjust=1.2, size=20, face="bold"),
  #panel.background = element_rect(colour = "black", size=1.2)) +
  scale_colour_brewer(palette = "Set1") + 
  xlab("Age (days)") + 
  ylab("Estimated daily survival rate (± 95% CI)")# +
#ggtitle("Figure 3")
Sex_Time_plot

Survival_to_Fledge_F <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("estimate")])
Survival_to_Fledge_F_lcl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("lcl")])
Survival_to_Fledge_F_ucl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("ucl")])
Survival_to_Fledge_M <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("estimate")])
Survival_to_Fledge_M_lcl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("lcl")])
Survival_to_Fledge_M_ucl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("ucl")])
estimate <- c(Survival_to_Fledge_F, Survival_to_Fledge_M)
lcl <- c(Survival_to_Fledge_F_lcl, Survival_to_Fledge_M_lcl)
ucl <- c(Survival_to_Fledge_F_ucl, Survival_to_Fledge_M_ucl)
Sex <- c("Female", "Male")
age <- c("Chick", "Chick")
Maio_Sex_chick_survival <- data.frame(Sex, age, estimate, lcl, ucl)
Maio_Sex_chick_survival$species <- "Kentish (Maio)"
###############################################################################
# Tuzla
# Plot of Phi(~Sex * Quadratic)p(~Year * Quadratic)
# Assess how chick survival varies by sex and age
# ----------------------------------------------------------------
Sex_Time_reals <- Tuzla_chicks_full_analysis_run[[504]]$results$real
Groups <- data.frame(str_split_fixed(rownames(Sex_Time_reals), " ", n = 5))
Sex_Time_reals <- cbind(Groups, Sex_Time_reals)
Sex_Time_reals <- Sex_Time_reals[which(Sex_Time_reals$X1 == "Phi"),]
Sex_Time_reals$age <- as.integer(unlist(str_extract_all(Sex_Time_reals$X4,"[0-9]+")))
Sex_Time_reals$Sex <- unlist(str_extract_all(Sex_Time_reals$X2,"[FM]"))
Sex_Time_reals$Sex <- as.factor(ifelse(Sex_Time_reals$Sex == "F","Female","Male"))
Sex_Time_reals$Year <- unlist(substr(Sex_Time_reals$X2, 5, 8))
Sex_Time_plot <- ggplot(Sex_Time_reals, aes(x = age, y = estimate, group = Sex)) + 
  theme_bw() +
  geom_line(size = 1.5, aes(colour = Sex)) +                                      
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) + #+
  theme(text=element_text(size=16, family="Candara"), # set the font as Candara
        legend.position = c(0.25, 0.75), 
        legend.justification = c(1, 0), panel.grid = element_blank(), 
        
        legend.text=element_text(size=30),
        legend.title=element_text(size=30),
        axis.title.x = element_text(size=35, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=35, vjust=1.2),
        axis.text.y  = element_text(size=30), 
        plot.title = element_text(vjust=1.2, size=20, face="bold"),
        panel.grid.major = element_blank()) +
  
  #legend.title = element_text(size=17),
  #legend.text  = element_text(size=17),
  #axis.title.x = element_text(size=20, vjust=-0.3),
  #axis.text.x  = element_text(size=17), 
  #axis.title.y = element_text(size=20, vjust=1.2),
  #axis.text.y  = element_text(size=17), 
  #plot.title = element_text(vjust=1.2, size=20, face="bold"),
  #panel.background = element_rect(colour = "black", size=1.2)) +
  scale_colour_brewer(palette = "Set1") + 
  xlab("Age (days)") + 
  ylab("Estimated daily survival rate (± 95% CI)")# +
#ggtitle("Figure 3")
Sex_Time_plot

Survival_to_Fledge_F <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("estimate")])
Survival_to_Fledge_F_lcl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("lcl")])
Survival_to_Fledge_F_ucl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Female"),
                      c("ucl")])
Survival_to_Fledge_M <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("estimate")])
Survival_to_Fledge_M_lcl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("lcl")])
Survival_to_Fledge_M_ucl <- 
  prod(Sex_Time_reals[which(Sex_Time_reals$Sex == "Male"),
                      c("ucl")])
estimate <- c(Survival_to_Fledge_F, Survival_to_Fledge_M)
lcl <- c(Survival_to_Fledge_F_lcl, Survival_to_Fledge_M_lcl)
ucl <- c(Survival_to_Fledge_F_ucl, Survival_to_Fledge_M_ucl)
Sex <- c("Female", "Male")
age <- c("Chick", "Chick")
Tuzla_Sex_chick_survival <- data.frame(Sex, age, estimate, lcl, ucl)
Tuzla_Sex_chick_survival$species <- "Kentish (Tuzla)"

# Join all chick estimates together
All_pops_sex_chick_survival <- rbind(KiP_Sex_chick_survival,
                                     WfP_Sex_chick_survival,
                                     MP_Sex_chick_survival,
                                     Maio_Sex_chick_survival,
                                     Tuzla_Sex_chick_survival)

write.table(All_pops_sex_chick_survival,
            file = "/home/luke/comparative_ASR/Chick_survival_analysis/Results/All_pops_sex_chick_survival.txt", sep = "\t")