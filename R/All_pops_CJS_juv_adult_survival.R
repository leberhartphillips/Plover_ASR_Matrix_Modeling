#-----------------------------------------------------------------------------#
#       RMark CJS Survival Analysis of Charadrius Juveniles and Adults        #
#                Includes 2015 data from Madagascar and Maio                  #
#                         Luke Eberhart-Phillips                              #
#-----------------------------------------------------------------------------#

# Set the working directory
setwd("/home/luke/comparative_ASR/Juv_Adult_survival_analysis/") 

# Reference RMark functions and other packages needed
library(RMark) 
library(stringr)
library(ggplot2)
library(stringr)
library(extrafont)

# Find fonts from computer that are candara or Candara
font_import(pattern="[C/c]andara") 
loadfonts(device = "win") # load these into R

# Adds the location of MARK the RProfile.site file
MarkPath <- "/usr/local/bin/mark"
MarkViewer<-"nano"

# Import prepared capture history file (three columns: "ch" ... concatenated annual presence/absence, 
# "sex"... male or female, or NA for unknown sex of chicks, "age"... Adult or juvenile [NOTE: this describes 
# the stage that the individual was initially ringed])
Andava_KiP <- read.table("Data_files/Andava_KiP_juv_adult_capture_history_2009-2015_MARK.txt",header=T,colClasses=c("character","factor","factor","factor"))
Andava_MP <- read.table("Data_files/Andava_MP_juv_adult_capture_history_2009-2015_MARK.txt",header=T,colClasses=c("character","factor","factor","factor"))
Andava_WfP <- read.table("Data_files/Andava_WfP_juv_adult_capture_history_2009-2015_MARK.txt",header=T,colClasses=c("character","factor","factor","factor"))
Tuzla_kp <- read.table("Data_files/Tuzla_juv_adult_capture_history_1996-2000_MARK.txt",header=T,colClasses=c("character","factor","factor"))
Maio_kp <- read.table("Data_files/Maio_juv_adult_capture_history_2009-2015_MARK.txt",header=T,colClasses=c("character","factor","factor"))

# Create processed RMARK data format as CJS with 2 groups (sex and age initally ringed) and starting at year 2007
Andava_KiP.proc=process.data(Andava_KiP,model="CJS",groups=c("sex","age"),begin.time=2009,age.var=2,initial.age=c(1,0))
Andava_MP.proc=process.data(Andava_MP,model="CJS",groups=c("sex","age"),begin.time=2009,age.var=2,initial.age=c(1,0))
Andava_WfP.proc=process.data(Andava_WfP,model="CJS",groups=c("sex","age"),begin.time=2009,age.var=2,initial.age=c(1,0))
Tuzla_kp.proc=process.data(Tuzla_kp,model="CJS",groups=c("sex","age"),begin.time=1996,age.var=2,initial.age=c(1,0))
Maio_kp.proc=process.data(Maio_kp,model="CJS",groups=c("sex","age"),begin.time=2007,age.var=2,initial.age=c(1,0))

# Create the design data
Andava_KiP.ddl=make.design.data(Andava_KiP.proc)
Andava_MP.ddl=make.design.data(Andava_MP.proc)
Andava_WfP.ddl=make.design.data(Andava_WfP.proc)
Tuzla_kp.ddl=make.design.data(Tuzla_kp.proc)
Maio_kp.ddl=make.design.data(Maio_kp.proc)

# adds firstyear/adult age field to design data in column "age"
Andava_KiP.ddl=add.design.data(Andava_KiP.proc,Andava_KiP.ddl,"Phi","age",bins=c(0,1,6),right=F,name="age",replace=T)
Andava_MP.ddl=add.design.data(Andava_MP.proc,Andava_MP.ddl,"Phi","age",bins=c(0,1,6),right=F,name="age",replace=T)
Andava_WfP.ddl=add.design.data(Andava_WfP.proc,Andava_WfP.ddl,"Phi","age",bins=c(0,1,6),right=F,name="age",replace=T)
Tuzla_kp.ddl=add.design.data(Tuzla_kp.proc,Tuzla_kp.ddl,"Phi","age",bins=c(0,1,6),right=F,name="age",replace=T)
Maio_kp.ddl=add.design.data(Maio_kp.proc,Maio_kp.ddl,"Phi","age",bins=c(0,1,8),right=F,name="age",replace=T)

# create a dummy field called marked.as.adult which is 0 for the group initally ringed as juvenile and 1 for the group marked as adults.
Andava_KiP.ddl$Phi$marked.as.adult=0
Andava_KiP.ddl$Phi$marked.as.adult[Andava_KiP.ddl$Phi$initial.age.class=="A"]=1 
Andava_KiP.ddl$p$marked.as.adult=0
Andava_KiP.ddl$p$marked.as.adult[Andava_KiP.ddl$p$initial.age.class=="A"]=1

Andava_MP.ddl$Phi$marked.as.adult=0
Andava_MP.ddl$Phi$marked.as.adult[Andava_MP.ddl$Phi$initial.age.class=="A"]=1 
Andava_MP.ddl$p$marked.as.adult=0
Andava_MP.ddl$p$marked.as.adult[Andava_MP.ddl$p$initial.age.class=="A"]=1

Andava_WfP.ddl$Phi$marked.as.adult=0
Andava_WfP.ddl$Phi$marked.as.adult[Andava_WfP.ddl$Phi$initial.age.class=="A"]=1 
Andava_WfP.ddl$p$marked.as.adult=0
Andava_WfP.ddl$p$marked.as.adult[Andava_WfP.ddl$p$initial.age.class=="A"]=1

Tuzla_kp.ddl$Phi$marked.as.adult=0
Tuzla_kp.ddl$Phi$marked.as.adult[Tuzla_kp.ddl$Phi$initial.age.class=="A"]=1 
Tuzla_kp.ddl$p$marked.as.adult=0
Tuzla_kp.ddl$p$marked.as.adult[Tuzla_kp.ddl$p$initial.age.class=="A"]=1

Maio_kp.ddl$Phi$marked.as.adult=0
Maio_kp.ddl$Phi$marked.as.adult[Maio_kp.ddl$Phi$initial.age.class=="A"]=1 
Maio_kp.ddl$p$marked.as.adult=0
Maio_kp.ddl$p$marked.as.adult[Maio_kp.ddl$p$initial.age.class=="A"]=1

# check parameter matrices to see if groups were binned correctly
# PIMS(mark(Andava_KiP.proc,Andava_KiP.ddl,model.parameters=list(Phi=list(formula=~age+sex)),output=F),"Phi")
# PIMS(mark(Andava_MP.proc,Andava_MP.ddl,model.parameters=list(Phi=list(formula=~age+sex)),output=F),"Phi")
# PIMS(mark(Andava_WfP.proc,Andava_WfP.ddl,model.parameters=list(Phi=list(formula=~age+sex)),output=F),"Phi")
# PIMS(mark(Tuzla_kp.proc,Tuzla_kp.ddl,model.parameters=list(Phi=list(formula=~age+sex)),output=F),"Phi")
# PIMS(mark(Maio_kp.proc,Maio_kp.ddl,model.parameters=list(Phi=list(formula=~age+sex)),output=F),"Phi")

# import quadratic time file so that it can be tested along side the annual models
time <- c(0:(Andava_KiP.proc$nocc[1]-1))
Quadratic <- time^2
Cubic <- time^3
quad_time <- data.frame(time, Quadratic, Cubic)
quad_time$time <- c(2009:2015)
Andava_KiP.ddl$p=merge_design.covariates(Andava_KiP.ddl$Phi,quad_time,bygroup=F,bytime=T)
Andava_KiP.ddl$Phi=merge_design.covariates(Andava_KiP.ddl$Phi,quad_time,bygroup=F,bytime=T)

time <- c(0:(Andava_MP.proc$nocc[1]-1))
Quadratic <- time^2
Cubic <- time^3
quad_time <- data.frame(time, Quadratic, Cubic)
quad_time$time <- c(2009:2015)
Andava_MP.ddl$p=merge_design.covariates(Andava_MP.ddl$Phi,quad_time,bygroup=F,bytime=T)
Andava_MP.ddl$Phi=merge_design.covariates(Andava_MP.ddl$Phi,quad_time,bygroup=F,bytime=T)

time <- c(0:(Andava_WfP.proc$nocc[1]-1))
Quadratic <- time^2
Cubic <- time^3
quad_time <- data.frame(time, Quadratic, Cubic)
quad_time$time <- c(2009:2015)
Andava_WfP.ddl$p=merge_design.covariates(Andava_WfP.ddl$Phi,quad_time,bygroup=F,bytime=T)
Andava_WfP.ddl$Phi=merge_design.covariates(Andava_WfP.ddl$Phi,quad_time,bygroup=F,bytime=T)

time <- c(0:(Tuzla_kp.proc$nocc[1]-1))
Quadratic <- time^2
Cubic <- time^3
quad_time <- data.frame(time, Quadratic, Cubic)
quad_time$time <- c(1996:2000,2014)
Tuzla_kp.ddl$p=merge_design.covariates(Tuzla_kp.ddl$Phi,quad_time,bygroup=F,bytime=T)
Tuzla_kp.ddl$Phi=merge_design.covariates(Tuzla_kp.ddl$Phi,quad_time,bygroup=F,bytime=T)

time <- c(0:(Maio_kp.proc$nocc[1]-1))
Quadratic <- time^2
Cubic <- time^3
quad_time <- data.frame(time, Quadratic, Cubic)
quad_time$time <- c(2007:2015)
Maio_kp.ddl$p=merge_design.covariates(Maio_kp.ddl$Phi,quad_time,bygroup=F,bytime=T)
Maio_kp.ddl$Phi=merge_design.covariates(Maio_kp.ddl$Phi,quad_time,bygroup=F,bytime=T)

# # Goodness of fit testing on the most general model
# Phi.sex.age.time=list(formula=~sex*age*time) # general Phi model
# p.sex.age.time=list(formula=~age*sex*time) # general p model
# setwd("H:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/RMARK_model_output/Maio/General_Model") # set wd so that results go to the correct folder
# general.model=mark(Maio_kp.proc,Maio_kp.ddl, model.parameters=list(Phi=Phi.sex.age.time,p=p.sex.age.time)) # run model
# export.MARK(Maio_kp.proc, project.name="Maio_General_Model",general.model,replace=T) # export results to folder in wd
# # Open MARK gui, File > "RMark Import", open ".Rinp" file and run median c-hat diagnostics on the model.
# cleanup() # cleanup old MARK files
# 
# release.gof(Maio_kp.proc)

# Create function of a priori models

KiP_juv_adult_full_analysis=function() 
{
  setwd("/home/luke/comparative_ASR/Juv_Adult_survival_analysis/RMark_output_files/KiP/") 
  Phi.age.time=list(formula=~age+time)
  Phi.sex.time=list(formula=~sex+time)
  Phi.sex.age.time=list(formula=~sex+age+time)
  Phi.sexxtime=list(formula=~sex*time)
  Phi.agextime=list(formula=~age*time)
  Phi.agexsexxtime=list(formula=~age*sex*time)
  Phi.age.Quadratic=list(formula=~age+Quadratic)
  Phi.sex.Quadratic=list(formula=~sex+Quadratic)
  Phi.sexxQuadratic=list(formula=~sex*Quadratic)
  Phi.agexQuadratic=list(formula=~age*Quadratic)
  Phi.sex=list(formula=~sex)
  Phi.age=list(formula=~age)
  Phi.age.sex=list(formula=~age+sex) 
  Phi.agexsex=list(formula=~age*sex) 
  Phi.dot=list(formula=~1) 
  Phi.time=list(formula=~time)
  Phi.Time=list(formula=~Time)
  Phi.Quadratic=list(formula=~Quadratic)
  Phi.age.Time=list(formula=~age+Time)
  Phi.sex.Time=list(formula=~sex+Time)
  Phi.sexxTime=list(formula=~sex*Time)
  Phi.agexTime=list(formula=~age*Time)
  Phi.sex.age.Time=list(formula=~sex+age+Time)
  Phi.sexXageXTime=list(formula=~sex*age*Time)
  p.sex=list(formula=~sex)
  p.age=list(formula=~age) 
  p.dot=list(formula=~1)
  p.time=list(formula=~time)
  p.Time=list(formula=~Time)
  p.Quadratic=list(formula=~Quadratic)
  p.sexxtime=list(formula=~sex*time)
  p.agextime=list(formula=~age*time)
  p.sexxTime=list(formula=~sex*Time)
  p.agexTime=list(formula=~age*Time)
  p.agexsex=list(formula=~age*sex)
  p.Quadraticxsex=list(formula=~Quadratic*sex)
  p.Quadraticxage=list(formula=~Quadratic*age)
  p.Quadraticxagexsex=list(formula=~Quadratic*age*sex)
  p.Timexagexsex=list(formula=~Time*age*sex)
  p.timexagexsex=list(formula=~time*age*sex)
  p.sex_time=list(formula=~sex+time)
  p.age_time=list(formula=~age+time)
  p.sex_Time=list(formula=~sex+Time)
  p.age_Time=list(formula=~age+Time)
  p.age_sex=list(formula=~age+sex)
  p.Quadratic_sex=list(formula=~Quadratic+sex)
  p.Quadratic_age=list(formula=~Quadratic+age)
  p.Quadratic_age_sex=list(formula=~Quadratic+age+sex)
  p.Time_age_sex=list(formula=~Time+age+sex)
  p.time_age_sex=list(formula=~time+age+sex)
  cml=create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
  model.list=mark.wrapper(cml,data=Andava_KiP.proc,ddl=Andava_KiP.ddl,
                          threads = 8, output = FALSE) #runs model list in MARK
  export.MARK(Andava_KiP.proc, project.name="KiP_juv_adult_full_analysis",model.list,replace=T)
  return(model.list) #stores completed model list
}

MP_juv_adult_full_analysis=function() 
{
  setwd("/home/luke/comparative_ASR/Juv_Adult_survival_analysis/RMark_output_files/MP/") 
  Phi.age.time=list(formula=~age+time)
  Phi.sex.time=list(formula=~sex+time)
  Phi.sex.age.time=list(formula=~sex+age+time)
  Phi.sexxtime=list(formula=~sex*time)
  Phi.agextime=list(formula=~age*time)
  Phi.agexsexxtime=list(formula=~age*sex*time)
  Phi.age.Quadratic=list(formula=~age+Quadratic)
  Phi.sex.Quadratic=list(formula=~sex+Quadratic)
  Phi.sexxQuadratic=list(formula=~sex*Quadratic)
  Phi.agexQuadratic=list(formula=~age*Quadratic)
  Phi.sex=list(formula=~sex)
  Phi.age=list(formula=~age)
  Phi.age.sex=list(formula=~age+sex) 
  Phi.agexsex=list(formula=~age*sex) 
  Phi.dot=list(formula=~1) 
  Phi.time=list(formula=~time)
  Phi.Time=list(formula=~Time)
  Phi.Quadratic=list(formula=~Quadratic)
  Phi.age.Time=list(formula=~age+Time)
  Phi.sex.Time=list(formula=~sex+Time)
  Phi.sexxTime=list(formula=~sex*Time)
  Phi.agexTime=list(formula=~age*Time)
  Phi.sex.age.Time=list(formula=~sex+age+Time)
  Phi.sexXageXTime=list(formula=~sex*age*Time)
  p.sex=list(formula=~sex)
  p.age=list(formula=~age) 
  p.dot=list(formula=~1)
  p.time=list(formula=~time)
  p.Time=list(formula=~Time)
  p.Quadratic=list(formula=~Quadratic)
  p.sexxtime=list(formula=~sex*time)
  p.agextime=list(formula=~age*time)
  p.sexxTime=list(formula=~sex*Time)
  p.agexTime=list(formula=~age*Time)
  p.agexsex=list(formula=~age*sex)
  p.Quadraticxsex=list(formula=~Quadratic*sex)
  p.Quadraticxage=list(formula=~Quadratic*age)
  p.Quadraticxagexsex=list(formula=~Quadratic*age*sex)
  p.Timexagexsex=list(formula=~Time*age*sex)
  p.timexagexsex=list(formula=~time*age*sex)
  p.sex_time=list(formula=~sex+time)
  p.age_time=list(formula=~age+time)
  p.sex_Time=list(formula=~sex+Time)
  p.age_Time=list(formula=~age+Time)
  p.age_sex=list(formula=~age+sex)
  p.Quadratic_sex=list(formula=~Quadratic+sex)
  p.Quadratic_age=list(formula=~Quadratic+age)
  p.Quadratic_age_sex=list(formula=~Quadratic+age+sex)
  p.Time_age_sex=list(formula=~Time+age+sex)
  p.time_age_sex=list(formula=~time+age+sex)
  cml=create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
  model.list=mark.wrapper(cml,data=Andava_MP.proc,ddl=Andava_MP.ddl,
                          threads = 8, output = FALSE) #runs model list in MARK
  export.MARK(Andava_MP.proc, project.name="MP_juv_adult_full_analysis",model.list,replace=T)
  return(model.list) #stores completed model list
}

WfP_juv_adult_full_analysis=function() 
{
  setwd("/home/luke/comparative_ASR/Juv_Adult_survival_analysis/RMark_output_files/WfP") 
  Phi.age.time=list(formula=~age+time)
  Phi.sex.time=list(formula=~sex+time)
  Phi.sex.age.time=list(formula=~sex+age+time)
  Phi.sexxtime=list(formula=~sex*time)
  Phi.agextime=list(formula=~age*time)
  Phi.agexsexxtime=list(formula=~age*sex*time)
  Phi.age.Quadratic=list(formula=~age+Quadratic)
  Phi.sex.Quadratic=list(formula=~sex+Quadratic)
  Phi.sexxQuadratic=list(formula=~sex*Quadratic)
  Phi.agexQuadratic=list(formula=~age*Quadratic)
  Phi.sex=list(formula=~sex)
  Phi.age=list(formula=~age)
  Phi.age.sex=list(formula=~age+sex) 
  Phi.agexsex=list(formula=~age*sex) 
  Phi.dot=list(formula=~1) 
  Phi.time=list(formula=~time)
  Phi.Time=list(formula=~Time)
  Phi.Quadratic=list(formula=~Quadratic)
  Phi.age.Time=list(formula=~age+Time)
  Phi.sex.Time=list(formula=~sex+Time)
  Phi.sexxTime=list(formula=~sex*Time)
  Phi.agexTime=list(formula=~age*Time)
  Phi.sex.age.Time=list(formula=~sex+age+Time)
  Phi.sexXageXTime=list(formula=~sex*age*Time)
  p.sex=list(formula=~sex)
  p.age=list(formula=~age) 
  p.dot=list(formula=~1)
  p.time=list(formula=~time)
  p.Time=list(formula=~Time)
  p.Quadratic=list(formula=~Quadratic)
  p.sexxtime=list(formula=~sex*time)
  p.agextime=list(formula=~age*time)
  p.sexxTime=list(formula=~sex*Time)
  p.agexTime=list(formula=~age*Time)
  p.agexsex=list(formula=~age*sex)
  p.Quadraticxsex=list(formula=~Quadratic*sex)
  p.Quadraticxage=list(formula=~Quadratic*age)
  p.Quadraticxagexsex=list(formula=~Quadratic*age*sex)
  p.Timexagexsex=list(formula=~Time*age*sex)
  p.timexagexsex=list(formula=~time*age*sex)
  p.sex_time=list(formula=~sex+time)
  p.age_time=list(formula=~age+time)
  p.sex_Time=list(formula=~sex+Time)
  p.age_Time=list(formula=~age+Time)
  p.age_sex=list(formula=~age+sex)
  p.Quadratic_sex=list(formula=~Quadratic+sex)
  p.Quadratic_age=list(formula=~Quadratic+age)
  p.Quadratic_age_sex=list(formula=~Quadratic+age+sex)
  p.Time_age_sex=list(formula=~Time+age+sex)
  p.time_age_sex=list(formula=~time+age+sex)
  cml=create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
  model.list=mark.wrapper(cml,data=Andava_WfP.proc,ddl=Andava_WfP.ddl,
                          threads = 8, output = FALSE) #runs model list in MARK
  export.MARK(Andava_WfP.proc, project.name="WfP_juv_adult_full_analysis",model.list,replace=T)
  return(model.list) #stores completed model list
}

Maio_juv_adult_full_analysis=function() 
{
  setwd("/home/luke/comparative_ASR/Juv_Adult_survival_analysis/RMark_output_files/Maio/") 
  Phi.age.time=list(formula=~age+time)
  Phi.sex.time=list(formula=~sex+time)
  Phi.sex.age.time=list(formula=~sex+age+time)
  Phi.sexxtime=list(formula=~sex*time)
  Phi.agextime=list(formula=~age*time)
  Phi.agexsexxtime=list(formula=~age*sex*time)
  Phi.age.Quadratic=list(formula=~age+Quadratic)
  Phi.sex.Quadratic=list(formula=~sex+Quadratic)
  Phi.sexxQuadratic=list(formula=~sex*Quadratic)
  Phi.agexQuadratic=list(formula=~age*Quadratic)
  Phi.sex=list(formula=~sex)
  Phi.age=list(formula=~age)
  Phi.age.sex=list(formula=~age+sex) 
  Phi.agexsex=list(formula=~age*sex) 
  Phi.dot=list(formula=~1) 
  Phi.time=list(formula=~time)
  Phi.Time=list(formula=~Time)
  Phi.Quadratic=list(formula=~Quadratic)
  Phi.age.Time=list(formula=~age+Time)
  Phi.sex.Time=list(formula=~sex+Time)
  Phi.sexxTime=list(formula=~sex*Time)
  Phi.agexTime=list(formula=~age*Time)
  Phi.sex.age.Time=list(formula=~sex+age+Time)
  Phi.sexXageXTime=list(formula=~sex*age*Time)
  p.sex=list(formula=~sex)
  p.age=list(formula=~age) 
  p.dot=list(formula=~1)
  p.time=list(formula=~time)
  p.Time=list(formula=~Time)
  p.Quadratic=list(formula=~Quadratic)
  p.sexxtime=list(formula=~sex*time)
  p.agextime=list(formula=~age*time)
  p.sexxTime=list(formula=~sex*Time)
  p.agexTime=list(formula=~age*Time)
  p.agexsex=list(formula=~age*sex)
  p.Quadraticxsex=list(formula=~Quadratic*sex)
  p.Quadraticxage=list(formula=~Quadratic*age)
  p.Quadraticxagexsex=list(formula=~Quadratic*age*sex)
  p.Timexagexsex=list(formula=~Time*age*sex)
  p.timexagexsex=list(formula=~time*age*sex)
  p.sex_time=list(formula=~sex+time)
  p.age_time=list(formula=~age+time)
  p.sex_Time=list(formula=~sex+Time)
  p.age_Time=list(formula=~age+Time)
  p.age_sex=list(formula=~age+sex)
  p.Quadratic_sex=list(formula=~Quadratic+sex)
  p.Quadratic_age=list(formula=~Quadratic+age)
  p.Quadratic_age_sex=list(formula=~Quadratic+age+sex)
  p.Time_age_sex=list(formula=~Time+age+sex)
  p.time_age_sex=list(formula=~time+age+sex)
  cml=create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
  model.list=mark.wrapper(cml,data=Maio_kp.proc,ddl=Maio_kp.ddl,
                          threads = 8, output = FALSE) #runs model list in MARK
  export.MARK(Maio_kp.proc, project.name="Maio_juv_adult_full_analysis",model.list,replace=T)
  return(model.list) #stores completed model list
}

Tuzla_juv_adult_full_analysis=function() 
{
  setwd("/home/luke/comparative_ASR/Juv_Adult_survival_analysis/RMark_output_files/Tuzla/") 
  Phi.age.time=list(formula=~age+time) #Phi(age(t))
  Phi.sex.time=list(formula=~sex+time) #Phi(sex(t))
  Phi.sex.age.time=list(formula=~sex+age+time) #Phi(sex(t))
  Phi.sexxtime=list(formula=~sex*time)
  Phi.agextime=list(formula=~age*time)
  Phi.agexsexxtime=list(formula=~age*sex*time)
  Phi.age.Quadratic=list(formula=~age+Quadratic) #Phi(age(t))
  Phi.sex.Quadratic=list(formula=~sex+Quadratic) #Phi(sex(t))
  Phi.sexxQuadratic=list(formula=~sex*Quadratic)
  Phi.agexQuadratic=list(formula=~age*Quadratic)
  Phi.sex=list(formula=~sex) #Phi(sex(.))
  Phi.age=list(formula=~age) #Phi(age(.))
  Phi.age.sex=list(formula=~age+sex) 
  Phi.agexsex=list(formula=~age*sex) 
  Phi.dot=list(formula=~1) #p=Phi(.)
  Phi.time=list(formula=~time) #Phi(t)
  Phi.Time=list(formula=~Time)
  Phi.Quadratic=list(formula=~Quadratic)
  Phi.age.Time=list(formula=~age+Time)
  Phi.sex.Time=list(formula=~sex+Time)
  Phi.sexxTime=list(formula=~sex*Time)
  Phi.agexTime=list(formula=~age*Time)
  Phi.sex.age.Time=list(formula=~sex+age+Time) #Phi(sex(t))
  Phi.sexXageXTime=list(formula=~sex*age*Time) #Phi(sex(t))
  p.sex=list(formula=~sex) #p(sex(.))
  p.age=list(formula=~age) #p(age(.))
  p.dot=list(formula=~1) #p(.)
  p.time=list(formula=~time) #p(t)
  p.Time=list(formula=~Time)
  p.Quadratic=list(formula=~Quadratic)
  p.sexxtime=list(formula=~sex*time)
  p.agextime=list(formula=~age*time)
  p.sexxTime=list(formula=~sex*Time)
  p.agexTime=list(formula=~age*Time)
  p.agexsex=list(formula=~age*sex)
  p.Quadraticxsex=list(formula=~Quadratic*sex)
  p.Quadraticxage=list(formula=~Quadratic*age)
  p.Quadraticxagexsex=list(formula=~Quadratic*age*sex)
  p.Timexagexsex=list(formula=~Time*age*sex)
  p.timexagexsex=list(formula=~time*age*sex)
  p.sex_time=list(formula=~sex+time)
  p.age_time=list(formula=~age+time)
  p.sex_Time=list(formula=~sex+Time)
  p.age_Time=list(formula=~age+Time)
  p.age_sex=list(formula=~age+sex)
  p.Quadratic_sex=list(formula=~Quadratic+sex)
  p.Quadratic_age=list(formula=~Quadratic+age)
  p.Quadratic_age_sex=list(formula=~Quadratic+age+sex)
  p.Time_age_sex=list(formula=~Time+age+sex)
  p.time_age_sex=list(formula=~time+age+sex)
  cml=create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
  model.list=mark.wrapper(cml,data=Tuzla_kp.proc,ddl=Tuzla_kp.ddl,
                          threads = 8, output = FALSE) #runs model list in MARK
  export.MARK(Tuzla_kp.proc, project.name="Tuzla_juv_adult_full_analysis",model.list,replace=T)
  return(model.list) #stores completed model list
}

# Run models
KiP_juv_adult_full_analysis_run <- KiP_juv_adult_full_analysis()
WfP_juv_adult_full_analysis_run <- WfP_juv_adult_full_analysis()
MP_juv_adult_full_analysis_run <- MP_juv_adult_full_analysis()
Maio_juv_adult_full_analysis_run <- Maio_juv_adult_full_analysis()
Tuzla_juv_adult_full_analysis_run <- Tuzla_juv_adult_full_analysis()

# Inspect model output table ranked by AIC
KiP_juv_adult_full_analysis_run
WfP_juv_adult_full_analysis_run
MP_juv_adult_full_analysis_run
Maio_juv_adult_full_analysis_run
Tuzla_juv_adult_full_analysis_run

write.table(KiP_juv_adult_full_analysis_run$model.table,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/KiP_juv_adult_full_analysis_run_model_table.txt", sep = "\t")
write.table(WfP_juv_adult_full_analysis_run$model.table,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/WfP_juv_adult_full_analysis_run_model_table.txt", sep = "\t")
write.table(MP_juv_adult_full_analysis_run$model.table,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/MP_juv_adult_full_analysis_run_model_table.txt", sep = "\t")
write.table(Maio_juv_adult_full_analysis_run$model.table,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/Maio_juv_adult_full_analysis_run_model_table.txt", sep = "\t")
write.table(Tuzla_juv_adult_full_analysis_run$model.table,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/Tuzla_juv_adult_full_analysis_run_model_table.txt", sep = "\t")

# find models only with the encounter structure "~1" and subset the
# AIC model list. This will reveal the best structure for p.
# KiP, best p structure is 1) p(~Year * Time) or 2) p(~Year * Quadratic)
KiP_not_constant_phi_models <- rownames(KiP_juv_adult_full_analysis_run$model.table[which(KiP_juv_adult_full_analysis_run$model.table$Phi != "~1"), ])
KiP_juv_adult_constant_phi_mod_list <- remove.mark(KiP_juv_adult_full_analysis_run, as.numeric(KiP_not_constant_phi_models))

# find models only with the encounter structure "~age * time" and subset the
# AIC model list. This wil reveal the best structure for phi, while using the
# best structure for p.
# KiP, best for Phi is Phi(~Sex * Time), therefore model 637 is best
KiP_not_age_time_models <- rownames(KiP_juv_adult_full_analysis_run$model.table[which(KiP_juv_adult_full_analysis_run$model.table$p != "~age * time"), ])
KiP_juv_adult_p_age_time_mod_list <- remove.mark(KiP_juv_adult_full_analysis_run, as.numeric(KiP_not_age_time_models))
rownames(KiP_juv_adult_full_analysis_run$model.table[which(KiP_juv_adult_full_analysis_run$model.table$Phi == "~age * sex" & 
                                                          KiP_juv_adult_full_analysis_run$model.table$p == "~age * time"), ])

write.table(KiP_juv_adult_p_age_time_mod_list$model.table,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/KiP_juv_adult_p_age_time_mod_list_model_table.txt", sep = "\t")

# find models only with the encounter structure "~1" and subset the
# AIC model list. This will reveal the best structure for p.
# WfP, best p structure is 1) p(~Year * Time) or 2) p(~Year * Quadratic)
WfP_not_constant_phi_models <- rownames(WfP_juv_adult_full_analysis_run$model.table[which(WfP_juv_adult_full_analysis_run$model.table$Phi != "~1"), ])
WfP_juv_adult_constant_phi_mod_list <- remove.mark(WfP_juv_adult_full_analysis_run, as.numeric(WfP_not_constant_phi_models))

# find models only with the encounter structure "~age * time" and subset the
# AIC model list. This wil reveal the best structure for phi, while using the
# best structure for p.
# WfP, best for Phi is Phi(~Sex * Time), therefore model 637 is best
WfP_not_time_age_sex_models <- rownames(WfP_juv_adult_full_analysis_run$model.table[which(WfP_juv_adult_full_analysis_run$model.table$p != "~time + age + sex"), ])
WfP_juv_adult_p_time_age_sex_mod_list <- remove.mark(WfP_juv_adult_full_analysis_run, as.numeric(WfP_not_time_age_sex_models))
rownames(WfP_juv_adult_full_analysis_run$model.table[which(WfP_juv_adult_full_analysis_run$model.table$Phi == "~age * sex" & 
                                                             WfP_juv_adult_full_analysis_run$model.table$p == "~time + age + sex"), ])

write.table(WfP_juv_adult_p_time_age_sex_mod_list$model.table,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/WfP_juv_adult_p_time_age_sex_mod_list_model_table.txt", sep = "\t")

# find models only with the encounter structure "~1" and subset the
# AIC model list. This will reveal the best structure for p.
# MP, best p structure is 1) p(~Year * Time) or 2) p(~Year * Quadratic)
MP_not_constant_phi_models <- rownames(MP_juv_adult_full_analysis_run$model.table[which(MP_juv_adult_full_analysis_run$model.table$Phi != "~1"), ])
MP_juv_adult_constant_phi_mod_list <- remove.mark(MP_juv_adult_full_analysis_run, as.numeric(MP_not_constant_phi_models))

# find models only with the encounter structure "~age * time" and subset the
# AIC model list. This wil reveal the best structure for phi, while using the
# best structure for p.
# MP, best for Phi is Phi(~Sex * Time), therefore model 637 is best
MP_not_age_time_models <- rownames(MP_juv_adult_full_analysis_run$model.table[which(MP_juv_adult_full_analysis_run$model.table$p != "~age + time"), ])
MP_juv_adult_p_age_time_mod_list <- remove.mark(MP_juv_adult_full_analysis_run, as.numeric(MP_not_age_time_models))
rownames(MP_juv_adult_full_analysis_run$model.table[which(MP_juv_adult_full_analysis_run$model.table$Phi == "~age * sex" & 
                                                             MP_juv_adult_full_analysis_run$model.table$p == "~age + time"), ])

write.table(MP_juv_adult_p_age_time_mod_list$model.table,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/MP_juv_adult_p_age_time_mod_list_model_table.txt", sep = "\t")

# find models only with the encounter structure "~1" and subset the
# AIC model list. This will reveal the best structure for p.
# Maio, best p structure is 1) p(~Year * Time) or 2) p(~Year * Quadratic)
Maio_not_constant_phi_models <- rownames(Maio_juv_adult_full_analysis_run$model.table[which(Maio_juv_adult_full_analysis_run$model.table$Phi != "~1"), ])
Maio_juv_adult_constant_phi_mod_list <- remove.mark(Maio_juv_adult_full_analysis_run, as.numeric(Maio_not_constant_phi_models))

# find models only with the encounter structure "~age * time" and subset the
# AIC model list. This wil reveal the best structure for phi, while using the
# best structure for p.
# Maio, best for Phi is Phi(~Sex * Time), therefore model 637 is best
Maio_not_age_x_time_models <- rownames(Maio_juv_adult_full_analysis_run$model.table[which(Maio_juv_adult_full_analysis_run$model.table$p != "~age * time"), ])
Maio_juv_adult_p_age_x_time_mod_list <- remove.mark(Maio_juv_adult_full_analysis_run, as.numeric(Maio_not_age_x_time_models))
rownames(Maio_juv_adult_full_analysis_run$model.table[which(Maio_juv_adult_full_analysis_run$model.table$Phi == "~age * sex" & 
                                                            Maio_juv_adult_full_analysis_run$model.table$p == "~age * time"), ])

write.table(Maio_juv_adult_p_age_x_time_mod_list$model.table,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/Maio_juv_adult_p_age_x_time_mod_list_model_table.txt", sep = "\t")

# find models only with the encounter structure "~1" and subset the
# AIC model list. This will reveal the best structure for p.
# Tuzla, best p structure is 1) p(~Year * Time) or 2) p(~Year * Quadratic)
Tuzla_not_constant_phi_models <- rownames(Tuzla_juv_adult_full_analysis_run$model.table[which(Tuzla_juv_adult_full_analysis_run$model.table$Phi != "~1"), ])
Tuzla_juv_adult_constant_phi_mod_list <- remove.mark(Tuzla_juv_adult_full_analysis_run, as.numeric(Tuzla_not_constant_phi_models))

# find models only with the encounter structure "~age * time" and subset the
# AIC model list. This wil reveal the best structure for phi, while using the
# best structure for p.
# Tuzla, best for Phi is Phi(~Sex * Time), therefore model 637 is best
Tuzla_not_time_x_age_x_sex_models <- rownames(Tuzla_juv_adult_full_analysis_run$model.table[which(Tuzla_juv_adult_full_analysis_run$model.table$p != "~time * age * sex"), ])
Tuzla_juv_adult_p_time_x_age_x_sex_mod_list <- remove.mark(Tuzla_juv_adult_full_analysis_run, as.numeric(Tuzla_not_time_x_age_x_sex_models))
rownames(Tuzla_juv_adult_full_analysis_run$model.table[which(Tuzla_juv_adult_full_analysis_run$model.table$Phi == "~age * sex" & 
                                                              Tuzla_juv_adult_full_analysis_run$model.table$p == "~time * age * sex"), ])

write.table(Tuzla_juv_adult_p_time_x_age_x_sex_mod_list$model.table,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/Tuzla_juv_adult_p_time_x_age_x_sex_mod_list_model_table.txt", sep = "\t")

# find models only with the encounter structure "~1" and subset the
# AIC model list. This will reveal the best structure for p.
not_constant_phi_models <- rownames(Ceuta_juv_adult_full_analysis_Run$model.table[which(Ceuta_juv_adult_full_analysis_Run$model.table$Phi != "~1"), ])
Ceuta_juv_adult_constant_phi_mod_list <- remove.mark(Ceuta_juv_adult_full_analysis_Run, as.numeric(not_constant_phi_models))

# find models only with the encounter structure "~Year * Quadratic" and subset the
# AIC model list. This wil reveal the best structure for phi, while using the
# best structure for p.
not_age_x_time_models <- rownames(Ceuta_juv_adult_full_analysis_Run$model.table[which(Ceuta_juv_adult_full_analysis_Run$model.table$p != "~age * time"), ])
Ceuta_juv_adult_p_age_x_time_mod_list <- remove.mark(Ceuta_juv_adult_full_analysis_Run, as.numeric(not_age_x_time_models))
rownames(Ceuta_juv_adult_full_analysis_Run$model.table[which(Ceuta_juv_adult_full_analysis_Run$model.table$Phi == "~age * sex" & 
                                                               Ceuta_juv_adult_full_analysis_Run$model.table$p == "~age * time"), ])

write.table(Ceuta_juv_adult_p_age_x_time_mod_list$model.table,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/Ceuta_juv_adult_p_age_x_time_mod_list_model_table.txt", sep = "\t")

###############################################################################
# KiP plotting and vital rate estimation
# Plot of Phi(~age * sex)p(~age * time)
# -------------------------------------------------------------
KiP_Sex_age_reals <- KiP_juv_adult_full_analysis_run[[162]]$results$real
Groups <- data.frame(str_split_fixed(rownames(KiP_Sex_age_reals), " ", n = 5))
KiP_Sex_age_reals <- cbind(Groups, KiP_Sex_age_reals)
KiP_Sex_age_reals <- KiP_Sex_age_reals[which(KiP_Sex_age_reals$X1 == "Phi"),]
KiP_Sex_age_reals$age <- unlist(str_extract_all(KiP_Sex_age_reals$X2,"[AJ]"))
KiP_Sex_age_reals$age <- as.factor(ifelse(KiP_Sex_age_reals$age == "A","Adult","Juvenile"))
KiP_Sex_age_reals$Sex <- unlist(str_extract_all(KiP_Sex_age_reals$X2,"[FM]"))
KiP_Sex_age_reals$Sex <- as.factor(ifelse(KiP_Sex_age_reals$Sex == "F","Female","Male"))
KiP_Sex_age_reals$Year <- as.factor(unlist(substr(KiP_Sex_age_reals$X5, 2, 5)))
KiP_Sex_age_reals$KiP_Sex_age <- paste(KiP_Sex_age_reals$Sex,KiP_Sex_age_reals$age,sep = "_")
limits <- aes(ymin = lcl, ymax = ucl)
dodge <- position_dodge(width=0.7)

KiP_Sex_age_reals$species <- "Kittlitz's"

KiP_Sex_age_survival <- 
  rbind(KiP_Sex_age_reals[,c("Sex", "age", "estimate", "lcl", "ucl", "species")],
        KiP_Sex_chick_survival)
KiP_Sex_age_survival$age <- 
  factor(KiP_Sex_age_survival$age,
         levels = c("Chick", "Juvenile", "Adult"))
KiP_Sex_age_survival$Sex <- 
  factor(KiP_Sex_age_survival$Sex,
         levels = c("Female", "Male"))
cbPalette <- c("#C0504D", "#558ED5")
KiP_Sex_age_plot2 <- ggplot(KiP_Sex_age_survival, aes(fill = Sex, x = age, y = estimate)) + 
  theme_bw() +
  geom_bar(width = 0.7, position="dodge", stat="identity") +
  geom_errorbar(limits, position=dodge, width=0.2, size = .6, shape = 1) +
  theme(text=element_text(size=16, family="Candara"),
        legend.position = "none",
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        axis.title.x = element_text(size=33, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=33, vjust=1.2),
        axis.text.y  = element_text(size=30),  
        plot.title = element_text(vjust=1.2, size=25, face="bold"),
        panel.grid.major = element_blank()) +
  xlab("Life stage") + 
  ylab("Apparent survival (± 95% CI)") +
  #ggtitle("Kittlitz's plover stage- and sex-specifc survival\n(2012 to 2014)") +
  scale_y_continuous(limits = c(0, 1)) +
  #scale_fill_hue(l=55)
  scale_fill_manual(values=cbPalette)
KiP_Sex_age_plot2

ggsave(filename = "Andava_KiP_SSD_ASR.tiff", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/RMARK_model_output/SSD_Plots",
       width = 8,
       height = 8, units = "in",
       dpi = 300)

###############################################################################
# WfP plotting and vital rate estimation
# Plot of Phi(~age * sex)p(~time + age + sex)
# -------------------------------------------------------------
WfP_Sex_age_reals <- WfP_juv_adult_full_analysis_run[[179]]$results$real
Groups <- data.frame(str_split_fixed(rownames(WfP_Sex_age_reals), " ", n = 5))
WfP_Sex_age_reals <- cbind(Groups, WfP_Sex_age_reals)
WfP_Sex_age_reals <- WfP_Sex_age_reals[which(WfP_Sex_age_reals$X1 == "Phi"),]
WfP_Sex_age_reals$age <- unlist(str_extract_all(WfP_Sex_age_reals$X2,"[AJ]"))
WfP_Sex_age_reals$age <- as.factor(ifelse(WfP_Sex_age_reals$age == "A","Adult","Juvenile"))
WfP_Sex_age_reals$Sex <- unlist(str_extract_all(WfP_Sex_age_reals$X2,"[FM]"))
WfP_Sex_age_reals$Sex <- as.factor(ifelse(WfP_Sex_age_reals$Sex == "F","Female","Male"))
WfP_Sex_age_reals$Year <- as.factor(unlist(substr(WfP_Sex_age_reals$X5, 2, 5)))
WfP_Sex_age_reals$WfP_Sex_age <- paste(WfP_Sex_age_reals$Sex,WfP_Sex_age_reals$age,sep = "_")
limits <- aes(ymin = lcl, ymax = ucl)
dodge <- position_dodge(width=0.7)

WfP_Sex_age_reals$species <- "White-fronted"

WfP_Sex_age_survival <- 
  rbind(WfP_Sex_age_reals[,c("Sex", "age", "estimate", "lcl", "ucl", "species")],
        WfP_Sex_chick_survival)
WfP_Sex_age_survival$age <- 
  factor(WfP_Sex_age_survival$age,
         levels = c("Chick", "Juvenile", "Adult"))
WfP_Sex_age_survival$Sex <- 
  factor(WfP_Sex_age_survival$Sex,
         levels = c("Female", "Male"))
cbPalette <- c("#C0504D", "#558ED5")
WfP_Sex_age_plot2 <- ggplot(WfP_Sex_age_survival, aes(fill = Sex, x = age, y = estimate)) + 
  theme_bw() +
  geom_bar(width = 0.7, position="dodge", stat="identity") +
  geom_errorbar(limits, position=dodge, width=0.2, size = .6, shape = 1) +
  theme(text=element_text(size=16, family="Candara"),
        legend.position = "none",
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        axis.title.x = element_text(size=33, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=33, vjust=1.2),
        axis.text.y  = element_text(size=30),  
        plot.title = element_text(vjust=1.2, size=25, face="bold"),
        panel.grid.major = element_blank()) +
  xlab("Life stage") + 
  ylab("Apparent survival (± 95% CI)") +
  #ggtitle("Kittlitz's plover stage- and sex-specifc survival\n(2012 to 2014)") +
  scale_y_continuous(limits = c(0, 1)) +
  #scale_fill_hue(l=55)
  scale_fill_manual(values=cbPalette)
WfP_Sex_age_plot2

###############################################################################
# MP plotting and vital rate estimation
# Plot of Phi(~age * sex)p(~age + time)
# -------------------------------------------------------------
MP_Sex_age_reals <- MP_juv_adult_full_analysis_run[[159]]$results$real
Groups <- data.frame(str_split_fixed(rownames(MP_Sex_age_reals), " ", n = 5))
MP_Sex_age_reals <- cbind(Groups, MP_Sex_age_reals)
MP_Sex_age_reals <- MP_Sex_age_reals[which(MP_Sex_age_reals$X1 == "Phi"),]
MP_Sex_age_reals$age <- unlist(str_extract_all(MP_Sex_age_reals$X2,"[AJ]"))
MP_Sex_age_reals$age <- as.factor(ifelse(MP_Sex_age_reals$age == "A","Adult","Juvenile"))
MP_Sex_age_reals$Sex <- unlist(str_extract_all(MP_Sex_age_reals$X2,"[FM]"))
MP_Sex_age_reals$Sex <- as.factor(ifelse(MP_Sex_age_reals$Sex == "F","Female","Male"))
MP_Sex_age_reals$Year <- as.factor(unlist(substr(MP_Sex_age_reals$X5, 2, 5)))
MP_Sex_age_reals$MP_Sex_age <- paste(MP_Sex_age_reals$Sex,MP_Sex_age_reals$age,sep = "_")
limits <- aes(ymin = lcl, ymax = ucl)
dodge <- position_dodge(width=0.7)

MP_Sex_age_reals$species <- "Madagascar"

MP_Sex_age_survival <- 
  rbind(MP_Sex_age_reals[,c("Sex", "age", "estimate", "lcl", "ucl", "species")],
        MP_Sex_chick_survival)
MP_Sex_age_survival$age <- 
  factor(MP_Sex_age_survival$age,
         levels = c("Chick", "Juvenile", "Adult"))
MP_Sex_age_survival$Sex <- 
  factor(MP_Sex_age_survival$Sex,
         levels = c("Female", "Male"))
cbPalette <- c("#C0504D", "#558ED5")
MP_Sex_age_plot2 <- ggplot(MP_Sex_age_survival, aes(fill = Sex, x = age, y = estimate)) + 
  theme_bw() +
  geom_bar(width = 0.7, position="dodge", stat="identity") +
  geom_errorbar(limits, position=dodge, width=0.2, size = .6, shape = 1) +
  theme(text=element_text(size=16, family="Candara"),
        legend.position = "none",
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        axis.title.x = element_text(size=33, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=33, vjust=1.2),
        axis.text.y  = element_text(size=30),  
        plot.title = element_text(vjust=1.2, size=25, face="bold"),
        panel.grid.major = element_blank()) +
  xlab("Life stage") + 
  ylab("Apparent survival (± 95% CI)") +
  #ggtitle("Kittlitz's plover stage- and sex-specifc survival\n(2012 to 2014)") +
  scale_y_continuous(limits = c(0, 1)) +
  #scale_fill_hue(l=55)
  scale_fill_manual(values=cbPalette)
MP_Sex_age_plot2

###############################################################################
# Maio plotting and vital rate estimation
# Plot of Phi(~age * sex)p(~age * time)
# -------------------------------------------------------------
Maio_Sex_age_reals <- Maio_juv_adult_full_analysis_run[[162]]$results$real
Groups <- data.frame(str_split_fixed(rownames(Maio_Sex_age_reals), " ", n = 5))
Maio_Sex_age_reals <- cbind(Groups, Maio_Sex_age_reals)
Maio_Sex_age_reals <- Maio_Sex_age_reals[which(Maio_Sex_age_reals$X1 == "Phi"),]
Maio_Sex_age_reals$age <- unlist(str_extract_all(Maio_Sex_age_reals$X2,"[AJ]"))
Maio_Sex_age_reals$age <- as.factor(ifelse(Maio_Sex_age_reals$age == "A","Adult","Juvenile"))
Maio_Sex_age_reals$Sex <- unlist(str_extract_all(Maio_Sex_age_reals$X2,"[FM]"))
Maio_Sex_age_reals$Sex <- as.factor(ifelse(Maio_Sex_age_reals$Sex == "F","Female","Male"))
Maio_Sex_age_reals$Year <- as.factor(unlist(substr(Maio_Sex_age_reals$X5, 2, 5)))
Maio_Sex_age_reals$Maio_Sex_age <- paste(Maio_Sex_age_reals$Sex,Maio_Sex_age_reals$age,sep = "_")
limits <- aes(ymin = lcl, ymax = ucl)
dodge <- position_dodge(width=0.7)

Maio_Sex_age_reals$species <- "Kentish (Maio)"

Maio_Sex_age_survival <- 
  rbind(Maio_Sex_age_reals[,c("Sex", "age", "estimate", "lcl", "ucl", "species")],
        Maio_Sex_chick_survival)
Maio_Sex_age_survival$age <- 
  factor(Maio_Sex_age_survival$age,
         levels = c("Chick", "Juvenile", "Adult"))
Maio_Sex_age_survival$Sex <- 
  factor(Maio_Sex_age_survival$Sex,
         levels = c("Female", "Male"))
cbPalette <- c("#C0504D", "#558ED5")
Maio_Sex_age_plot2 <- ggplot(Maio_Sex_age_survival, aes(fill = Sex, x = age, y = estimate)) + 
  theme_bw() +
  geom_bar(width = 0.7, position="dodge", stat="identity") +
  geom_errorbar(limits, position=dodge, width=0.2, size = .6, shape = 1) +
  theme(text=element_text(size=16, family="Candara"),
        legend.position = "none",
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        axis.title.x = element_text(size=33, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=33, vjust=1.2),
        axis.text.y  = element_text(size=30),  
        plot.title = element_text(vjust=1.2, size=25, face="bold"),
        panel.grid.major = element_blank()) +
  xlab("Life stage") + 
  ylab("Apparent survival (± 95% CI)") +
  #ggtitle("Kittlitz's plover stage- and sex-specifc survival\n(2012 to 2014)") +
  scale_y_continuous(limits = c(0, 1)) +
  #scale_fill_hue(l=55)
  scale_fill_manual(values=cbPalette)
Maio_Sex_age_plot2

###############################################################################
# Tuzla plotting and vital rate estimation
# Plot of Phi(~age * sex)p(~time * age * sex)
# -------------------------------------------------------------
Tuzla_Sex_age_reals <- Tuzla_juv_adult_full_analysis_run[[181]]$results$real
Groups <- data.frame(str_split_fixed(rownames(Tuzla_Sex_age_reals), " ", n = 5))
Tuzla_Sex_age_reals <- cbind(Groups, Tuzla_Sex_age_reals)
Tuzla_Sex_age_reals <- Tuzla_Sex_age_reals[which(Tuzla_Sex_age_reals$X1 == "Phi"),]
Tuzla_Sex_age_reals$age <- unlist(str_extract_all(Tuzla_Sex_age_reals$X2,"[AJ]"))
Tuzla_Sex_age_reals$age <- as.factor(ifelse(Tuzla_Sex_age_reals$age == "A","Adult","Juvenile"))
Tuzla_Sex_age_reals$Sex <- unlist(str_extract_all(Tuzla_Sex_age_reals$X2,"[FM]"))
Tuzla_Sex_age_reals$Sex <- as.factor(ifelse(Tuzla_Sex_age_reals$Sex == "F","Female","Male"))
Tuzla_Sex_age_reals$Year <- as.factor(unlist(substr(Tuzla_Sex_age_reals$X5, 2, 5)))
Tuzla_Sex_age_reals$Tuzla_Sex_age <- paste(Tuzla_Sex_age_reals$Sex,Tuzla_Sex_age_reals$age,sep = "_")
limits <- aes(ymin = lcl, ymax = ucl)
dodge <- position_dodge(width=0.7)

Tuzla_Sex_age_reals$species <- "Kentish (Tuzla)"

Tuzla_Sex_age_survival <- 
  rbind(Tuzla_Sex_age_reals[,c("Sex", "age", "estimate", "lcl", "ucl", "species")],
        Tuzla_Sex_chick_survival)
Tuzla_Sex_age_survival$age <- 
  factor(Tuzla_Sex_age_survival$age,
         levels = c("Chick", "Juvenile", "Adult"))
Tuzla_Sex_age_survival$Sex <- 
  factor(Tuzla_Sex_age_survival$Sex,
         levels = c("Female", "Male"))
cbPalette <- c("#C0504D", "#558ED5")
Tuzla_Sex_age_plot2 <- ggplot(Tuzla_Sex_age_survival, aes(fill = Sex, x = age, y = estimate)) + 
  theme_bw() +
  geom_bar(width = 0.7, position="dodge", stat="identity") +
  geom_errorbar(limits, position=dodge, width=0.2, size = .6, shape = 1) +
  theme(text=element_text(size=16, family="Candara"),
        legend.position = "none",
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        axis.title.x = element_text(size=33, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=33, vjust=1.2),
        axis.text.y  = element_text(size=30),  
        plot.title = element_text(vjust=1.2, size=25, face="bold"),
        panel.grid.major = element_blank()) +
  xlab("Life stage") + 
  ylab("Apparent survival (± 95% CI)") +
  #ggtitle("Kittlitz's plover stage- and sex-specifc survival\n(2012 to 2014)") +
  scale_y_continuous(limits = c(0, 1)) +
  #scale_fill_hue(l=55)
  scale_fill_manual(values=cbPalette)
Tuzla_Sex_age_plot2

###############################################################################
# Ceuta plotting and vital rate estimation
# Plot of Phi(~age * sex)p(~time * age * sex)
# -------------------------------------------------------------
Ceuta_Sex_age_reals <- Ceuta_juv_adult_full_analysis_Run[[162]]$results$real
Groups <- data.frame(str_split_fixed(rownames(Ceuta_Sex_age_reals), " ", n = 5))
Ceuta_Sex_age_reals <- cbind(Groups, Ceuta_Sex_age_reals)
Ceuta_Sex_age_reals <- Ceuta_Sex_age_reals[which(Ceuta_Sex_age_reals$X1 == "Phi"),]
Ceuta_Sex_age_reals$age <- unlist(str_extract_all(Ceuta_Sex_age_reals$X2,"[AJ]"))
Ceuta_Sex_age_reals$age <- as.factor(ifelse(Ceuta_Sex_age_reals$age == "A","Adult","Juvenile"))
Ceuta_Sex_age_reals$Sex <- unlist(str_extract_all(Ceuta_Sex_age_reals$X2,"[FM]"))
Ceuta_Sex_age_reals$Sex <- as.factor(ifelse(Ceuta_Sex_age_reals$Sex == "F","Female","Male"))
Ceuta_Sex_age_reals$Year <- as.factor(unlist(substr(Ceuta_Sex_age_reals$X5, 2, 5)))
Ceuta_Sex_age_reals$Ceuta_Sex_age <- paste(Ceuta_Sex_age_reals$Sex,Ceuta_Sex_age_reals$age,sep = "_")
limits <- aes(ymin = lcl, ymax = ucl)
dodge <- position_dodge(width=0.7)

Ceuta_Sex_age_reals$species <- "Snowy"

Ceuta_Sex_age_survival <- 
  rbind(Ceuta_Sex_age_reals[,c("Sex", "age", "estimate", "lcl", "ucl", "species")],
        Ceuta_Sex_chick_survival)
Ceuta_Sex_age_survival$age <- 
  factor(Ceuta_Sex_age_survival$age,
         levels = c("Chick", "Juvenile", "Adult"))
Ceuta_Sex_age_survival$Sex <- 
  factor(Ceuta_Sex_age_survival$Sex,
         levels = c("Female", "Male"))
cbPalette <- c("#C0504D", "#558ED5")
Ceuta_Sex_age_plot2 <- ggplot(Ceuta_Sex_age_survival, aes(fill = Sex, x = age, y = estimate)) + 
  theme_bw() +
  geom_bar(width = 0.7, position="dodge", stat="identity") +
  geom_errorbar(limits, position=dodge, width=0.2, size = .6, shape = 1) +
  theme(text=element_text(size=16, family="Candara"),
        legend.position = "none",
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        axis.title.x = element_text(size=33, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=33, vjust=1.2),
        axis.text.y  = element_text(size=30),  
        plot.title = element_text(vjust=1.2, size=25, face="bold"),
        panel.grid.major = element_blank()) +
  xlab("Life stage") + 
  ylab("Apparent survival (± 95% CI)") +
  #ggtitle("Kittlitz's plover stage- and sex-specifc survival\n(2012 to 2014)") +
  scale_y_continuous(limits = c(0, 1)) +
  #scale_fill_hue(l=55)
  scale_fill_manual(values=cbPalette)
Ceuta_Sex_age_plot2

# Join all vital rate estimates together
All_pops_sex_chick_juv_adult_survival <- rbind(KiP_Sex_age_survival,
                                               WfP_Sex_age_survival,
                                               MP_Sex_age_survival,
                                               Maio_Sex_age_survival,
                                               Tuzla_Sex_age_survival,
                                               Ceuta_Sex_age_survival)

rownames(All_pops_sex_chick_juv_adult_survival) <- NULL
write.table(All_pops_sex_chick_juv_adult_survival,
            file = "/home/luke/comparative_ASR/Juv_Adult_survival_analysis/Results/All_pops_sex_chick_juv_adult_survival.txt", sep = "\t")

All_pops_sex_age_plot <- 
  ggplot(All_pops_sex_chick_juv_adult_survival, aes(fill = Sex, x = age, y = estimate)) + 
  theme_bw() +
  geom_bar(width = 0.7, position="dodge", stat="identity") +
  geom_errorbar(limits, position=dodge, width=0.2, size = .6, shape = 1) +
  theme(text=element_text(size=16, family="Candara"),
        legend.position = "none",
        axis.title.x = element_text(size=11, vjust=-0.1),
        axis.text.x  = element_text(size=10), 
        axis.title.y = element_text(size=11, vjust=1.2),
        axis.text.y  = element_text(size=10),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Life stage") + 
  ylab("Apparent survival (± 95% CI)") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values=cbPalette) +
  facet_grid(species ~ .)

cleanup(ask=FALSE)