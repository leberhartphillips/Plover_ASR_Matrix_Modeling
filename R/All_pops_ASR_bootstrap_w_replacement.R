#-----------------------------------------------------------------------------#
#               Matrix Modelling of ASR in Ceuta Snowy Plovers                #
#                         Bootstrapping proceedure                            #
#                           (with replacement)                                #
#                 Luke Eberhart-Phillips and Martin Stoffel                   #
#-----------------------------------------------------------------------------#

# Reference RMark functions
library(RMark) 
library(stringr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
MarkPath <- "/usr/local/bin/mark"
MarkViewer<-"nano"

# Function that builds the two-sex Lefkovitch matrix using the vital rates
# specified in the VR_list object.  The chick_surv option notifies that
# the VR_list includes chick survival rates (TRUE) or not (FALSE).
plover_matrix <- 
  function(VR_list, chick_surv = TRUE)
  {
    if(chick_surv)
    {
      # Construct population projection matrix
      # Layout:
      #                  F_1st_yr     F_Adt      M_1st_yr     M_Adt
      #            -----------------------------------------------
      # F_1st_yr   |            0        NA             0        NA
      # F_Adt      |  F_CS * F_JS      F_AS             0         0
      # M_1st_yr   |            0        NA             0        NA
      # M_Adt      |            0         0   M_CS * M_JS      M_AS
      
      # Define plover life-stages of the Ceuta snowy plover matrix model
      stages <- c("F_1st_yr",  "F_Adt",  "M_1st_yr",  "M_Adt")
      result <- matrix(c(0, NA, 0, NA, 
                         (VR_list$F_Chk_survl*VR_list$F_Juv_survl),
                         VR_list$F_Adt_survl, 0, 0,
                         0, NA, 0, NA,
                         0, 0, (VR_list$M_Chk_survl*VR_list$M_Juv_survl),
                         VR_list$M_Adt_survl),
                       nrow = 4, byrow = TRUE,
                       dimnames = list(stages, stages))
    }
    else
    {
      # Construct population projection matrix
      # Layout:
      #                  F_1st_yr     F_Adt      M_1st_yr     M_Adt
      #            -----------------------------------------------
      # F_1st_yr   |            0        NA             0        NA
      # F_Adt      |         F_JS      F_AS             0         0
      # M_1st_yr   |            0        NA             0        NA
      # M_Adt      |            0         0          M_JS      M_AS
      
      # Define plover life-stages of the Ceuta snowy plover matrix model
      stages <- c("F_1st_yr",  "F_Adt",  "M_1st_yr",  "M_Adt")
      result <- matrix(c(0, NA, 0, NA, 
                         VR_list$F_Juv_survl,
                         VR_list$F_Adt_survl, 0, 0,
                         0, NA, 0, NA,
                         0, 0, VR_list$M_Juv_survl,
                         VR_list$M_Adt_survl),
                       nrow = 4, byrow = TRUE,
                       dimnames = list(stages, stages))
    }
    result
  }

# Deterministic projection function that calculates the ASR of the stable 
# stage distribution while incorporating frequency-dependent fecundity.
#
# Arguments in the function include:
# A ---> an two sex x by x projection matrix
# n ---> an x lengthed vector representing starting stage distribution (the
#        default is a vector with 10 individuals in each stage)
# h ---> the harem size for the species.  h > 1 is polgynous, h < 1 is
#        polyandrous, h = 1 is monogamous
# k ---> the clutch size for the species.
# PSR ---> the primary sex ratio (default is 0.5)
# iterations ---> the number of iterations to simulate (default is 20)
freq_dep_SSD_ASR <- 
  function (A, n = rep(10, nrow(A)), h = 1, k = 1, iterations = 30, PSR = 0.5) 
  {
    x <- length(n) # Number of stages in matrix
    t <- iterations # Number of time steps to simulate
    stage <- matrix(numeric(x * t), nrow = x) # an empty t by x matrix
    for (i in 1:t) { # for loop that goes through each of t time steps
      stage[,i] <- n # stage distribution at time t
      M2 <- stage[4, i] # number of male adults at time t
      F2 <- stage[2, i] # number of female adults at time t
      A[1,x/2]        <- (k*M2)/(M2+(F2/h))*PSR # F freq-dep fecundity of F
      A[(x/4)*3,x/2]  <- (k*M2)/(M2+(F2/h))*PSR # F freq-dep fecundity of M
      A[1,x]          <- (k*F2)/(M2+(F2/h))*PSR # M freq-dep fecundity of F
      A[(x/4)*3,x]    <- (k*F2)/(M2+(F2/h))*PSR # M freq-dep fecundity of M
      n <- A %*% n # define the new n (i.e., new stage distribution at time t)
    }
    rownames(stage) <- rownames(A) # define rownames of stage matrix
    colnames(stage) <- 0:(t - 1) # define colnames of stage matrix
    w <- stage[, t] # define stable stage as the last stage
    stable.stage <- w/sum(w)
    ASR <- stable.stage[x]/(stable.stage[x/2]+stable.stage[x])
    pop.proj <- list(ASR = ASR, # make a list of results
                     stable.stage = stable.stage, 
                     stage.vectors = stage,
                     SSD_M2 = stable.stage[4],
                     SSD_F2 = stable.stage[2])
    pop.proj # print the list as output to the function
  }

# Import prepared capture history file (three columns: "ch" ... concatenated 
# annual presence/absence, "sex"... male or female, or NA for unknown sex of 
# chicks, "age"... Adult or juvenile [NOTE: this describes the stage that the 
# individual was initially ringed])
setwd("/home/luke/comparative_ASR/Juv_Adult_survival_analysis/") 

Ceuta_F_A <- 
  read.table("Data_files/Ceuta_juv_adult_capture_history_MARK_w_rings.txt",
             header=T,colClasses=c("factor","character","factor","factor"))
KiP_F_A <- 
  read.table("Data_files/Andava_KiP_juv_adult_capture_history_2009-2015_MARK_w_rings.txt",
             header=T,colClasses=c("character","factor","factor","factor","factor"))
MP_F_A <- 
  read.table("Data_files/Andava_MP_juv_adult_capture_history_2009-2015_MARK_w_rings.txt",
             header=T,colClasses=c("character","factor","factor","factor","factor"))
WfP_F_A <- 
  read.table("Data_files/Andava_WfP_juv_adult_capture_history_2009-2015_MARK_w_rings.txt",
             header=T,colClasses=c("character","factor","factor","factor","factor"))
Tuzla_F_A <- 
  read.table("Data_files/Tuzla_juv_adult_capture_history_1996-2004_MARK_w_rings.txt",
             header=T,colClasses=c("factor","character","factor","factor"))
Maio_F_A <- 
  read.table("Data_files/Maio_juv_adult_capture_history_2007-2015_MARK_w_rings.txt",
             header=T,colClasses=c("factor","character","factor","factor"))

setwd("/home/luke/comparative_ASR/Chick_survival_analysis/") 

Ceuta_C <- 
  read.table("Data_files/Ceuta_chick_capture_history_2006-2012_MARK_w_rings.txt",
             header=T,colClasses=c("factor","character","factor","integer",
                                   "numeric","factor","factor","factor",
                                   "factor","integer","factor"))
KiP_C <- 
  read.table("Data_files/Andava_KiP_chick_capture_history_2009-2015_MARK_w_rings.txt",
             header=T,colClasses=c("character","factor","factor","factor","integer","numeric",
                                   "factor","factor","integer"))
MP_C <- 
  read.table("Data_files/Andava_MP_chick_capture_history_2009-2015_MARK_w_rings.txt",
             header=T,colClasses=c("character","factor","factor","factor","integer","numeric",
                                   "factor","factor","integer"))
WfP_C <- 
  read.table("Data_files/Andava_WfP_chick_capture_history_2009-2015_MARK_w_rings.txt",
             header=T,colClasses=c("character","factor","factor","factor","integer","numeric",
                                   "factor","factor","integer"))
Maio_C <- 
  read.table("Data_files/Maio_chick_capture_history_2007-2015_MARK_w_rings.txt",
             header=T,colClasses=c("character","factor","factor","integer","numeric",
                                   "factor","factor","factor","integer"))
Tuzla_C <- 
  read.table("Data_files/Tuzla_chick_capture_history_1996-2000_MARK_w_rings.txt",
             header=T,colClasses=c("character","factor","factor","integer","numeric",
                                   "factor","factor"))

names(Ceuta_C)[1] <- "Ring"

plover_boot_all_pops <- function(Ceuta_F_A, Ceuta_C,
                                 KiP_F_A, KiP_C,
                                 MP_F_A, MP_C,
                                 WfP_F_A, WfP_C,
                                 Maio_F_A, Maio_C,
                                 Tuzla_F_A, Tuzla_C) {
  #Ceuta_C_boot <- Ceuta_C[sample(1:nrow(Ceuta_C), size = round(0.5 * nrow(Ceuta_C))), ]
  Ceuta_C_boot <- Ceuta_C[sample(1:nrow(Ceuta_C), size = nrow(Ceuta_C), replace = TRUE), ]
  Ceuta_pres <- Ceuta_F_A$ring %in% Ceuta_C_boot$Ring
  Ceuta_F_A_boot1 <- Ceuta_F_A[Ceuta_pres, ]
  spare_Ceuta_F_A <- Ceuta_F_A[!Ceuta_pres, ]
  Ceuta_F_A_boot2 <- spare_Ceuta_F_A[sample(1:nrow(spare_Ceuta_F_A), size = nrow(Ceuta_F_A) - nrow(Ceuta_F_A_boot1), replace = TRUE), ]
  Ceuta_F_A_boot <- rbind(Ceuta_F_A_boot1, Ceuta_F_A_boot2)
  
  KiP_C_boot <- KiP_C[sample(1:nrow(KiP_C), size = nrow(KiP_C), replace = TRUE), ]
  KiP_pres <- KiP_F_A$ring %in% KiP_C_boot$Ring
  KiP_F_A_boot1 <- KiP_F_A[KiP_pres, ]
  spare_KiP_F_A <- KiP_F_A[!KiP_pres, ]
  KiP_F_A_boot2 <- spare_KiP_F_A[sample(1:nrow(spare_KiP_F_A), size = nrow(KiP_F_A) - nrow(KiP_F_A_boot1), replace = TRUE), ]
  KiP_F_A_boot <- rbind(KiP_F_A_boot1, KiP_F_A_boot2)
  
  WfP_C_boot <- WfP_C[sample(1:nrow(WfP_C), size = nrow(WfP_C), replace = TRUE), ]
  WfP_pres <- WfP_F_A$ring %in% WfP_C_boot$Ring
  WfP_F_A_boot1 <- WfP_F_A[WfP_pres, ]
  spare_WfP_F_A <- WfP_F_A[!WfP_pres, ]
  WfP_F_A_boot2 <- spare_WfP_F_A[sample(1:nrow(spare_WfP_F_A), size = nrow(WfP_F_A) - nrow(WfP_F_A_boot1), replace = TRUE), ]
  WfP_F_A_boot <- rbind(WfP_F_A_boot1, WfP_F_A_boot2)
  
  MP_C_boot <- MP_C[sample(1:nrow(MP_C), size = nrow(MP_C), replace = TRUE), ]
  MP_pres <- MP_F_A$ring %in% MP_C_boot$Ring
  MP_F_A_boot1 <- MP_F_A[MP_pres, ]
  spare_MP_F_A <- MP_F_A[!MP_pres, ]
  MP_F_A_boot2 <- spare_MP_F_A[sample(1:nrow(spare_MP_F_A), size = nrow(MP_F_A) - nrow(MP_F_A_boot1), replace = TRUE), ]
  MP_F_A_boot <- rbind(MP_F_A_boot1, MP_F_A_boot2)
  
  Maio_C_boot <- Maio_C[sample(1:nrow(Maio_C), size = nrow(Maio_C), replace = TRUE), ]
  Maio_pres <- Maio_F_A$ring %in% Maio_C_boot$Ring
  Maio_F_A_boot1 <- Maio_F_A[Maio_pres, ]
  spare_Maio_F_A <- Maio_F_A[!Maio_pres, ]
  Maio_F_A_boot2 <- spare_Maio_F_A[sample(1:nrow(spare_Maio_F_A), size = nrow(Maio_F_A) - nrow(Maio_F_A_boot1), replace = TRUE), ]
  Maio_F_A_boot <- rbind(Maio_F_A_boot1, Maio_F_A_boot2)
  
  Tuzla_C_boot <- Tuzla_C[sample(1:nrow(Tuzla_C), size = nrow(Tuzla_C), replace = TRUE), ]
  Tuzla_pres <- Tuzla_F_A$ring %in% Tuzla_C_boot$Ring
  Tuzla_F_A_boot1 <- Tuzla_F_A[Tuzla_pres, ]
  spare_Tuzla_F_A <- Tuzla_F_A[!Tuzla_pres, ]
  Tuzla_F_A_boot2 <- spare_Tuzla_F_A[sample(1:nrow(spare_Tuzla_F_A), size = nrow(Tuzla_F_A) - nrow(Tuzla_F_A_boot1), replace = TRUE), ]
  Tuzla_F_A_boot <- rbind(Tuzla_F_A_boot1, Tuzla_F_A_boot2)
  
  out <- list(Ceuta_C_boot = Ceuta_C_boot, Ceuta_F_A_boot = Ceuta_F_A_boot,
              KiP_C_boot = KiP_C_boot, KiP_F_A_boot = KiP_F_A_boot,
              WfP_C_boot = WfP_C_boot, WfP_F_A_boot = WfP_F_A_boot,
              MP_C_boot = MP_C_boot, MP_F_A_boot = MP_F_A_boot,
              Maio_C_boot = Maio_C_boot, Maio_F_A_boot = Maio_F_A_boot,
              Tuzla_C_boot = Tuzla_C_boot, Tuzla_F_A_boot = Tuzla_F_A_boot)
}

calc_ASR_all_pops <- function(plover_boot_list) {
  Ceuta_C <- plover_boot_list[["Ceuta_C_boot"]]
  Ceuta_F_A <- plover_boot_list[["Ceuta_F_A_boot"]]
  
  KiP_C <- plover_boot_list[["KiP_C_boot"]]
  KiP_F_A <- plover_boot_list[["KiP_F_A_boot"]]
  
  WfP_C <- plover_boot_list[["WfP_C_boot"]]
  WfP_F_A <- plover_boot_list[["WfP_F_A_boot"]]
  
  MP_C <- plover_boot_list[["MP_C_boot"]]
  MP_F_A <- plover_boot_list[["MP_F_A_boot"]]
  
  Maio_C <- plover_boot_list[["Maio_C_boot"]]
  Maio_F_A <- plover_boot_list[["Maio_F_A_boot"]]
  
  Tuzla_C <- plover_boot_list[["Tuzla_C_boot"]]
  Tuzla_F_A <- plover_boot_list[["Tuzla_F_A_boot"]]
  
  # remove ring column
  Ceuta_F_A <- Ceuta_F_A[,-1]
  Ceuta_C <- Ceuta_C[,-1]
  
  KiP_F_A <- KiP_F_A[,-3]
  KiP_C <- KiP_C[,-2]
  
  WfP_F_A <- WfP_F_A[,-3]
  WfP_C <- WfP_C[,-2]
  
  MP_F_A <- MP_F_A[,-3]
  MP_C <- MP_C[,-2]
  
  Maio_F_A <- Maio_F_A[,-1]
  Maio_C <- Maio_C[,-2]
  
  Tuzla_F_A <- Tuzla_F_A[,-1]
  Tuzla_C <- Tuzla_C[,-2]
  
  # Remove capture histories that have no resights (i.e., all zeros in the ch)
  Ceuta_C <- Ceuta_C[which(str_detect(Ceuta_C[,"ch"],"1") == TRUE),]
  KiP_C <- KiP_C[which(str_detect(KiP_C[,"ch"],"1") == TRUE),]
  WfP_C <- WfP_C[which(str_detect(WfP_C[,"ch"],"1") == TRUE),]
  MP_C <- MP_C[which(str_detect(MP_C[,"ch"],"1") == TRUE),]
  Maio_C <- Maio_C[which(str_detect(Maio_C[,"ch"],"1") == TRUE),]
  Tuzla_C <- Tuzla_C[which(str_detect(Tuzla_C[,"ch"],"1") == TRUE),]
  
  # Create processed RMARK data format as CJS with 2 groups (sex and age 
  # initally ringed) and starting at year 2006
  Ceuta_F_A.proc=process.data(Ceuta_F_A,model="CJS",groups=c("sex","age"),
                              begin.time=2006,age.var=2,initial.age=c(1,0))
  KiP_F_A.proc=process.data(KiP_F_A,model="CJS",groups=c("sex","age"),
                               begin.time=2009,age.var=2,initial.age=c(1,0))
  WfP_F_A.proc=process.data(WfP_F_A,model="CJS",groups=c("sex","age"),
                              begin.time=2009,age.var=2,initial.age=c(1,0))
  MP_F_A.proc=process.data(MP_F_A,model="CJS",groups=c("sex","age"),
                               begin.time=2009,age.var=2,initial.age=c(1,0))
  Tuzla_F_A.proc=process.data(Tuzla_F_A,model="CJS",groups=c("sex","age"),
                             begin.time=1996,age.var=2,initial.age=c(1,0))
  Maio_F_A.proc=process.data(Maio_F_A,model="CJS",groups=c("sex","age"),
                            begin.time=2007,age.var=2,initial.age=c(1,0))
  # Create processed RMARK data format as CJS with 3 groups (sex, site, and 
  # cross-fostering treatment).
  Ceuta_C.proc=process.data(Ceuta_C,model="CJS",
                            groups=c("Sex","Care_site","CF","Year"))
  KiP_C.proc=process.data(KiP_C,model="CJS",
                               groups=c("Sex","Year"))
  MP_C.proc=process.data(MP_C,model="CJS",
                              groups=c("Sex","Year"))
  WfP_C.proc=process.data(WfP_C,model="CJS",
                               groups=c("Sex","Year"))
  Maio_C.proc=process.data(Maio_C,model="CJS",
                                groups=c("Sex","Site","Year"))
  Tuzla_C.proc=process.data(Tuzla_C,model="CJS",
                                 groups=c("Sex","Year"))
  # Create the design data
  Ceuta_F_A.ddl=make.design.data(Ceuta_F_A.proc)
  Ceuta_C.ddl=make.design.data(Ceuta_C.proc)
  
  KiP_F_A.ddl=make.design.data(KiP_F_A.proc)
  KiP_C.ddl=make.design.data(KiP_C.proc)
  
  WfP_F_A.ddl=make.design.data(WfP_F_A.proc)
  WfP_C.ddl=make.design.data(WfP_C.proc)
  
  MP_F_A.ddl=make.design.data(MP_F_A.proc)
  MP_C.ddl=make.design.data(MP_C.proc)
  
  Maio_F_A.ddl=make.design.data(Maio_F_A.proc)
  Maio_C.ddl=make.design.data(Maio_C.proc)
  
  Tuzla_F_A.ddl=make.design.data(Tuzla_F_A.proc)
  Tuzla_C.ddl=make.design.data(Tuzla_C.proc)

  # adds firstyear/adult age field to design data in column "age"
  Ceuta_F_A.ddl=add.design.data(Ceuta_F_A.proc,Ceuta_F_A.ddl,"Phi","age",bins=c(0,1,7),right=F,name="age",replace=T)
  KiP_F_A.ddl=add.design.data(KiP_F_A.proc,KiP_F_A.ddl,"Phi","age",bins=c(0,1,6),right=F,name="age",replace=T)
  WfP_F_A.ddl=add.design.data(WfP_F_A.proc,WfP_F_A.ddl,"Phi","age",bins=c(0,1,6),right=F,name="age",replace=T)
  MP_F_A.ddl=add.design.data(MP_F_A.proc,MP_F_A.ddl,"Phi","age",bins=c(0,1,6),right=F,name="age",replace=T)
  Tuzla_F_A.ddl=add.design.data(Tuzla_F_A.proc,Tuzla_F_A.ddl,"Phi","age",bins=c(0,1,6),right=F,name="age",replace=T)
  Maio_F_A.ddl=add.design.data(Maio_F_A.proc,Maio_F_A.ddl,"Phi","age",bins=c(0,1,8),right=F,name="age",replace=T)
  
  
  # create a dummy field called marked.as.adult which is 0 for the group initally ringed as juvenile and 1 for the group marked as adults.
  Ceuta_F_A.ddl$Phi$marked.as.adult=0
  Ceuta_F_A.ddl$Phi$marked.as.adult[Ceuta_F_A.ddl$Phi$initial.age.class=="A"]=1 
  Ceuta_F_A.ddl$p$marked.as.adult=0
  Ceuta_F_A.ddl$p$marked.as.adult[Ceuta_F_A.ddl$p$initial.age.class=="A"]=1
  
  KiP_F_A.ddl$Phi$marked.as.adult=0
  KiP_F_A.ddl$Phi$marked.as.adult[KiP_F_A.ddl$Phi$initial.age.class=="A"]=1 
  KiP_F_A.ddl$p$marked.as.adult=0
  KiP_F_A.ddl$p$marked.as.adult[KiP_F_A.ddl$p$initial.age.class=="A"]=1
  
  WfP_F_A.ddl$Phi$marked.as.adult=0
  WfP_F_A.ddl$Phi$marked.as.adult[WfP_F_A.ddl$Phi$initial.age.class=="A"]=1 
  WfP_F_A.ddl$p$marked.as.adult=0
  WfP_F_A.ddl$p$marked.as.adult[WfP_F_A.ddl$p$initial.age.class=="A"]=1
  
  MP_F_A.ddl$Phi$marked.as.adult=0
  MP_F_A.ddl$Phi$marked.as.adult[MP_F_A.ddl$Phi$initial.age.class=="A"]=1 
  MP_F_A.ddl$p$marked.as.adult=0
  MP_F_A.ddl$p$marked.as.adult[MP_F_A.ddl$p$initial.age.class=="A"]=1
  
  Maio_F_A.ddl$Phi$marked.as.adult=0
  Maio_F_A.ddl$Phi$marked.as.adult[Maio_F_A.ddl$Phi$initial.age.class=="A"]=1 
  Maio_F_A.ddl$p$marked.as.adult=0
  Maio_F_A.ddl$p$marked.as.adult[Maio_F_A.ddl$p$initial.age.class=="A"]=1
  
  Tuzla_F_A.ddl$Phi$marked.as.adult=0
  Tuzla_F_A.ddl$Phi$marked.as.adult[Tuzla_F_A.ddl$Phi$initial.age.class=="A"]=1 
  Tuzla_F_A.ddl$p$marked.as.adult=0
  Tuzla_F_A.ddl$p$marked.as.adult[Tuzla_F_A.ddl$p$initial.age.class=="A"]=1
  
  # check parameter matrices to see if groups were binned correctly
  #PIMS(mark(Ceuta_F_A.proc,Ceuta_F_A.ddl,model.parameters=list(Phi=list(formula=~age+sex)),output=F),"Phi")
  
  # Create quadratic time variable so that it can be tested along side the annual models
  time <- c(0:(Ceuta_F_A.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  quad_time$time <- c(2006:2012)
  Ceuta_F_A.ddl$p=merge_design.covariates(Ceuta_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  Ceuta_F_A.ddl$Phi=merge_design.covariates(Ceuta_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  time <- c(0:(KiP_F_A.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  quad_time$time <- c(2009:2015)
  KiP_F_A.ddl$p=merge_design.covariates(KiP_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  KiP_F_A.ddl$Phi=merge_design.covariates(KiP_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  time <- c(0:(WfP_F_A.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  quad_time$time <- c(2009:2015)
  WfP_F_A.ddl$p=merge_design.covariates(WfP_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  WfP_F_A.ddl$Phi=merge_design.covariates(WfP_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  time <- c(0:(MP_F_A.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  quad_time$time <- c(2009:2015)
  MP_F_A.ddl$p=merge_design.covariates(MP_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  MP_F_A.ddl$Phi=merge_design.covariates(MP_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  time <- c(0:(Maio_F_A.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  quad_time$time <- c(2007:2015)
  Maio_F_A.ddl$p=merge_design.covariates(Maio_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  Maio_F_A.ddl$Phi=merge_design.covariates(Maio_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  time <- c(0:(Tuzla_F_A.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  quad_time$time <- c(1996:2000,2014)
  Tuzla_F_A.ddl$p=merge_design.covariates(Tuzla_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  Tuzla_F_A.ddl$Phi=merge_design.covariates(Tuzla_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  

  time <- c(0:(Ceuta_C.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  Ceuta_C.ddl$p=merge_design.covariates(Ceuta_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  Ceuta_C.ddl$Phi=merge_design.covariates(Ceuta_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  time <- c(0:(KiP_C.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  KiP_C.ddl$p=merge_design.covariates(KiP_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  KiP_C.ddl$Phi=merge_design.covariates(KiP_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  time <- c(0:(MP_C.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  MP_C.ddl$p=merge_design.covariates(MP_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  MP_C.ddl$Phi=merge_design.covariates(MP_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  time <- c(0:(WfP_C.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  WfP_C.ddl$p=merge_design.covariates(WfP_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  WfP_C.ddl$Phi=merge_design.covariates(WfP_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  time <- c(0:(Maio_C.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  Maio_C.ddl$p=merge_design.covariates(Maio_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  Maio_C.ddl$Phi=merge_design.covariates(Maio_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  time <- c(0:(Tuzla_C.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  Tuzla_C.ddl$p=merge_design.covariates(Tuzla_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  Tuzla_C.ddl$Phi=merge_design.covariates(Tuzla_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  # Fledgling and adult survival models:
  # Ceuta:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/Ceuta/Juv_adult/") # set wd so that results go to the correct folder
  Ceuta_Phi_F_A <- list(formula=~age*sex) 
  Ceuta_p_F_A <- list(formula=~age*time)
  Ceuta_model_F_A <- 
    mark(Ceuta_F_A.proc,Ceuta_F_A.ddl, 
         model.parameters = list(Phi = Ceuta_Phi_F_A, p = Ceuta_p_F_A),
         threads = 8, output = FALSE)
  
  # KiP:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/KiP/Juv_adult/") # set wd so that results go to the correct folder
  KiP_Phi_F_A <- list(formula=~age*sex) 
  KiP_p_F_A <- list(formula=~age*time)
  KiP_model_F_A <- 
    mark(KiP_F_A.proc,KiP_F_A.ddl, 
         model.parameters = list(Phi = KiP_Phi_F_A, p = KiP_p_F_A),
         threads = 8, output = FALSE)
  
  # WfP:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/WfP/Juv_adult/") # set wd so that results go to the correct folder
  WfP_Phi_F_A <- list(formula=~age*sex) 
  WfP_p_F_A <- list(formula=~~time + age + sex)
  WfP_model_F_A <- 
    mark(WfP_F_A.proc,WfP_F_A.ddl, 
         model.parameters = list(Phi = WfP_Phi_F_A, p = WfP_p_F_A),
         threads = 8, output = FALSE)
  
  # MP:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/MP/Juv_adult/") # set wd so that results go to the correct folder
  MP_Phi_F_A <- list(formula=~age*sex) 
  MP_p_F_A <- list(formula=~~age + time)
  MP_model_F_A <- 
    mark(MP_F_A.proc,MP_F_A.ddl, 
         model.parameters = list(Phi = MP_Phi_F_A, p = MP_p_F_A),
         threads = 8, output = FALSE)
  
  # Maio:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/Maio/Juv_adult/") # set wd so that results go to the correct folder
  Maio_Phi_F_A <- list(formula=~age*sex) 
  Maio_p_F_A <- list(formula=~age * time)
  Maio_model_F_A <- 
    mark(Maio_F_A.proc,Maio_F_A.ddl, 
         model.parameters = list(Phi = Maio_Phi_F_A, p = Maio_p_F_A),
         threads = 8, output = FALSE)
  
  # Tuzla:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/Tuzla/Juv_adult/") # set wd so that results go to the correct folder
  Tuzla_Phi_F_A <- list(formula=~age*sex) 
  Tuzla_p_F_A <- list(formula=~time * age * sex)
  Tuzla_model_F_A <- 
    mark(Tuzla_F_A.proc,Tuzla_F_A.ddl, 
         model.parameters = list(Phi = Tuzla_Phi_F_A, p = Tuzla_p_F_A),
         threads = 8, output = FALSE)
  
  # Chick survival model:
  # Ceuta:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/Ceuta/Chick/") # set wd so that results go to the correct folder
  Ceuta_Phi_C=list(formula=~Sex*Quadratic)
  Ceuta_p_C=list(formula=~Year*Quadratic)
  Ceuta_model_C <- 
    mark(Ceuta_C.proc,Ceuta_C.ddl, 
         model.parameters = list(Phi = Ceuta_Phi_C, p = Ceuta_p_C),
         threads = 8, output = FALSE)
  
  # KiP:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/KiP/Chick/") # set wd so that results go to the correct folder
  KiP_Phi_C=list(formula=~~Sex * Time)
  KiP_p_C=list(formula=~Year*Quadratic)
  KiP_model_C <- 
    mark(KiP_C.proc,KiP_C.ddl, 
         model.parameters = list(Phi = KiP_Phi_C, p = KiP_p_C),
         threads = 8, output = FALSE)
  
  # WfP:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/WfP/Chick/") # set wd so that results go to the correct folder
  WfP_Phi_C=list(formula=~Sex * Cubic)
  WfP_p_C=list(formula=~Year * Time)
  WfP_model_C <- 
    mark(WfP_C.proc,WfP_C.ddl, 
         model.parameters = list(Phi = WfP_Phi_C, p = WfP_p_C),
         threads = 8, output = FALSE)
  
  # MP:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/MP/Chick/") # set wd so that results go to the correct folder
  MP_Phi_C=list(formula=~Sex)
  MP_p_C=list(formula=~~Year + Quadratic)
  MP_model_C <- 
    mark(MP_C.proc,MP_C.ddl, 
         model.parameters = list(Phi = MP_Phi_C, p = MP_p_C),
         threads = 8, output = FALSE)
  
  # Maio:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/Maio/Chick/") # set wd so that results go to the correct folder
  Maio_Phi_C=list(formula=~Sex * Cubic)
  Maio_p_C=list(formula=~Sex * Year * Cubic)
  Maio_model_C <- 
    mark(Maio_C.proc,Maio_C.ddl, 
         model.parameters = list(Phi = Maio_Phi_C, p = Maio_p_C),
         threads = 8, output = FALSE)
  
  # Tuzla:
  setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/Tuzla/Chick/") # set wd so that results go to the correct folder
  Tuzla_Phi_C=list(formula=~Sex * Quadratic)
  Tuzla_p_C=list(formula=~Year + Time)
  Tuzla_model_C <- 
    mark(Tuzla_C.proc,Tuzla_C.ddl, 
         model.parameters = list(Phi = Tuzla_Phi_C, p = Tuzla_p_C),
         threads = 8, output = FALSE)
  
  # extract and format survival rates from chick model output
  # Ceuta:
  Ceuta_C_reals <- Ceuta_model_C$results$real
  Groups <- data.frame(str_split_fixed(rownames(Ceuta_C_reals), " ", n = 5))
  Ceuta_C_reals <- cbind(Groups, Ceuta_C_reals)
  Ceuta_C_reals <- Ceuta_C_reals[which(Ceuta_C_reals$X1 == "Phi"),]
  Ceuta_C_reals$Sex <- unlist(str_extract_all(Ceuta_C_reals$X2,"[FM]"))
  Ceuta_C_reals$Sex <- as.factor(ifelse(Ceuta_C_reals$Sex == "F","Female","Male"))
  Ceuta_Survival_to_Fledge_F <- 
    prod(Ceuta_C_reals[which(Ceuta_C_reals$Sex == "Female"),
                        c("estimate")])
  Ceuta_Survival_to_Fledge_M <- 
    prod(Ceuta_C_reals[which(Ceuta_C_reals$Sex == "Male"),
                        c("estimate")])
  
  estimate <- c(Ceuta_Survival_to_Fledge_F, Ceuta_Survival_to_Fledge_M)
  Sex <- c("Female", "Male")
  age <- c("Chick", "Chick")
  Sex_Age <- paste(Sex,age,sep = "_")
  Ceuta_sex_chick_survival <- data.frame(Sex_Age, estimate)
  Ceuta_sex_chick_survival$species <- "Snowy"
  
  # extract and format survival rates from fledgling and adult model output
  Ceuta_F_A_reals <- Ceuta_model_F_A$results$real
  Groups <- data.frame(str_split_fixed(rownames(Ceuta_F_A_reals), " ", n = 5))
  Ceuta_F_A_reals <- cbind(Groups, Ceuta_F_A_reals)
  Ceuta_F_A_reals <- Ceuta_F_A_reals[which(Ceuta_F_A_reals$X1 == "Phi"),]
  Ceuta_F_A_reals$age <- unlist(str_extract_all(Ceuta_F_A_reals$X2,"[AJ]"))
  Ceuta_F_A_reals$age <- as.factor(ifelse(Ceuta_F_A_reals$age == "A","Adult","Juvenile"))
  Ceuta_F_A_reals$Sex <- unlist(str_extract_all(Ceuta_F_A_reals$X2,"[FM]"))
  Ceuta_F_A_reals$Sex <- as.factor(ifelse(Ceuta_F_A_reals$Sex == "F","Female","Male"))
  Ceuta_F_A_reals$Sex_Age <- paste(Ceuta_F_A_reals$Sex,Ceuta_F_A_reals$age,sep = "_")
  Ceuta_Survival_rates <- Ceuta_F_A_reals[,c("Sex_Age", "estimate")]
  Ceuta_Survival_rates$species <- "Snowy"
  row.names(Ceuta_Survival_rates) <- NULL
  Ceuta_Survival_rates <- rbind(Ceuta_Survival_rates, Ceuta_sex_chick_survival)
  
  # KiP:
  KiP_C_reals <- KiP_model_C$results$real
  Groups <- data.frame(str_split_fixed(rownames(KiP_C_reals), " ", n = 5))
  KiP_C_reals <- cbind(Groups, KiP_C_reals)
  KiP_C_reals <- KiP_C_reals[which(KiP_C_reals$X1 == "Phi"),]
  KiP_C_reals$Sex <- unlist(str_extract_all(KiP_C_reals$X2,"[FM]"))
  KiP_C_reals$Sex <- as.factor(ifelse(KiP_C_reals$Sex == "F","Female","Male"))
  KiP_Survival_to_Fledge_F <- 
    prod(KiP_C_reals[which(KiP_C_reals$Sex == "Female"),
                       c("estimate")])
  KiP_Survival_to_Fledge_M <- 
    prod(KiP_C_reals[which(KiP_C_reals$Sex == "Male"),
                       c("estimate")])
  
  estimate <- c(KiP_Survival_to_Fledge_F, KiP_Survival_to_Fledge_M)
  Sex <- c("Female", "Male")
  age <- c("Chick", "Chick")
  Sex_Age <- paste(Sex,age,sep = "_")
  KiP_sex_chick_survival <- data.frame(Sex_Age, estimate)
  KiP_sex_chick_survival$species <- "Kittlitz's"
  
  # extract and format survival rates from fledgling and adult model output
  KiP_F_A_reals <- KiP_model_F_A$results$real
  Groups <- data.frame(str_split_fixed(rownames(KiP_F_A_reals), " ", n = 5))
  KiP_F_A_reals <- cbind(Groups, KiP_F_A_reals)
  KiP_F_A_reals <- KiP_F_A_reals[which(KiP_F_A_reals$X1 == "Phi"),]
  KiP_F_A_reals$age <- unlist(str_extract_all(KiP_F_A_reals$X2,"[AJ]"))
  KiP_F_A_reals$age <- as.factor(ifelse(KiP_F_A_reals$age == "A","Adult","Juvenile"))
  KiP_F_A_reals$Sex <- unlist(str_extract_all(KiP_F_A_reals$X2,"[FM]"))
  KiP_F_A_reals$Sex <- as.factor(ifelse(KiP_F_A_reals$Sex == "F","Female","Male"))
  KiP_F_A_reals$Sex_Age <- paste(KiP_F_A_reals$Sex,KiP_F_A_reals$age,sep = "_")
  KiP_Survival_rates <- KiP_F_A_reals[,c("Sex_Age", "estimate")]
  KiP_Survival_rates$species <- "Kittlitz's"
  row.names(KiP_Survival_rates) <- NULL
  KiP_Survival_rates <- rbind(KiP_Survival_rates, KiP_sex_chick_survival)
  
  # WfP:
  WfP_C_reals <- WfP_model_C$results$real
  Groups <- data.frame(str_split_fixed(rownames(WfP_C_reals), " ", n = 5))
  WfP_C_reals <- cbind(Groups, WfP_C_reals)
  WfP_C_reals <- WfP_C_reals[which(WfP_C_reals$X1 == "Phi"),]
  WfP_C_reals$Sex <- unlist(str_extract_all(WfP_C_reals$X2,"[FM]"))
  WfP_C_reals$Sex <- as.factor(ifelse(WfP_C_reals$Sex == "F","Female","Male"))
  WfP_Survival_to_Fledge_F <- 
    prod(WfP_C_reals[which(WfP_C_reals$Sex == "Female"),
                       c("estimate")])
  WfP_Survival_to_Fledge_M <- 
    prod(WfP_C_reals[which(WfP_C_reals$Sex == "Male"),
                       c("estimate")])
  
  estimate <- c(WfP_Survival_to_Fledge_F, WfP_Survival_to_Fledge_M)
  Sex <- c("Female", "Male")
  age <- c("Chick", "Chick")
  Sex_Age <- paste(Sex,age,sep = "_")
  WfP_sex_chick_survival <- data.frame(Sex_Age, estimate)
  WfP_sex_chick_survival$species <- "White-fronted"
  
  # extract and format survival rates from fledgling and adult model output
  WfP_F_A_reals <- WfP_model_F_A$results$real
  Groups <- data.frame(str_split_fixed(rownames(WfP_F_A_reals), " ", n = 5))
  WfP_F_A_reals <- cbind(Groups, WfP_F_A_reals)
  WfP_F_A_reals <- WfP_F_A_reals[which(WfP_F_A_reals$X1 == "Phi"),]
  WfP_F_A_reals$age <- unlist(str_extract_all(WfP_F_A_reals$X2,"[AJ]"))
  WfP_F_A_reals$age <- as.factor(ifelse(WfP_F_A_reals$age == "A","Adult","Juvenile"))
  WfP_F_A_reals$Sex <- unlist(str_extract_all(WfP_F_A_reals$X2,"[FM]"))
  WfP_F_A_reals$Sex <- as.factor(ifelse(WfP_F_A_reals$Sex == "F","Female","Male"))
  WfP_F_A_reals$Sex_Age <- paste(WfP_F_A_reals$Sex,WfP_F_A_reals$age,sep = "_")
  WfP_Survival_rates <- WfP_F_A_reals[,c("Sex_Age", "estimate")]
  WfP_Survival_rates$species <- "White-fronted"
  row.names(WfP_Survival_rates) <- NULL
  WfP_Survival_rates <- rbind(WfP_Survival_rates, WfP_sex_chick_survival)
  
  # MP:
  MP_C_reals <- MP_model_C$results$real
  Groups <- data.frame(str_split_fixed(rownames(MP_C_reals), " ", n = 5))
  MP_C_reals <- cbind(Groups, MP_C_reals)
  MP_C_reals <- MP_C_reals[which(MP_C_reals$X1 == "Phi"),]
  MP_C_reals$Sex <- unlist(str_extract_all(MP_C_reals$X2,"[FM]"))
  MP_C_reals$Sex <- as.factor(ifelse(MP_C_reals$Sex == "F","Female","Male"))
  MP_Survival_to_Fledge_F <- 
    MP_C_reals[which(MP_C_reals$Sex == "Female"),
                       c("estimate")]^30
  MP_Survival_to_Fledge_M <- 
    MP_C_reals[which(MP_C_reals$Sex == "Male"),
                       c("estimate")]^30
  
  estimate <- c(MP_Survival_to_Fledge_F, MP_Survival_to_Fledge_M)
  Sex <- c("Female", "Male")
  age <- c("Chick", "Chick")
  Sex_Age <- paste(Sex,age,sep = "_")
  MP_sex_chick_survival <- data.frame(Sex_Age, estimate)
  MP_sex_chick_survival$species <- "Madagascar"
  
  # extract and format survival rates from fledgling and adult model output
  MP_F_A_reals <- MP_model_F_A$results$real
  Groups <- data.frame(str_split_fixed(rownames(MP_F_A_reals), " ", n = 5))
  MP_F_A_reals <- cbind(Groups, MP_F_A_reals)
  MP_F_A_reals <- MP_F_A_reals[which(MP_F_A_reals$X1 == "Phi"),]
  MP_F_A_reals$age <- unlist(str_extract_all(MP_F_A_reals$X2,"[AJ]"))
  MP_F_A_reals$age <- as.factor(ifelse(MP_F_A_reals$age == "A","Adult","Juvenile"))
  MP_F_A_reals$Sex <- unlist(str_extract_all(MP_F_A_reals$X2,"[FM]"))
  MP_F_A_reals$Sex <- as.factor(ifelse(MP_F_A_reals$Sex == "F","Female","Male"))
  MP_F_A_reals$Sex_Age <- paste(MP_F_A_reals$Sex,MP_F_A_reals$age,sep = "_")
  MP_Survival_rates <- MP_F_A_reals[,c("Sex_Age", "estimate")]
  MP_Survival_rates$species <- "Madagascar"
  row.names(MP_Survival_rates) <- NULL
  MP_Survival_rates <- rbind(MP_Survival_rates, MP_sex_chick_survival)
  
  # Maio:
  Maio_C_reals <- Maio_model_C$results$real
  Groups <- data.frame(str_split_fixed(rownames(Maio_C_reals), " ", n = 5))
  Maio_C_reals <- cbind(Groups, Maio_C_reals)
  Maio_C_reals <- Maio_C_reals[which(Maio_C_reals$X1 == "Phi"),]
  Maio_C_reals$Sex <- unlist(str_extract_all(Maio_C_reals$X2,"[FM]"))
  Maio_C_reals$Sex <- as.factor(ifelse(Maio_C_reals$Sex == "F","Female","Male"))
  Maio_Survival_to_Fledge_F <- 
    prod(Maio_C_reals[which(Maio_C_reals$Sex == "Female"),
                       c("estimate")])
  Maio_Survival_to_Fledge_M <- 
    prod(Maio_C_reals[which(Maio_C_reals$Sex == "Male"),
                       c("estimate")])
  
  estimate <- c(Maio_Survival_to_Fledge_F, Maio_Survival_to_Fledge_M)
  Sex <- c("Female", "Male")
  age <- c("Chick", "Chick")
  Sex_Age <- paste(Sex,age,sep = "_")
  Maio_sex_chick_survival <- data.frame(Sex_Age, estimate)
  Maio_sex_chick_survival$species <- "Kentish (Maio)"
  
  # extract and format survival rates from fledgling and adult model output
  Maio_F_A_reals <- Maio_model_F_A$results$real
  Groups <- data.frame(str_split_fixed(rownames(Maio_F_A_reals), " ", n = 5))
  Maio_F_A_reals <- cbind(Groups, Maio_F_A_reals)
  Maio_F_A_reals <- Maio_F_A_reals[which(Maio_F_A_reals$X1 == "Phi"),]
  Maio_F_A_reals$age <- unlist(str_extract_all(Maio_F_A_reals$X2,"[AJ]"))
  Maio_F_A_reals$age <- as.factor(ifelse(Maio_F_A_reals$age == "A","Adult","Juvenile"))
  Maio_F_A_reals$Sex <- unlist(str_extract_all(Maio_F_A_reals$X2,"[FM]"))
  Maio_F_A_reals$Sex <- as.factor(ifelse(Maio_F_A_reals$Sex == "F","Female","Male"))
  Maio_F_A_reals$Sex_Age <- paste(Maio_F_A_reals$Sex,Maio_F_A_reals$age,sep = "_")
  Maio_Survival_rates <- Maio_F_A_reals[,c("Sex_Age", "estimate")]
  Maio_Survival_rates$species <- "Kentish (Maio)"
  row.names(Maio_Survival_rates) <- NULL
  Maio_Survival_rates <- rbind(Maio_Survival_rates, Maio_sex_chick_survival)
  
  # Tuzla:
  Tuzla_C_reals <- Tuzla_model_C$results$real
  Groups <- data.frame(str_split_fixed(rownames(Tuzla_C_reals), " ", n = 5))
  Tuzla_C_reals <- cbind(Groups, Tuzla_C_reals)
  Tuzla_C_reals <- Tuzla_C_reals[which(Tuzla_C_reals$X1 == "Phi"),]
  Tuzla_C_reals$Sex <- unlist(str_extract_all(Tuzla_C_reals$X2,"[FM]"))
  Tuzla_C_reals$Sex <- as.factor(ifelse(Tuzla_C_reals$Sex == "F","Female","Male"))
  Tuzla_Survival_to_Fledge_F <- 
    prod(Tuzla_C_reals[which(Tuzla_C_reals$Sex == "Female"),
                       c("estimate")])
  Tuzla_Survival_to_Fledge_M <- 
    prod(Tuzla_C_reals[which(Tuzla_C_reals$Sex == "Male"),
                       c("estimate")])
  
  estimate <- c(Tuzla_Survival_to_Fledge_F, Tuzla_Survival_to_Fledge_M)
  Sex <- c("Female", "Male")
  age <- c("Chick", "Chick")
  Sex_Age <- paste(Sex,age,sep = "_")
  Tuzla_sex_chick_survival <- data.frame(Sex_Age, estimate)
  Tuzla_sex_chick_survival$species <- "Kentish (Tuzla)"
  
  # extract and format survival rates from fledgling and adult model output
  Tuzla_F_A_reals <- Tuzla_model_F_A$results$real
  Groups <- data.frame(str_split_fixed(rownames(Tuzla_F_A_reals), " ", n = 5))
  Tuzla_F_A_reals <- cbind(Groups, Tuzla_F_A_reals)
  Tuzla_F_A_reals <- Tuzla_F_A_reals[which(Tuzla_F_A_reals$X1 == "Phi"),]
  Tuzla_F_A_reals$age <- unlist(str_extract_all(Tuzla_F_A_reals$X2,"[AJ]"))
  Tuzla_F_A_reals$age <- as.factor(ifelse(Tuzla_F_A_reals$age == "A","Adult","Juvenile"))
  Tuzla_F_A_reals$Sex <- unlist(str_extract_all(Tuzla_F_A_reals$X2,"[FM]"))
  Tuzla_F_A_reals$Sex <- as.factor(ifelse(Tuzla_F_A_reals$Sex == "F","Female","Male"))
  Tuzla_F_A_reals$Sex_Age <- paste(Tuzla_F_A_reals$Sex,Tuzla_F_A_reals$age,sep = "_")
  Tuzla_Survival_rates <- Tuzla_F_A_reals[,c("Sex_Age", "estimate")]
  Tuzla_Survival_rates$species <- "Kentish (Tuzla)"
  row.names(Tuzla_Survival_rates) <- NULL
  Tuzla_Survival_rates <- rbind(Tuzla_Survival_rates, Tuzla_sex_chick_survival)
  
  # Define Ceuta vital rates estimated from mark-recapture analysis:
  Ceuta_VR <- list(F_Chk_survl = Ceuta_Survival_rates[5,2],
                   F_Juv_survl = Ceuta_Survival_rates[3,2],
                   F_Adt_survl = Ceuta_Survival_rates[1,2],
                   M_Chk_survl = Ceuta_Survival_rates[6,2],
                   M_Juv_survl = Ceuta_Survival_rates[4,2],
                   M_Adt_survl = Ceuta_Survival_rates[2,2],
                   # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                   h = 1,
                   k = 3,
                   # Define primary sex ratio (assumed to be 0.5)
                   PSR = 0.5)
  
  # Ceuta matrix:
  Ceuta_matrix <- plover_matrix(Ceuta_VR)
  
  # Determine the ASR at the stable stage distribution (assume h is 1, which 
  # is monogamy)
  Ceuta_ASR_h_1 <- freq_dep_SSD_ASR(A = Ceuta_matrix, h = 1, k = 3)
  
  # Extract ASR
  Ceuta_ASR <- Ceuta_ASR_h_1$ASR
  
  # Define KiP vital rates estimated from mark-recapture analysis:
  KiP_VR <- list(#F_Chk_survl = KiP_Survival_rates[5,2],
                   F_Juv_survl = KiP_Survival_rates[3,2],
                   F_Adt_survl = KiP_Survival_rates[1,2],
                   #M_Chk_survl = KiP_Survival_rates[6,2],
                   M_Juv_survl = KiP_Survival_rates[4,2],
                   M_Adt_survl = KiP_Survival_rates[2,2],
                   # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                   h = 1,
                   k = 2,
                   # Define primary sex ratio (assumed to be 0.5)
                   PSR = 0.5)
  
  # KiP matrix:
  KiP_matrix <- plover_matrix(KiP_VR, chick_surv = FALSE)
  
  # Determine the ASR at the stable stage distribution (assume h is 1, which 
  # is monogamy)
  KiP_ASR_h_1 <- freq_dep_SSD_ASR(A = KiP_matrix, h = 1, k = 2)
  
  # Extract ASR
  KiP_ASR <- KiP_ASR_h_1$ASR
  
  # Define WfP vital rates estimated from mark-recapture analysis:
  WfP_VR <- list(#F_Chk_survl = WfP_Survival_rates[5,2],
                   F_Juv_survl = WfP_Survival_rates[3,2],
                   F_Adt_survl = WfP_Survival_rates[1,2],
                   #M_Chk_survl = WfP_Survival_rates[6,2],
                   M_Juv_survl = WfP_Survival_rates[4,2],
                   M_Adt_survl = WfP_Survival_rates[2,2],
                   # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                   h = 1,
                   k = 2,
                   # Define primary sex ratio (assumed to be 0.5)
                   PSR = 0.5)
  
  # WfP matrix:
  WfP_matrix <- plover_matrix(WfP_VR, chick_surv = FALSE)
  
  # Determine the ASR at the stable stage distribution (assume h is 1, which 
  # is monogamy)
  WfP_ASR_h_1 <- freq_dep_SSD_ASR(A = WfP_matrix, h = 1, k = 2)
  
  # Extract ASR
  WfP_ASR <- WfP_ASR_h_1$ASR
  
  # Define MP vital rates estimated from mark-recapture analysis:
  MP_VR <- list(#F_Chk_survl = MP_Survival_rates[5,2],
                   F_Juv_survl = MP_Survival_rates[3,2],
                   F_Adt_survl = MP_Survival_rates[1,2],
                   #M_Chk_survl = MP_Survival_rates[6,2],
                   M_Juv_survl = MP_Survival_rates[4,2],
                   M_Adt_survl = MP_Survival_rates[2,2],
                   # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                   h = 1,
                   k = 2,
                   # Define primary sex ratio (assumed to be 0.5)
                   PSR = 0.5)
  
  # MP matrix:
  MP_matrix <- plover_matrix(MP_VR, chick_surv = FALSE)
  
  # Determine the ASR at the stable stage distribution (assume h is 1, which 
  # is monogamy)
  MP_ASR_h_1 <- freq_dep_SSD_ASR(A = MP_matrix, h = 1, k = 2)
  
  # Extract ASR
  MP_ASR <- MP_ASR_h_1$ASR
  
  # Define Maio vital rates estimated from mark-recapture analysis:
  Maio_VR <- list(#F_Chk_survl = Maio_Survival_rates[5,2],
                   F_Juv_survl = Maio_Survival_rates[3,2],
                   F_Adt_survl = Maio_Survival_rates[1,2],
                   #M_Chk_survl = Maio_Survival_rates[6,2],
                   M_Juv_survl = Maio_Survival_rates[4,2],
                   M_Adt_survl = Maio_Survival_rates[2,2],
                   # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                   h = 1,
                   k = 3,
                   # Define primary sex ratio (assumed to be 0.5)
                   PSR = 0.5)
  
  # Maio matrix:
  Maio_matrix <- plover_matrix(Maio_VR, chick_surv = FALSE)
  
  # Determine the ASR at the stable stage distribution (assume h is 1, which 
  # is monogamy)
  Maio_ASR_h_1 <- freq_dep_SSD_ASR(A = Maio_matrix, h = 1, k = 3)
  
  # Extract ASR
  Maio_ASR <- Maio_ASR_h_1$ASR
  
  # Define Tuzla vital rates estimated from mark-recapture analysis:
  Tuzla_VR <- list(#F_Chk_survl = Tuzla_Survival_rates[5,2],
                   F_Juv_survl = Tuzla_Survival_rates[3,2],
                   F_Adt_survl = Tuzla_Survival_rates[1,2],
                   #M_Chk_survl = Tuzla_Survival_rates[6,2],
                   M_Juv_survl = Tuzla_Survival_rates[4,2],
                   M_Adt_survl = Tuzla_Survival_rates[2,2],
                   # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                   h = 1,
                   k = 3,
                   # Define primary sex ratio (assumed to be 0.5)
                   PSR = 0.5)
  
  # Tuzla matrix:
  Tuzla_matrix <- plover_matrix(Tuzla_VR, chick_surv = FALSE)
  
  # Determine the ASR at the stable stage distribution (assume h is 1, which 
  # is monogamy)
  Tuzla_ASR_h_1 <- freq_dep_SSD_ASR(A = Tuzla_matrix, h = 1, k = 3)
  
  # Extract ASR
  Tuzla_ASR <- Tuzla_ASR_h_1$ASR
  
  #
  ASR <- c(Ceuta_ASR, KiP_ASR, WfP_ASR, MP_ASR, Maio_ASR, Tuzla_ASR)
}

# bootstrapping
run_boot_all_pops <- function(Ceuta_F_A, Ceuta_C,
                              KiP_F_A, KiP_C,
                              MP_F_A, MP_C,
                              WfP_F_A, WfP_C,
                              Maio_F_A, Maio_C,
                              Tuzla_F_A, Tuzla_C){
  plover_list_all_pops <- plover_boot_all_pops(Ceuta_F_A, Ceuta_C,
                                               KiP_F_A, KiP_C,
                                               MP_F_A, MP_C,
                                               WfP_F_A, WfP_C,
                                               Maio_F_A, Maio_C,
                                               Tuzla_F_A, Tuzla_C)
  ASR <- calc_ASR_all_pops(plover_list_all_pops)
}

# All_pops_boot_test_w_2 <- replicate(2, run_boot_all_pops(Ceuta_F_A, Ceuta_C,
#                                                         KiP_F_A, KiP_C,
#                                                         MP_F_A, MP_C,
#                                                         WfP_F_A, WfP_C,
#                                                         Maio_F_A, Maio_C,
#                                                         Tuzla_F_A, Tuzla_C))

# 1000 bootstraps! 
# All_pops_boot_1000_March_11_2016 <- replicate(1000, run_boot_all_pops(Ceuta_F_A, Ceuta_C,
#                                                                       KiP_F_A, KiP_C,
#                                                                       MP_F_A, MP_C,
#                                                                       WfP_F_A, WfP_C,
#                                                                       Maio_F_A, Maio_C,
#                                                                       Tuzla_F_A, Tuzla_C))
All_pops_boot_1000_March_13_2016 <- replicate(1000, run_boot_all_pops(Ceuta_F_A, Ceuta_C,
                                                                      KiP_F_A, KiP_C,
                                                                      MP_F_A, MP_C,
                                                                      WfP_F_A, WfP_C,
                                                                      Maio_F_A, Maio_C,
                                                                      Tuzla_F_A, Tuzla_C))
write.table(All_pops_boot_1000_March_13_2016, file = "/home/luke/ceuta_ASR/bootstrap_products/All_pops_boot_1000_March_13_2016.txt", sep = "\t")
all_ASRs_w_replace_1000 <- as.data.frame(t(All_pops_boot_1000_March_13_2016))
colnames(all_ASRs_w_replace_1000) <- c("Snowy", "Kittlitz's", "Madagascar", "White-fronted", "Kentish (Maio)", "Kentish (Tuzla)")
all_ASRs_w_replace_1000_melt <- melt(all_ASRs_w_replace_1000, value.name = "ASR")
colnames(all_ASRs_w_replace_1000_melt) <- c("Population", "ASR")
ASR_w_repl_1000_JA_p_timexage_Phi_agexsex_CH_Phi_Quadraticxsex_p_yearxQuadratic <- replicate(1000, run_boot(Ceuta_F_A, Ceuta_C))

write.table(all_ASRs_w_replace_1000, file = "/home/luke/ceuta_ASR/bootstrap_products/1000_ASR_w_replace.txt", sep = "\t")

all_ASRs_w_replace_1000 <- as.data.frame(ASR_w_repl_1000_JA_p_timexage_Phi_agexsex_CH_Phi_Quadraticxsex_p_yearxQuadratic)
colnames(all_ASRs_w_replace_1000) <- "ASR"
all_ASRs_w_replace_1000 <- sort(all_ASRs_w_replace_1000$ASR)
all_ASRs_w_replace_1000[950]
all_ASRs_w_replace_1000[50]
all_ASRs_w_replace_1000 <- as.data.frame(all_ASRs_w_replace_1000)

all_ASRs_100 <- sort(all_ASRs_100$all_ASRs)
all_ASRs_100[95]
all_ASRs_100[5]

Bootstrap_histogram_all <- 
  ggplot(data = all_ASRs_w_replace_1000_melt, aes( x = ASR)) +
    geom_histogram(binwidth = 0.01) +
    facet_grid(Population ~ .) +
    theme_bw()

Bootstrap_histogram_all <- 
  ggplot(data = all_ASRs_w_replace_1000_melt, aes(x = ASR)) +
  geom_histogram(binwidth = 0.01) +
  facet_grid(Population ~ .) +
  theme_bw() +
  theme(text=element_text(family="Arial"),
        legend.position="none",
        legend.position = c(0, 1), 
        legend.justification = c(0, 1),
        legend.text=element_text(size=11),
        legend.title=element_blank(),
        legend.key.height=unit(0.8,"line"),
        legend.key.width=unit(0.8,"line"),
        legend.background = element_rect(fill=NA),
        axis.title.x = element_text(size=12, vjust=-0.1),
        axis.text.x  = element_text(size=10), 
        axis.title.y = element_text(size=12, vjust=-0.1),
        axis.text.y  = element_text(size=10),#, angle=90, hjust = 0.5, vjust = 0.5),
        #axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "grey")) +
  ylab("Frequency") +
  xlab("ASR estimate") #+
#   annotate("segment", x = 0.7597657, xend = 0.7597657, y = -Inf, yend = Inf, colour = "red") +
#   annotate("segment", x = 0.7989411, xend = 0.7989411, y = -Inf, yend = Inf, colour = "blue") +
#   annotate("segment", x = 0.5934889, xend = 0.5934889, y = -Inf, yend = Inf, colour = "blue")


ggsave(Bootstrap_histogram, 
       filename = "ASR_bootstrap_histogram_All_pops.jpg", 
       path = "/Users/Luke/Documents/Academic_Projects/PhD/Plover_Matrix_Modelling/Ceuta_Analysis/Figures",
       width = 4,
       height = 4, units = "in",
       dpi = 300)