#-----------------------------------------------------------------------------#
#               Matrix Modelling of ASR in Kittlitz's plover                  #
#                         Bootstrapping proceedure                            #
#                           (with replacement)                                #
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
setwd("/home/luke/comparative_ASR/") 

KiP_F_A <- 
  read.table("Data_files/Andava_KiP_juv_adult_capture_history_2009-2015_MARK_w_rings.txt",
             header=T,colClasses=c("character","factor","factor","factor","factor"))

KiP_C <- 
  read.table("Data_files/Andava_KiP_chick_capture_history_2009-2015_MARK_w_rings.txt",
             header=T,colClasses=c("character","factor","factor","factor","integer","numeric",
                                   "factor","factor","integer"))

plover_boot_KiP <- function(KiP_F_A, KiP_C) {
  
  KiP_C_boot <- KiP_C[sample(1:nrow(KiP_C), size = nrow(KiP_C), replace = TRUE), ]
  KiP_pres <- KiP_F_A$ring %in% KiP_C_boot$Ring
  KiP_F_A_boot1 <- KiP_F_A[KiP_pres, ]
  spare_KiP_F_A <- KiP_F_A[!KiP_pres, ]
  KiP_F_A_boot2 <- spare_KiP_F_A[sample(1:nrow(spare_KiP_F_A), size = nrow(KiP_F_A) - nrow(KiP_F_A_boot1), replace = TRUE), ]
  KiP_F_A_boot <- rbind(KiP_F_A_boot1, KiP_F_A_boot2)
  
  out <- list(KiP_C_boot = KiP_C_boot, KiP_F_A_boot = KiP_F_A_boot)
}

calc_ASR_KiP <- function(plover_boot_list, num_boot) {
  KiP_C <- plover_boot_list[["KiP_C_boot"]]
  KiP_F_A <- plover_boot_list[["KiP_F_A_boot"]]
  
  # remove ring column
  KiP_F_A <- KiP_F_A[,-3]
  KiP_C <- KiP_C[,-2]
  
  # Remove capture histories that have no resights (i.e., all zeros in the ch)
  KiP_C <- KiP_C[which(str_detect(KiP_C[,"ch"],"1") == TRUE),]
  
  # Create processed RMARK data format as CJS with 2 groups (sex and age 
  # initally ringed) and starting at year 2006
  KiP_F_A.proc=process.data(KiP_F_A,model="CJS",groups=c("sex","age"),
                            begin.time=2009,age.var=2,initial.age=c(1,0))

  # Create processed RMARK data format as CJS with 3 groups (sex, site, and 
  # cross-fostering treatment).
  KiP_C.proc=process.data(KiP_C,model="CJS",
                          groups=c("Sex","Year"))

  # Create the design data
  KiP_F_A.ddl=make.design.data(KiP_F_A.proc)
  KiP_C.ddl=make.design.data(KiP_C.proc)

  # adds firstyear/adult age field to design data in column "age"
  KiP_F_A.ddl=add.design.data(KiP_F_A.proc,KiP_F_A.ddl,"Phi","age",bins=c(0,1,6),right=F,name="age",replace=T)

  # create a dummy field called marked.as.adult which is 0 for the group initally ringed as juvenile and 1 for the group marked as adults.
  KiP_F_A.ddl$Phi$marked.as.adult=0
  KiP_F_A.ddl$Phi$marked.as.adult[KiP_F_A.ddl$Phi$initial.age.class=="A"]=1 
  KiP_F_A.ddl$p$marked.as.adult=0
  KiP_F_A.ddl$p$marked.as.adult[KiP_F_A.ddl$p$initial.age.class=="A"]=1

  # check parameter matrices to see if groups were binned correctly
  #PIMS(mark(Ceuta_F_A.proc,Ceuta_F_A.ddl,model.parameters=list(Phi=list(formula=~age+sex)),output=F),"Phi")
  
  # Create quadratic time variable so that it can be tested along side the annual models
  time <- c(0:(KiP_F_A.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  quad_time$time <- c(2009:2015)
  KiP_F_A.ddl$p=merge_design.covariates(KiP_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  KiP_F_A.ddl$Phi=merge_design.covariates(KiP_F_A.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  time <- c(0:(KiP_C.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  KiP_C.ddl$p=merge_design.covariates(KiP_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  KiP_C.ddl$Phi=merge_design.covariates(KiP_C.ddl$Phi,quad_time,bygroup=F,bytime=T)

  # Fledgling and adult survival models:
  # KiP:
  KiP_F_A_full_p_boot=function() 
  {
    setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/KiP/") 
    Phi.agexsex=list(formula=~age*sex) 
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
    model.list=mark.wrapper(cml,data=KiP_F_A.proc,ddl=KiP_F_A.ddl,
                                     threads = 4, brief = TRUE, delete = TRUE) #runs model list in MARK
    return(model.list) #stores completed model list
  }
  KiP_F_A_full_p_boot_run <- KiP_F_A_full_p_boot()
  KiP_AIC_table_F_A <- KiP_F_A_full_p_boot_run$model.table
  KiP_AIC_table_F_A$species <- "Kittlitz's"
  KiP_model_F_A_num <- as.numeric(rownames(KiP_F_A_full_p_boot_run$model.table[1,]))
  KiP_model_F_A_str <- KiP_F_A_full_p_boot_run$model.table[1,]

  # Chick survival model:
  # KiP:
  KiP_C_full_p_boot=function() 
  {
    setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/KiP/") 
    Phi.Time.x.sex=list(formula=~Sex*Time) #Phi(sex(T))
    Phi.Quadratic.x.sex=list(formula=~Sex*Quadratic) #Phi(sex(Quadratic))
    Phi.Cubic.x.sex=list(formula=~Sex*Cubic) #Phi(sex(Quadratic))
    Phi.Sex=list(formula=~Sex)
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
    model.list=mark.wrapper(cml, data=KiP_C.proc, ddl=KiP_C.ddl,
                                     threads = 4, brief = TRUE, delete = TRUE) #runs model list in MARK
    return(model.list) #stores completed model list
  }
  KiP_C_full_p_boot_run <- KiP_C_full_p_boot()
  KiP_AIC_table_C <- KiP_C_full_p_boot_run$model.table
  KiP_AIC_table_C$species <- "Kittlitz's"
  KiP_model_C_num <- as.numeric(rownames(KiP_C_full_p_boot_run$model.table[1,]))
  KiP_model_C_str <- KiP_C_full_p_boot_run$model.table[1,]
  
  # Extract reals
  # extract and format survival rates from chick model output
  # KiP:
  KiP_C_reals <- KiP_C_full_p_boot_run[[KiP_model_C_num]]$results$real
  Groups <- data.frame(str_split_fixed(rownames(KiP_C_reals), " ", n = 5))
  KiP_C_reals <- cbind(Groups, KiP_C_reals)
  KiP_C_reals <- KiP_C_reals[which(KiP_C_reals$X1 == "Phi"),]
  KiP_C_reals$Sex <- unlist(str_extract_all(KiP_C_reals$X2,"[FM]"))
  KiP_C_reals$Sex <- as.factor(ifelse(KiP_C_reals$Sex == "F","Female","Male"))
  if(nrow(KiP_C_reals) == 2)
  {
    KiP_Survival_to_Fledge_F <- 
      KiP_C_reals[which(KiP_C_reals$Sex == "Female"),
                  c("estimate")]^25
    KiP_Survival_to_Fledge_M <- 
      KiP_C_reals[which(KiP_C_reals$Sex == "Male"),
                  c("estimate")]^25
  }
  if(nrow(KiP_C_reals) != 2){
    KiP_Survival_to_Fledge_F <- 
      prod(KiP_C_reals[which(KiP_C_reals$Sex == "Female"),
                       c("estimate")][c(1:26)])
    KiP_Survival_to_Fledge_M <- 
      prod(KiP_C_reals[which(KiP_C_reals$Sex == "Male"),
                       c("estimate")][c(1:26)])
  }
  estimate <- c(KiP_Survival_to_Fledge_F, KiP_Survival_to_Fledge_M)
  Sex <- c("Female", "Male")
  age <- c("Chick", "Chick")
  Sex_Age <- paste(Sex,age,sep = "_")
  KiP_sex_chick_survival <- data.frame(Sex_Age, estimate)
  KiP_sex_chick_survival$species <- "Kittlitz's"
  
  # extract and format survival rates from fledgling and adult model output
  KiP_F_A_reals <- KiP_F_A_full_p_boot_run[[KiP_model_F_A_num]]$results$real
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

  ###### ASR estimation #######################################################
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
  
  KiP_boot_list <- 
    list(KiP_AIC_table_C, 
         KiP_AIC_table_F_A, 
         KiP_Survival_rates,
         KiP_ASR)
  
}

# bootstrapping
run_boot_KiP <- function(num_boot, KiP_F_A, KiP_C){
  plover_list_KiP <- plover_boot_KiP(KiP_F_A, KiP_C)
  Result <- calc_ASR_KiP(plover_list_KiP, num_boot)
}

niter <- 1000

KiP_boot_1000_full_p <- sapply(1:niter, run_boot_KiP, KiP_F_A, KiP_C)

# extract AIC_table_C
KiP_AIC_table_C <- do.call(rbind, lapply(seq(from = 1, to = niter * 4, by = 4), function(x) KiP_boot_1000_full_p[[x]]))
num_mods <- nrow(KiP_AIC_table_C)/niter
KiP_AIC_table_C$iter <- rep(1:niter, each = num_mods)

# extract AIC_table_F_A
KiP_AIC_table_F_A <- do.call(rbind, lapply(seq(from = 2, to = niter * 4, by = 4), function(x) KiP_boot_1000_full_p[[x]]))
num_mods <- nrow(KiP_AIC_table_F_A)/niter
KiP_AIC_table_F_A$iter <- rep(1:niter, each = num_mods)

# extract KiP_Survival_rates
KiP_Survival_rates <- do.call(rbind, lapply(seq(from = 3, to = niter * 4, by = 4), function(x) KiP_boot_1000_full_p[[x]]))
KiP_Survival_rates$iter <- rep(1:niter, each = 6)

# extract KiP_ASR
KiP_ASR <- sapply(seq(from = 4, to = niter * 4, by = 4), function(x) KiP_boot_1000_full_p[[x]])
KiP_ASR <- data.frame(ASR = unname(KiP_ASR), iter = 1:niter)

write.table(KiP_AIC_table_C, 
            file = "/home/luke/comparative_ASR/Bootstrap/Total_output/KiP_AIC_table_C.txt",
            sep = "\t", row.names = FALSE, 
            col.names = TRUE)

write.table(KiP_AIC_table_F_A, 
            file = "/home/luke/comparative_ASR/Bootstrap/Total_output/KiP_AIC_table_F_A.txt",
            sep = "\t", row.names = FALSE, 
            col.names = TRUE)

write.table(KiP_Survival_rates, 
            file = "/home/luke/comparative_ASR/Bootstrap/Total_output/KiP_Survival_rates.txt",
            sep = "\t", row.names = FALSE, 
            col.names = TRUE)

write.table(KiP_ASR, 
            file = "/home/luke/comparative_ASR/Bootstrap/Total_output/KiP_ASR.txt",
            sep = "\t", row.names = FALSE, 
            col.names = TRUE)