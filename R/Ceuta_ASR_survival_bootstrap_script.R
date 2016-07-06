#-----------------------------------------------------------------------------#
#              Matrix Modelling of ASR in Snowy plover                        #
#                         Bootstrapping proceedure                            #
#                           (with replacement)                                #
#-----------------------------------------------------------------------------#

# Reference RMark functions
library(RMark) 
library(stringr)
library(dplyr)
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
# HSR ---> the hatching sex ratio (default is 0.5)
# iterations ---> the number of iterations to simulate (default is 20)
freq_dep_SSD_SR <- 
  function (A, n = rep(10, nrow(A)), h = 1, k = 1, iterations = 30, HSR = 0.5) 
  {
    x <- length(n) # Number of stages in matrix
    t <- iterations # Number of time steps to simulate
    stage <- matrix(numeric(x * t), nrow = x) # an empty t by x matrix
    for (i in 1:t) { # for loop that goes through each of t time steps
      stage[,i] <- n # stage distribution at time t
      M2 <- stage[4, i] # number of male adults at time t
      F2 <- stage[2, i] # number of female adults at time t
      A[1,x/2]        <- (k*M2)/(M2+(F2/h))*HSR # F freq-dep fecundity of F
      A[(x/4)*3,x/2]  <- (k*M2)/(M2+(F2/h))*HSR # F freq-dep fecundity of M
      A[1,x]          <- (k*F2)/(M2+(F2/h))*HSR # M freq-dep fecundity of F
      A[(x/4)*3,x]    <- (k*F2)/(M2+(F2/h))*HSR # M freq-dep fecundity of M
      n <- A %*% n # define the new n (i.e., new stage distribution at time t)
    }
    rownames(stage) <- rownames(A) # define rownames of stage matrix
    colnames(stage) <- 0:(t - 1) # define colnames of stage matrix
    w <- stage[, t] # define stable stage as the last stage
    stable.stage <- w/sum(w)
    ASR <- stable.stage[x]/(stable.stage[x/2]+stable.stage[x])
    FYSR <- stable.stage[x-1]/(stable.stage[(x/2)-1]+stable.stage[x-1])
    pop.proj <- list(ASR = ASR, # make a list of results
                     FYSR = FYSR,
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

Ceuta_F_A <- 
  read.table("Data_files/Ceuta_juv_adult_capture_history_MARK_w_rings.txt",
             header=T,colClasses=c("factor","character","factor","factor"))

Ceuta_C <- 
  read.table("Data_files/Ceuta_chick_capture_history_2006-2012_MARK_w_rings.txt",
             header=T,colClasses=c("factor","character","factor","integer",
                                   "numeric","factor","factor","factor",
                                   "factor","integer","factor"))

names(Ceuta_C)[1] <- "Ring"

plover_boot_Ceuta <- function(Ceuta_F_A, Ceuta_C) {
  
  Ceuta_C_boot <- Ceuta_C[sample(1:nrow(Ceuta_C), size = nrow(Ceuta_C), replace = TRUE), ]
  Ceuta_pres <- Ceuta_F_A$ring %in% Ceuta_C_boot$Ring
  Ceuta_F_A_boot1 <- Ceuta_F_A[Ceuta_pres, ]
  spare_Ceuta_F_A <- Ceuta_F_A[!Ceuta_pres, ]
  Ceuta_F_A_boot2 <- spare_Ceuta_F_A[sample(1:nrow(spare_Ceuta_F_A), size = nrow(Ceuta_F_A) - nrow(Ceuta_F_A_boot1), replace = TRUE), ]
  Ceuta_F_A_boot <- rbind(Ceuta_F_A_boot1, Ceuta_F_A_boot2)
  
  out <- list(Ceuta_C_boot = Ceuta_C_boot, Ceuta_F_A_boot = Ceuta_F_A_boot)
}

calc_ASR_Ceuta <- function(plover_boot_list, num_boot) {
  Ceuta_C <- plover_boot_list[["Ceuta_C_boot"]]
  Ceuta_F_A <- plover_boot_list[["Ceuta_F_A_boot"]]
  
  # remove ring column
  Ceuta_F_A <- Ceuta_F_A[,-1]
  Ceuta_C <- Ceuta_C[,-1]
  
  # Remove capture histories that have no resights (i.e., all zeros in the ch)
  Ceuta_C <- Ceuta_C[which(str_detect(Ceuta_C[,"ch"],"1") == TRUE),]
  
  # Create processed RMARK data format as CJS with 2 groups (sex and age 
  # initally ringed) and starting at year 2006
  Ceuta_F_A.proc=process.data(Ceuta_F_A,model="CJS",groups=c("sex","age"),
                              begin.time=2006,age.var=2,initial.age=c(1,0))
  
  # Create processed RMARK data format as CJS with 3 groups (sex, site, and 
  # cross-fostering treatment).
  Ceuta_C.proc=process.data(Ceuta_C,model="CJS",
                            groups=c("Sex","Care_site","CF","Year"))
  
  # Create the design data
  Ceuta_F_A.ddl=make.design.data(Ceuta_F_A.proc)
  Ceuta_C.ddl=make.design.data(Ceuta_C.proc)
  
  # adds firstyear/adult age field to design data in column "age"
  Ceuta_F_A.ddl=add.design.data(Ceuta_F_A.proc,Ceuta_F_A.ddl,"Phi","age",bins=c(0,1,7),right=F,name="age",replace=T)
  
  # create a dummy field called marked.as.adult which is 0 for the group initally ringed as juvenile and 1 for the group marked as adults.
  Ceuta_F_A.ddl$Phi$marked.as.adult=0
  Ceuta_F_A.ddl$Phi$marked.as.adult[Ceuta_F_A.ddl$Phi$initial.age.class=="A"]=1 
  Ceuta_F_A.ddl$p$marked.as.adult=0
  Ceuta_F_A.ddl$p$marked.as.adult[Ceuta_F_A.ddl$p$initial.age.class=="A"]=1
  
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
  
  time <- c(0:(Ceuta_C.proc$nocc[1]-1))
  Quadratic <- time^2
  Cubic <- time^3
  quad_time <- data.frame(time, Quadratic, Cubic)
  Ceuta_C.ddl$p=merge_design.covariates(Ceuta_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  Ceuta_C.ddl$Phi=merge_design.covariates(Ceuta_C.ddl$Phi,quad_time,bygroup=F,bytime=T)
  
  # Fledgling and adult survival models:
  # Ceuta:
  Ceuta_F_A_full_p_boot = function() 
  {
    setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/Ceuta/") # set wd so that results go to the correct folder
    Phi.agexsex=list(formula = ~ age * sex)
    p.sex = list(formula =  ~ sex) #p(sex(.))
    p.age = list(formula =  ~ age) #p(age(.))
    p.dot = list(formula =  ~ 1) #p(.)
    p.time = list(formula =  ~ time) #p(t)
    p.Time = list(formula =  ~ Time)
    p.Quadratic = list(formula =  ~ Quadratic)
    p.sexxtime = list(formula =  ~ sex * time)
    p.agextime = list(formula =  ~ age * time)
    p.sexxTime = list(formula =  ~ sex * Time)
    p.agexTime = list(formula =  ~ age * Time)
    p.agexsex = list(formula =  ~ age * sex)
    p.Quadraticxsex = list(formula =  ~ Quadratic * sex)
    p.Quadraticxage = list(formula =  ~ Quadratic * age)
    #p.Quadraticxagexsex = list(formula =  ~ Quadratic * age * sex)
    #p.Timexagexsex = list(formula =  ~ Time * age * sex)
    #p.timexagexsex = list(formula =  ~ time * age * sex)
    p.Quadratic_agexsex = list(formula =  ~ Quadratic + age * sex)
    p.Time_agexsex = list(formula =  ~ Time + age * sex)
    p.time_agexsex = list(formula =  ~ time + age * sex)
    p.sex_time = list(formula =  ~ sex + time)
    p.age_time = list(formula =  ~ age + time)
    p.sex_Time = list(formula =  ~ sex + Time)
    p.age_Time = list(formula =  ~ age + Time)
    p.age_sex = list(formula =  ~ age + sex)
    p.Quadratic_sex = list(formula =  ~ Quadratic + sex)
    p.Quadratic_age = list(formula =  ~ Quadratic + age)
    p.Quadratic_age_sex = list(formula =  ~ Quadratic + age + sex)
    p.Time_age_sex = list(formula =  ~ Time + age + sex)
    p.time_age_sex = list(formula =  ~ time + age + sex)
    cml = create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
    model.list = mark.wrapper(cml, data = Ceuta_F_A.proc, ddl = Ceuta_F_A.ddl,
                              threads = 4, brief = TRUE, delete = TRUE) #runs model list in MARK
    return(model.list) #stores completed model list
  }
  Ceuta_F_A_full_p_boot_run <- Ceuta_F_A_full_p_boot()
  Ceuta_AIC_table_F_A <- Ceuta_F_A_full_p_boot_run$model.table
  Ceuta_AIC_table_F_A$species <- "Snowy"
  Ceuta_model_F_A_num <- as.numeric(rownames(Ceuta_F_A_full_p_boot_run$model.table[1,]))

  # Chick survival model:
  # Ceuta:
  Ceuta_C_full_p_boot=function() 
  {
    setwd("/home/luke/comparative_ASR/Bootstrap/RMark_output_files/Ceuta/") 
    #Phi.Time.x.sex=list(formula=~Sex*Time) #Phi(sex(T))
    Phi.Quadratic.x.sex=list(formula=~Sex*Quadratic) #Phi(sex(Quadratic))
    #Phi.Cubic.x.sex=list(formula=~Sex*Cubic) #Phi(sex(Quadratic))
    #Phi.Sex=list(formula=~Sex)
    p.dot=list(formula=~1) #p(.)
    p.Time=list(formula=~Time)
    p.Quadratic=list(formula=~Quadratic)
    p.Cubic=list(formula=~Cubic)
    p.year=list(formula=~Year)
    p.sex=list(formula=~Sex)
    p.year.x.Time=list(formula=~Year*Time)
    p.year.x.Quadratic=list(formula=~Year*Quadratic)
    p.year.x.Time=list(formula=~Sex*Time)
    p.sex.x.Quadratic=list(formula=~Sex*Quadratic)
    # p.sex.x.Time.x.Sex=list(formula=~Year*Time*Sex)
    # p.year.x.Quadratic.x.Sex=list(formula=~Year*Quadratic*Sex)
    p.year_Time.x.Sex=list(formula=~Year+Time*Sex)
    p.year_Quadratic.x.Sex=list(formula=~Year+Quadratic*Sex)
    p.sex.Time=list(formula=~Sex+Time)
    p.sex.Quadratic=list(formula=~Sex+Quadratic)
    p.year.Time=list(formula=~Year+Time)
    p.year.Quadratic=list(formula=~Year+Quadratic)
    p.year.Time.Sex=list(formula=~Year+Time+Sex)
    p.year.Quadratic.Sex=list(formula=~Year+Quadratic+Sex)
    cml=create.model.list("CJS") #creates model list for all a priori models above that begin with Phi. or p.
    model.list=mark.wrapper(cml, data=Ceuta_C.proc, ddl=Ceuta_C.ddl,
                            threads = 4, brief = TRUE, delete = TRUE) #runs model list in MARK
    return(model.list) #stores completed model list
  }
  Ceuta_C_full_p_boot_run <- Ceuta_C_full_p_boot()
  Ceuta_AIC_table_C <- Ceuta_C_full_p_boot_run$model.table
  Ceuta_AIC_table_C$species <- "Snowy"
  Ceuta_model_C_num <- as.numeric(rownames(Ceuta_C_full_p_boot_run$model.table[1,]))

  # Extract reals
  # extract and format survival rates from chick model output
  # Ceuta:
  Ceuta_C_reals <- Ceuta_C_full_p_boot_run[[Ceuta_model_C_num]]$results$real
  Groups <- data.frame(str_split_fixed(rownames(Ceuta_C_reals), " ", n = 5))
  Ceuta_C_reals <- cbind(Groups, Ceuta_C_reals)
  Ceuta_C_reals <- Ceuta_C_reals[which(Ceuta_C_reals$X1 == "Phi"),]
  Ceuta_C_reals$Sex <- unlist(str_extract_all(Ceuta_C_reals$X2,"[FM]"))
  Ceuta_C_reals$Sex <- as.factor(ifelse(Ceuta_C_reals$Sex == "F","Female","Male"))
  if(nrow(Ceuta_C_reals) == 2)
  {
    Ceuta_Survival_to_Fledge_F <- 
      Ceuta_C_reals[which(Ceuta_C_reals$Sex == "Female"),
                    c("estimate")]^25
    Ceuta_Survival_to_Fledge_M <- 
      Ceuta_C_reals[which(Ceuta_C_reals$Sex == "Male"),
                    c("estimate")]^25
  }
  if(nrow(Ceuta_C_reals) != 2){
    Ceuta_Survival_to_Fledge_F <- 
      prod(Ceuta_C_reals[which(Ceuta_C_reals$Sex == "Female"),
                         c("estimate")][c(1:26)])
    Ceuta_Survival_to_Fledge_M <- 
      prod(Ceuta_C_reals[which(Ceuta_C_reals$Sex == "Male"),
                         c("estimate")][c(1:26)])
  }
  estimate <- c(Ceuta_Survival_to_Fledge_F, Ceuta_Survival_to_Fledge_M)
  Sex <- c("Female", "Male")
  age <- c("Chick", "Chick")
  Sex_Age <- paste(Sex,age,sep = "_")
  Ceuta_sex_chick_survival <- data.frame(Sex_Age, estimate)
  Ceuta_sex_chick_survival$species <- "Snowy"
  
  # extract and format survival rates from fledgling and adult model output
  Ceuta_F_A_reals <- Ceuta_F_A_full_p_boot_run[[Ceuta_model_F_A_num]]$results$real
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
  
  # Extract reals
  # extract and format capture rates from chick model output
  Ceuta_C_reals_p <- Ceuta_C_full_p_boot_run[[Ceuta_model_C_num]]$results$real
  Groups <- data.frame(str_split_fixed(rownames(Ceuta_C_reals_p), " ", n = 5))
  Ceuta_C_reals_p <- cbind(Groups, Ceuta_C_reals_p)
  Ceuta_C_reals_p <- Ceuta_C_reals_p[which(Ceuta_C_reals_p$X1 == "p"),]
  Ceuta_C_reals_p$Sex <- unlist(str_extract_all(Ceuta_C_reals_p$X2,"[FM]"))
  Ceuta_C_reals_p$Sex <- as.factor(ifelse(Ceuta_C_reals_p$Sex == "F","Female","Male"))
  if(nrow(Ceuta_C_reals_p) == 2)
  {
    Ceuta_Survival_to_Fledge_F_p <- 
      Ceuta_C_reals_p[which(Ceuta_C_reals_p$Sex == "Female"),
                    c("estimate")]
    Ceuta_Survival_to_Fledge_M_p <- 
      Ceuta_C_reals_p[which(Ceuta_C_reals_p$Sex == "Male"),
                    c("estimate")]
  }
  if(nrow(Ceuta_C_reals_p) != 2){
    Ceuta_Survival_to_Fledge_F_p <- 
      mean(Ceuta_C_reals_p[which(Ceuta_C_reals_p$Sex == "Female"),
                         c("estimate")][c(1:26)])
    Ceuta_Survival_to_Fledge_M_p <- 
      mean(Ceuta_C_reals_p[which(Ceuta_C_reals_p$Sex == "Male"),
                         c("estimate")][c(1:26)])
  }
  estimate <- c(Ceuta_Survival_to_Fledge_F_p, Ceuta_Survival_to_Fledge_M_p)
  Sex <- c("Female", "Male")
  age <- c("Chick", "Chick")
  Sex_Age <- paste(Sex,age,sep = "_")
  Ceuta_sex_chick_capture <- data.frame(Sex_Age, estimate)
  Ceuta_sex_chick_capture$species <- "Snowy"
  
  # extract and format capture rates from fledgling and adult model output
  Ceuta_F_A_reals_p <- Ceuta_F_A_full_p_boot_run[[Ceuta_model_F_A_num]]$results$real
  Groups <- data.frame(str_split_fixed(rownames(Ceuta_F_A_reals_p), " ", n = 5))
  Ceuta_F_A_reals_p <- cbind(Groups, Ceuta_F_A_reals_p)
  Ceuta_F_A_reals_p <- Ceuta_F_A_reals_p[which(Ceuta_F_A_reals_p$X1 == "p"),]
  Ceuta_F_A_reals_p$age <- unlist(str_extract_all(Ceuta_F_A_reals_p$X2,"[AJ]"))
  Ceuta_F_A_reals_p$age <- as.factor(ifelse(Ceuta_F_A_reals_p$age == "A","Adult","Juvenile"))
  Ceuta_F_A_reals_p$Sex <- unlist(str_extract_all(Ceuta_F_A_reals_p$X2,"[FM]"))
  Ceuta_F_A_reals_p$Sex <- as.factor(ifelse(Ceuta_F_A_reals_p$Sex == "F","Female","Male"))
  Ceuta_F_A_reals_p$Sex_Age <- paste(Ceuta_F_A_reals_p$Sex,Ceuta_F_A_reals_p$age,sep = "_")
  
  Ceuta_Survival_A_F_p <- 
    mean(Ceuta_F_A_reals_p[which(Ceuta_F_A_reals_p$Sex_Age == "Female_Adult"),
                         c("estimate")])
  Ceuta_Survival_A_M_p <- 
    mean(Ceuta_F_A_reals_p[which(Ceuta_F_A_reals_p$Sex_Age == "Male_Adult"),
                         c("estimate")])
  
  Ceuta_Survival_F_F_p <- 
    mean(Ceuta_F_A_reals_p[which(Ceuta_F_A_reals_p$Sex_Age == "Female_Juvenile"),
                         c("estimate")])
  Ceuta_Survival_F_M_p <- 
    mean(Ceuta_F_A_reals_p[which(Ceuta_F_A_reals_p$Sex_Age == "Male_Juvenile"),
                         c("estimate")])
  
  estimate <- c(Ceuta_Survival_A_F_p, Ceuta_Survival_A_M_p, Ceuta_Survival_F_F_p, Ceuta_Survival_F_M_p)
  Sex_Age <- c("Female_Adult", "Male_Adult", "Female_Juvenile", "Male_Juvenile")
  Ceuta_F_A_sex_capture <- data.frame(Sex_Age, estimate)
  Ceuta_F_A_sex_capture$species <- "Snowy"
  Ceuta_Capture_rates <- rbind(Ceuta_F_A_sex_capture, Ceuta_sex_chick_capture)
  row.names(Ceuta_Capture_rates) <- NULL

  ###### ASR estimation #######################################################
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
                   # Define HSR sex ratio
                   HSR = 0.486)
  
  # Ceuta matrix:
  Ceuta_matrix <- plover_matrix(Ceuta_VR)
  
  # Determine the ASR at the stable stage distribution (assume h is 1, which 
  # is monogamy)
  Ceuta_SR_h_1 <- freq_dep_SSD_SR(A = Ceuta_matrix, h = 1, k = 3)
  
  # Extract ASR
  Ceuta_ASR <- Ceuta_SR_h_1$ASR
  
  # Extract FYSR
  Ceuta_FYSR <- Ceuta_SR_h_1$FYSR
  
  Ceuta_boot_list <- 
    list(Ceuta_AIC_table_C, 
        Ceuta_AIC_table_F_A, 
        Ceuta_Survival_rates,
        Ceuta_Capture_rates,
        Ceuta_ASR,
        Ceuta_FYSR)
    
}

# bootstrapping
run_boot_Ceuta <- function(num_boot, Ceuta_F_A, Ceuta_C){
  plover_list_Ceuta <- plover_boot_Ceuta(Ceuta_F_A, Ceuta_C)
  Result <- calc_ASR_Ceuta(plover_list_Ceuta, num_boot)
}

niter <- 1000

Ceuta_boot_1000_full_p <- sapply(1:niter, run_boot_Ceuta, Ceuta_F_A, Ceuta_C)

# extract AIC_table_C
Ceuta_AIC_table_C <- do.call(rbind, lapply(seq(from = 1, to = niter * 6, by = 6), function(x) Ceuta_boot_1000_full_p[[x]]))
num_mods <- nrow(Ceuta_AIC_table_C)/niter
Ceuta_AIC_table_C$iter <- rep(1:niter, each = num_mods)

# extract AIC_table_F_A
Ceuta_AIC_table_F_A <- do.call(rbind, lapply(seq(from = 2, to = niter * 6, by = 6), function(x) Ceuta_boot_1000_full_p[[x]]))
num_mods <- nrow(Ceuta_AIC_table_F_A)/niter
Ceuta_AIC_table_F_A$iter <- rep(1:niter, each = num_mods)

# extract Ceuta_Survival_rates
Ceuta_Survival_rates <- do.call(rbind, lapply(seq(from = 3, to = niter * 6, by = 6), function(x) Ceuta_boot_1000_full_p[[x]]))
Ceuta_Survival_rates$iter <- rep(1:niter, each = 6)


# extract Ceuta_Capture_rates
Ceuta_Capture_rates <- do.call(rbind, lapply(seq(from = 4, to = niter * 6, by = 6), function(x) Ceuta_boot_1000_full_p[[x]]))
Ceuta_Capture_rates$iter <- rep(1:niter, each = 6)

# extract Ceuta_ASR
Ceuta_ASR <- sapply(seq(from = 5, to = niter * 6, by = 6), function(x) Ceuta_boot_1000_full_p[[x]])
Ceuta_ASR <- data.frame(ASR = unname(Ceuta_ASR), iter = 1:niter)

# extract Ceuta_FYSR
Ceuta_FYSR <- sapply(seq(from = 6, to = niter * 6, by = 6), function(x) Ceuta_boot_1000_full_p[[x]])
Ceuta_FYSR <- data.frame(ASR = unname(Ceuta_FYSR), iter = 1:niter)

write.table(Ceuta_AIC_table_C, 
            file = "/home/luke/comparative_ASR/Bootstrap/Total_output/June/Ceuta_AIC_table_C_June.txt",
            sep = "\t", row.names = FALSE, 
            col.names = TRUE, append = FALSE)

write.table(Ceuta_AIC_table_F_A, 
            file = "/home/luke/comparative_ASR/Bootstrap/Total_output/June/Ceuta_AIC_table_F_A_June.txt",
            sep = "\t", row.names = FALSE, 
            col.names = TRUE, append = FALSE)

write.table(Ceuta_Survival_rates, 
            file = "/home/luke/comparative_ASR/Bootstrap/Total_output/June/Ceuta_Survival_rates_June.txt",
            sep = "\t", row.names = FALSE, 
            col.names = TRUE, append = FALSE)

write.table(Ceuta_Capture_rates, 
            file = "/home/luke/comparative_ASR/Bootstrap/Total_output/June/Ceuta_Capture_rates_June.txt",
            sep = "\t", row.names = FALSE, 
            col.names = TRUE, append = FALSE)

write.table(Ceuta_ASR, 
            file = "/home/luke/comparative_ASR/Bootstrap/Total_output/June/Ceuta_ASR_June.txt",
            sep = "\t", row.names = FALSE, 
            col.names = TRUE, append = FALSE)

write.table(Ceuta_FYSR, 
            file = "/home/luke/comparative_ASR/Bootstrap/Total_output/June/Ceuta_FYSR_June.txt",
            sep = "\t", row.names = FALSE, 
            col.names = TRUE, append = FALSE)