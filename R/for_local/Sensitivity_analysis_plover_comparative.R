#-----------------------------------------------------------------------------#
#         Deterministic Projection Matrix Model for Snowy Plover ASR with     #
#                     frequency dependent fecundity function                  #
#                                                                             #
#                           Luke Eberhart-Phillips                            #
#                               June 29, 2016                                 #
#-----------------------------------------------------------------------------#

# The following packages are needed to run and visualize analyses
library(reshape2)
library(dplyr)
library(extrafont)
library(ggplot2)
library(gridExtra)
library(popdemo)
library(popbio)
library(RColorBrewer)
library(grid)
library(stringr)

# Find fonts from computer that are candara or Candara
font_import(pattern="[M/m]enlo", prompt = FALSE) 

###############################################################################
# The following functions are needed to run analyses. Note: Some are dependent 
# on others and need to be loaded first.

# NOTE: many of the following functions were modified from Chris Stubben's 
# "popbio" package which were all written for one-sex matrix models estimating
# lambda.  I've modified them to accomodate two-sex matrix models with ASR as 
# the response instead of lambda.

# Function to compile the two-sex projection matrix from a list of vital rates
# "VR_list"
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

freq_dep_SSD_ASR <- 
  function (A, n = rep(10, nrow(A)), h = 1, k = 1, iterations = 30, HSR = 0.5) 
    {
    x <- length(n) # Number of transitional stages in matrix
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
    pop.proj <- list(ASR = ASR, # make a list of results
                     stable.stage = stable.stage, 
                     stage.vectors = stage,
                     SSD_M2 = stable.stage[4],
                     SSD_F2 = stable.stage[2])
    pop.proj # print the list as output to the function
  }

# Functions to run a sensitivity analysis on all model components
lower_level_sens_analysis <- 
  function(freq_dep_ASR, VR_list, chick_surv = TRUE)
  {
    # check to see if valid estimates of sex=specific chick survival are available
    if(!chick_surv)
    {
      # make a list of all parameters
      vr <- list(
        F_Juv_survl = VR_list$F_Juv_survl,
        F_Adt_survl = VR_list$F_Adt_survl,
        M_Juv_survl = VR_list$M_Juv_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        h = VR_list$h,
        k = VR_list$k,
        PSR = VR_list$PSR,
        M2 = unname(freq_dep_ASR$SSD_M2),
        F2 = unname(freq_dep_ASR$SSD_F2))
      # make a matrix of the elements
      el <- expression(0, ((k*M2)/(M2+(F2/h)))*(1-PSR), 0, ((k*F2)/(M2+(F2/h)))*(1-PSR),
                       F_Juv_survl, F_Adt_survl, 0, 0,
                       0, ((k*M2)/(M2+(F2/h)))*PSR, 0, ((k*F2)/(M2+(F2/h)))*PSR,
                       0, 0, M_Juv_survl, M_Adt_survl)
    }
    else
    {
      # make a list of all parameters
      vr <- list(
        F_Chk_survl = VR_list$F_Chk_survl,
        F_Juv_survl = VR_list$F_Juv_survl,
        F_Adt_survl = VR_list$F_Adt_survl,
        M_Chk_survl = VR_list$M_Chk_survl,
        M_Juv_survl = VR_list$M_Juv_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        h = VR_list$h,
        k = VR_list$k,
        PSR = VR_list$PSR,
        M2 = unname(freq_dep_ASR$SSD_M2),
        F2 = unname(freq_dep_ASR$SSD_F2))
      # make a matrix of the elements
      el <- expression(0, ((k*M2)/(M2+(F2/h)))*(1-PSR), 0, ((k*F2)/(M2+(F2/h)))*(1-PSR),
                       (F_Chk_survl*F_Juv_survl), F_Adt_survl, 0, 0,
                       0, ((k*M2)/(M2+(F2/h)))*PSR, 0, ((k*F2)/(M2+(F2/h)))*PSR,
                       0, 0, (M_Chk_survl*M_Juv_survl), M_Adt_survl)
    }
    # calculate the effect of proportional changes in vital rates
    n <- length(vr)
    vr_nums <- seq(0,2,.1)
    vrsen <- matrix(numeric(n*length(vr_nums)), ncol=n, dimnames=list(vr_nums, names(vr)))
    for (h in 1:n)
    {
      vr2 <- vr
      for (i in 1:length(vr_nums))
      {
        vr2[[h]]<-vr_nums[i]
        A<-matrix(sapply(el, eval, vr2, NULL), nrow=sqrt(length(el)), byrow=TRUE)
        vrsen[i,h] <- eigen(A)$vectors[2,1]/(eigen(A)$vectors[2,1]+eigen(A)$vectors[4,1])
      }
    }
    vrelas <- matrix(numeric(n*length(vr_nums)), ncol=n, dimnames=list(vr_nums, names(vr)))
    for (h in 1:n)
    {
      for (i in 1:length(vr_nums))
      {
        vr2 <- vr
        vr2[[h]]<-vr_nums[i]*vr2[[h]]
        A<-matrix(sapply(el, eval, vr2 , NULL), nrow=sqrt(length(el)), byrow=TRUE)
        vrelas[i,h] <- (eigen(A)$vectors[2,1]/eigen(A)$vectors[4,1])/unname(freq_dep_ASR_f_m$ASR_f_m)
      }
    }
    if(!chick_surv)
    {
      colnames(vrsen) <- c("Female fledgling survival", 
                           "Female adult survival", 
                           "Male fledgling survival", 
                           "Male adult survival",
                           "Polygamy index (h)", 
                           "Clutch size",
                           "Primary sex ratio",
                           "Breeding males",
                           "Breeding females")
      colnames(vrelas) <- c("Female fledgling survival", 
                            "Female adult survival", 
                            "Male fledgling survival", 
                            "Male adult survival",
                            "Polygamy index (h)", 
                            "Clutch size",
                            "Primary sex ratio",
                            "Breeding males",
                            "Breeding females")
    }
    else
    {
      colnames(vrsen) <- c("Female chick survival", 
                           "Female fledgling survival", 
                           "Female adult survival", 
                           "Male chick survival", 
                           "Male fledgling survival", 
                           "Male adult survival",
                           "Polygamy index (h)", 
                           "Clutch size",
                           "Primary sex ratio",
                           "Breeding males",
                           "Breeding females")
      colnames(vrelas) <- c("Female chick survival", 
                            "Female fledgling survival", 
                            "Female adult survival", 
                            "Male chick survival", 
                            "Male fledgling survival", 
                            "Male adult survival",
                            "Polygamy index (h)", 
                            "Clutch size",
                            "Primary sex ratio",
                            "Breeding males",
                            "Breeding females")
    }
    Sensitivities <- melt(vrsen)
    colnames(Sensitivities) <- c("Perturbation", "Vitalrate", "Sensitivity")
    Elasticities <- melt(vrelas)
    colnames(Elasticities) <- c("Perturbation", "Vitalrate", "Elasticity")
    results <- list(Sensitivities = Sensitivities,
                    Elasticities = Elasticities,
                    Element_expression = el)
    results
  }

ASR_analysis <- 
  function (A, zero = TRUE) 
    {
  ev <- eigen(A) # makes list of the eigen values and eigen vectors of A
  lmax <- which.max(Re(ev$values)) # index of dominant eigen value
  W <- ev$vectors # Eigen vectors
  w <- abs(Re(W[, lmax])) # dominant eigen vector
  stable.stage = w/sum(w) # stable stage distribution
  ASR <- stable.stage[4]/(stable.stage[2]+stable.stage[4]) # SSD ASR
  V <- try(Conj(solve(W)), silent = TRUE) # check if possible to proceed
  if (class(V) == "try-error") {
    ASR.analysis <- list(ASR = ASR, stable.stage = stable.stage, 
                         sensitivities = A * NA, elasticities = A * NA)
  }
  else {
    v <- abs(Re(V[lmax, ])) # solve matrix to find ...(?)
    s <- v %o% w # outer product of v and w
    if (zero) {
      s[A == 0] <- 0
    }
    e <- s * A/ASR # calculate elasticities
    x <- dimnames(A) # get vital rate names
    dimnames(s) <- x # assign vital rate names to s
    names(w) <- x[[1]]
    names(v) <- x[[1]]
    ASR.analysis <- list(ASR = ASR, stable.stage = stable.stage, 
                         sensitivities = s, elasticities = e)
  }
  ASR.analysis
}

vitalsens_ASR_f_m <- 
  function (elements, VR_list, freq_dep_ASR_f_m, chick_surv = TRUE) 
  {
    if(!chick_surv)
    {
      vitalrates<-list(
        F_Juv_survl = VR_list$F_Juv_survl,
        F_Adt_survl = VR_list$F_Adt_survl,
        M_Juv_survl = VR_list$M_Juv_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        h = VR_list$h,
        k = VR_list$k,
        PSR = VR_list$PSR,
        M2 = unname(freq_dep_ASR_f_m$SSD_M2),
        F2 = unname(freq_dep_ASR_f_m$SSD_F2))
      
      vitalrates_LTRE<-list(
        F_Juv_survl = VR_list$M_Juv_survl,
        F_Adt_survl = VR_list$M_Adt_survl,
        M_Juv_survl = VR_list$M_Juv_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        h = VR_list$h,
        k = VR_list$k,
        HSR = VR_list$HSR,
        M2 = unname(freq_dep_ASR$SSD_M2),
        F2 = unname(freq_dep_ASR$SSD_F2))
    }
    else
    {
      # make a list of all parameters
      vitalrates<-list(
        F_Chk_survl = VR_list$F_Chk_survl,
        F_Juv_survl = VR_list$F_Juv_survl,
        F_Adt_survl = VR_list$F_Adt_survl,
        M_Chk_survl = VR_list$M_Chk_survl,
        M_Juv_survl = VR_list$M_Juv_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        h = VR_list$h,
        k = VR_list$k,
        PSR = VR_list$PSR,
        M2 = unname(freq_dep_ASR_f_m$SSD_M2),
        F2 = unname(freq_dep_ASR_f_m$SSD_F2))
      
      vitalrates_LTRE<-list(
        F_Chk_survl = VR_list$M_Chk_survl,
        F_Juv_survl = VR_list$M_Juv_survl,
        F_Adt_survl = VR_list$M_Adt_survl,
        M_Chk_survl = VR_list$M_Chk_survl,
        M_Juv_survl = VR_list$M_Juv_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        h = VR_list$h,
        k = VR_list$k,
        HSR = VR_list$HSR,
        M2 = unname(freq_dep_ASR$SSD_M2),
        F2 = unname(freq_dep_ASR$SSD_F2))
    }
    if (is.vector(vitalrates)) {
      vitalrates <- as.list(vitalrates)
    }
    if (!is.list(vitalrates)) {
      stop("Vital rates should be a vector or list")
    }
    if (class(elements) != "expression") {
      stop("Matrix elements should be an expression")
    }
    n <- sqrt(length(elements))
    if (n%%1 != 0) {
      stop(paste("Length of element expression is", length(elements), 
                 "- Expecting power of 2 like 4,9,16 to form a square matrix"))
    }
    vrs <- try(sapply(elements, eval, vitalrates, NULL), silent = TRUE)
    if (class(vrs) == "try-error") {
      vrs <- sub("Error in eval\\(expr, envir, enclos\\) :", 
                 "", vrs[1])
      stop(paste("Cannot evaluate element expression using given vital rates:", 
                 vrs))
    }
    res <- data.frame(estimate = unlist(vitalrates), sensitivity = 0, 
                      elasticity = 0)
    A <- matrix(vrs, nrow = n, byrow = TRUE)
    ASR_f_m <- ASR_analysis_f_m(A)
    deriv.funcs <- sapply(elements, deriv, namevec = names(vitalrates), 
                          function.arg = TRUE)
    devs <- lapply(deriv.funcs, function(x) do.call(x, vitalrates))
    for (i in 1:length(vitalrates)) {
      derivs <- matrix(as.numeric(lapply(devs, function(x) attr(x, 
                                                                "gradient")[i])), nrow = n, byrow = TRUE)
      res[i, 2] <- sum(derivs * ASR_f_m$sensitivities)
      res[i, 3] <- vitalrates[[i]]/ASR_f_m$ASR_f_m * sum(derivs * ASR_f_m$sensitivities)
    }
    x <- res
    x$Vital_rate <- as.factor(rownames(x))
    if(!chick_surv)
    {
      y <- x[which(x$Vital_rate =="F_Juv_survl"|
                     x$Vital_rate =="M_Juv_survl"|
                     x$Vital_rate =="F_Adt_survl"|
                     x$Vital_rate =="M_Adt_survl"|
                     x$Vital_rate =="h"|
                     x$Vital_rate =="PSR"),]
      colnames(y) <- c("Estimate", "Sensitivity", "Elasticity", "Vital_rate")
      y_melt <- melt(y[,c(2:4)])
      y_melt$VR <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"Juv"), "Fledgling",
                                    ifelse(str_detect(y_melt$Vital_rate,"Adt"), "Adult",
                                           ifelse(str_detect(y_melt$Vital_rate,"PSR"), "1° sex ratio", "Fecundity"))))
      y_melt$VR <- factor(y_melt$VR, levels = c("Adult",
                                                "Fledgling",
                                                "Chick",
                                                " ",
                                                "1° sex ratio",
                                                "Fecundity"))
    }
    else
    {
      y <- x[which(x$Vital_rate == "F_Chk_survl"|
                     x$Vital_rate == "M_Chk_survl"|
                     x$Vital_rate == "F_Juv_survl"|
                     x$Vital_rate == "M_Juv_survl"|
                     x$Vital_rate == "F_Adt_survl"|
                     x$Vital_rate == "M_Adt_survl"|
                     x$Vital_rate == "h"|
                     x$Vital_rate == "PSR"),]
      colnames(y) <- c("Estimate", "Sensitivity", "Elasticity", "Vital_rate")
      y_melt <- melt(y[,c(2:4)])
      y_melt$VR <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"Chk"), "Chick",
                                    ifelse(str_detect(y_melt$Vital_rate,"Juv"), "Fledgling",
                                           ifelse(str_detect(y_melt$Vital_rate,"Adt"), "Adult",
                                                  ifelse(str_detect(y_melt$Vital_rate,"PSR"), "1° sex ratio", "Fecundity")))))
      y_melt$VR <- factor(y_melt$VR, levels = c("Adult",
                                                "Fledgling",
                                                "Chick",
                                                " ",
                                                "1° sex ratio",
                                                "Fecundity"))
    }
    y_melt$Sex <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"F_"), "Female", 
                                   ifelse(str_detect(y_melt$Vital_rate,"M_"),"Male", "Other")))
    y_melt$Sex <- 
      factor(y_melt$Sex,
             levels = c("Female","Male", "Other"))
    y_melt$value_trans <- ifelse(y_melt$Sex == "Female", abs(y_melt$value)*-1, y_melt$value)
    males <- subset(y_melt, (Sex == "Male" | Sex == "Other"))
    females <- subset(y_melt, Sex == "Female" )
    results <- list(males = males,
                    females = females)
    row.names(results$males) <- NULL
    row.names(results$females) <- NULL
    results$males[nrow(results$males)+1,] <- c("Dummy","Elasticity",0," ","Other",0)
    results$females[nrow(results$females)+1,] <- c("Dummy","Elasticity",0," ","Other",0)
    results$males[nrow(results$males)+1,] <- c("Dummy","Sensitivity",0," ","Other",0)
    results$females[nrow(results$females)+1,] <- c("Dummy","Sensitivity",0," ","Other",0)
    results$males$value_trans <- as.numeric(results$males$value_trans)
    results$females$value_trans <- as.numeric(results$females$value_trans)
    if(!chick_surv)
    {
      results$males[nrow(results$males)+1,] <- c("Dummy","Elasticity",0,"Chick","Male",0)
      results$females[nrow(results$females)+1,] <- c("Dummy","Elasticity",0,"Chick","Female",0)
      results$males[nrow(results$males)+1,] <- c("Dummy","Sensitivity",0,"Chick","Male",0)
      results$females[nrow(results$females)+1,] <- c("Dummy","Sensitivity",0,"Chick","Female",0)
      results$males$value_trans <- as.numeric(results$males$value_trans)
      results$females$value_trans <- as.numeric(results$females$value_trans)
    }
    results
  }

# LTRE analysis of ASR
vitalsens_ASR <- 
  function (elements, VR_list, freq_dep_ASR) 
    {
    # make a list of all parameters
    vitalrates<-list(
      F_Chk_survl = VR_list$F_Chk_survl,
      F_Juv_survl = VR_list$F_Juv_survl,
      F_Adt_survl = VR_list$F_Adt_survl,
      M_Chk_survl = VR_list$M_Chk_survl,
      M_Juv_survl = VR_list$M_Juv_survl,
      M_Adt_survl = VR_list$M_Adt_survl,
      h = VR_list$h,
      k = VR_list$k,
      HSR = VR_list$HSR,
      M2 = unname(freq_dep_ASR$SSD_M2),
      F2 = unname(freq_dep_ASR$SSD_F2))
    
    vitalrates_LTRE<-list(
      F_Chk_survl = VR_list$M_Chk_survl,
      F_Juv_survl = VR_list$M_Juv_survl,
      F_Adt_survl = VR_list$M_Adt_survl,
      M_Chk_survl = VR_list$M_Chk_survl,
      M_Juv_survl = VR_list$M_Juv_survl,
      M_Adt_survl = VR_list$M_Adt_survl,
      h = VR_list$h,
      k = VR_list$k,
      HSR = VR_list$HSR,
      M2 = unname(freq_dep_ASR$SSD_M2),
      F2 = unname(freq_dep_ASR$SSD_F2))
    
  if (is.vector(vitalrates)) {
    vitalrates <- as.list(vitalrates)
  }
  if (!is.list(vitalrates)) {
    stop("Vital rates should be a vector or list")
  }
  if (class(elements) != "expression") {
    stop("Matrix elements should be an expression")
  }
    if (is.vector(vitalrates_LTRE)) {
      vitalrates_LTRE <- as.list(vitalrates_LTRE)
    }
    if (!is.list(vitalrates_LTRE)) {
      stop("Vital rates should be a vector or list")
    }
  n <- sqrt(length(elements))
  if (n%%1 != 0) {
    stop(paste("Length of element expression is", length(elements), 
               "- Expecting power of 2 like 4,9,16 to form a square matrix"))
  }
  vrs <- try(sapply(elements, eval, vitalrates, NULL), silent = TRUE)
  vrs_LTRE <- try(sapply(elements, eval, vitalrates_LTRE, NULL), silent = TRUE)
  if (class(vrs) == "try-error") {
    vrs <- sub("Error in eval\\(expr, envir, enclos\\) :",
               "", vrs[1])
    stop(paste("Cannot evaluate element expression using given vital rates:",
               vrs))
  }
  if (class(vrs_LTRE) == "try-error") {
    vrs_LTRE <- sub("Error in eval\\(expr, envir, enclos\\) :",
               "", vrs_LTRE[1])
    stop(paste("Cannot evaluate element expression using given vital rates:",
               vrs_LTRE))
  }
  res <- data.frame(estimate = unlist(vitalrates), sensitivity = 0, 
                    elasticity = 0, LTRE = 0)
  A <- matrix(vrs, nrow = n, byrow = TRUE)
  A_LTRE <- matrix(vrs_LTRE, nrow = n, byrow = TRUE)
  Ac <- (A + A_LTRE)/2
  SAc <- ASR_analysis(Ac)
  ASR <- ASR_analysis(A)
  deriv.funcs <- sapply(elements, deriv, namevec = names(vitalrates), 
                        function.arg = TRUE)
  devs <- lapply(deriv.funcs, function(x) do.call(x, vitalrates))
  for (i in 1:length(vitalrates)) {
    derivs <- matrix(as.numeric(lapply(devs, function(x) attr(x, 
                                                              "gradient")[i])), nrow = n, byrow = TRUE)
    res[i, 2] <- sum(derivs * ASR$sensitivities)
    res[i, 3] <- vitalrates[[i]]/ASR$ASR * sum(derivs * ASR$sensitivities)
    res[i, 4] <- sum(derivs * SAc$sensitivities)
  }
  x <- res
  x$Vital_rate <- as.factor(rownames(x))
  if(!chick_surv)
  {
    y <- x[which(x$Vital_rate == "F_Juv_survl"|
                   x$Vital_rate == "M_Juv_survl"|
                   x$Vital_rate == "F_Adt_survl"|
                   x$Vital_rate == "M_Adt_survl"|
                   x$Vital_rate == "h"|
                   x$Vital_rate == "k"|
                   x$Vital_rate == "M2"|
                   x$Vital_rate == "F2"|
                   x$Vital_rate == "HSR"),]
    colnames(y) <- c("Estimate", "Sensitivity", "Elasticity", "LTRE", "Vital_rate")
    y_melt <- melt(y[,c(2:5)])
    y_melt$VR <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"Juv"), "Juvenile",
                                   ifelse(str_detect(y_melt$Vital_rate,"Adt"), "Adult",
                                          ifelse(str_detect(y_melt$Vital_rate,"HSR"), "Hatching sex ratio", 
                                                 ifelse(str_detect(y_melt$Vital_rate,"h"), "Mating System",
                                                        ifelse(str_detect(y_melt$Vital_rate,"k"), "Clutch size", "No. breeding adults"))))))
    y_melt$VR <- factor(y_melt$VR, levels = c("Adult",
                                              "Fledgling",
                                              "Chick",
                                              " ",
                                              "No. breeding adults",
                                              "Mating System",
                                              "Clutch size",
                                              "Hatching sex ratio"))
  }
  else
  {
  y <- x[which(x$Vital_rate == "F_Chk_survl"|
               x$Vital_rate == "M_Chk_survl"|
               x$Vital_rate == "F_Juv_survl"|
               x$Vital_rate == "M_Juv_survl"|
               x$Vital_rate == "F_Adt_survl"|
               x$Vital_rate == "M_Adt_survl"|
               x$Vital_rate == "h"|
               x$Vital_rate == "k"|
               x$Vital_rate == "M2"|
               x$Vital_rate == "F2"|
               x$Vital_rate == "HSR"),]
  colnames(y) <- c("Estimate", "Sensitivity", "Elasticity", "LTRE", "Vital_rate")
  y_melt <- melt(y[,c(2:5)])
  y_melt$VR <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"Chk"), "Chick",
                                ifelse(str_detect(y_melt$Vital_rate,"Juv"), "Fledgling",
                                       ifelse(str_detect(y_melt$Vital_rate,"Adt"), "Adult",
                                              ifelse(str_detect(y_melt$Vital_rate,"HSR"), "Hatching sex ratio", 
                                                     ifelse(str_detect(y_melt$Vital_rate,"h"), "Mating System",
                                                            ifelse(str_detect(y_melt$Vital_rate,"k"), "Clutch size", "No. breeding adults")))))))
  y_melt$VR <- factor(y_melt$VR, levels = c("Adult",
                                            "Fledgling",
                                            "Chick",
                                            " ",
                                            "No. breeding adults",
                                            "Mating System",
                                            "Clutch size",
                                            "Hatching sex ratio"))
  }
  y_melt$Sex <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"F"), "Female", 
                                 ifelse(str_detect(y_melt$Vital_rate,"M"),"Male", "Other")))
  y_melt$Sex <- 
    factor(y_melt$Sex,
           levels = c("Female","Male", "Other"))
  #y_melt$prop_value <- y_melt$value/sum(y_melt$value)
  y_melt$value_trans <- ifelse(y_melt$Sex == "Female", abs(y_melt$value)*-1, 
                               ifelse(y_melt$Vital_rate == "M2", abs(y_melt$value), y_melt$value))
  y_melt$value_trans <- ifelse(y_melt$Sex == "Other", abs(y_melt$value), y_melt$value_trans)
  males <- subset(y_melt, (Sex == "Male" | Sex == "Other"))
  females <- subset(y_melt, Sex == "Female" )
  results <- list(males = males,
                  females = females)
  row.names(results$males) <- NULL
  row.names(results$females) <- NULL
  results$males$value_trans <- as.numeric(results$males$value_trans)
  results$females$value_trans <- as.numeric(results$females$value_trans)
  results
  }

# LSA analysis of ASR
LSA_analysis <-
  function(ASR, survival_rates)
  {
    Female_Adult <- filter(survival_rates, Sex_Age == "Female_Adult")
    Female_Chick <- filter(survival_rates, Sex_Age == "Female_Chick")
    Female_Juvenile <- filter(survival_rates, Sex_Age == "Female_Juvenile")
    Male_Adult <- filter(survival_rates, Sex_Age == "Male_Adult")
    Male_Chick <- filter(survival_rates, Sex_Age == "Male_Chick")
    Male_Juvenile <- filter(survival_rates, Sex_Age == "Male_Juvenile")
    
    r_Female_Adult <- summary(lm(Female_Adult$estimate ~ ASR$ASR))$r.squared
    r_Female_Chick <- summary(lm(Female_Chick$estimate ~ ASR$ASR))$r.squared
    r_Female_Juvenile <- summary(lm(Female_Juvenile$estimate ~ ASR$ASR))$r.squared
    r_Male_Adult <- summary(lm(Male_Adult$estimate ~ ASR$ASR))$r.squared
    r_Male_Chick <- summary(lm(Male_Chick$estimate ~ ASR$ASR))$r.squared
    r_Male_Juvenile <- summary(lm(Male_Juvenile$estimate ~ ASR$ASR))$r.squared
    
    r_estimates <- c(r_Female_Adult, r_Female_Chick, r_Female_Juvenile,
                     r_Male_Adult, r_Male_Chick, r_Male_Juvenile)
    Stages <- levels(survival_rates$Sex_Age)
    LSA_table <- as.data.frame(cbind(Stages, r_estimates))
    LSA_table$r_estimates <- as.numeric(as.character(LSA_table$r_estimates))
    LSA_table$VR <- as.factor(ifelse(str_detect(LSA_table$Stages,"Chick"), "Chick",
                                  ifelse(str_detect(LSA_table$Stages,"Juvenile"), "Fledgling", "Adult")))
    LSA_table$Sex <- as.factor(ifelse(str_detect(LSA_table$Stages,"F"), "Female", "Male"))
    LSA_table$value_trans <- ifelse(LSA_table$Sex == "Female", abs(LSA_table$r_estimates)*-1, LSA_table$r_estimates)
    LSA_table$VR <- factor(LSA_table$VR, levels = c("Adult",
                                                    "Fledgling",
                                                    "Chick"))
    LSA_table
  }

# plotting functions
Elasticity_line_plot <- 
  function(lower_level_sens_analysis_result)
  {
    Elasticity_line_plot <- ggplot(lower_level_sens_analysis_result$Elasticities,
                                   aes(x = Perturbation, y = Elasticity, group = Vitalrate)) +  
      theme_bw() +
      geom_line(size = 1, aes(colour = Vitalrate)) +
      theme(text=element_text(family="Arial"),
            legend.text=element_text(size=11),
            legend.title=element_text(size=12),
            legend.key.height=unit(0.8,"line"),
            legend.key.width=unit(0.8,"line"),
            legend.background = element_rect(fill=NA),
            axis.title.x = element_text(size=12, vjust=-0.1),
            axis.text.x  = element_text(size=11), 
            axis.title.y = element_text(size=12, vjust=1.2),
            axis.text.y  = element_text(size=11), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
            axis.ticks.length = unit(0.2, "cm"),
            axis.ticks.x = element_line(size = 0.5, colour = "grey40")) +
      xlab("Proportion of current vital rate") +
      ylab("Proportion of current ASR") +
      scale_y_continuous(limits = c(0, 1.5)) +
      scale_colour_manual(values = cbPalette, name = "Model components")
    Elasticity_line_plot
  }

Sensitivity_line_plot <- 
  function(lower_level_sens_analysis_result)
  {
    Sensitivity_line_plot <- ggplot(lower_level_sens_analysis_result$Sensitivities, 
                                    aes(x = Perturbation, y = Sensitivity, group = Vitalrate)) +  
      theme_bw() +
      geom_line(size = 1, aes(colour = Vitalrate)) +
      theme(text=element_text(family="Arial"),
            legend.text=element_text(size=11),
            legend.title=element_text(size=12),
            legend.key.height=unit(0.8,"line"),
            legend.key.width=unit(0.8,"line"),
            legend.background = element_rect(fill=NA),
            axis.title.x = element_text(size=12, vjust=-0.1),
            axis.text.x  = element_text(size=11), 
            axis.title.y = element_text(size=12, vjust=1.2),
            axis.text.y  = element_text(size=11), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
            axis.ticks.length = unit(0.2, "cm"),
            axis.ticks.x = element_line(size = 0.5, colour = "grey40")) +
      xlab("Value of vital rate") +
      ylab("Adult sex ratio (proportion \u2642)") +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
#       annotate("pointrange", x = c(0.1541992, 0.2282648), y = c(0.7426747, 0.7426747), 
#                ymin = c(0.7426747, 0.7426747), ymax = c(0.7426747, 0.7426747),
#                colour = cbPalette, size = 1) +
      scale_colour_manual(values=cbPalette, name="Model components")
    Sensitivity_line_plot
  }

Elasticity_bar_plot <- 
  function(vitalsens_ASR_result)
    {
      Elasticity_bar_plot <- 
        ggplot() +  
        theme_bw() +
        geom_bar(data = vitalsens_ASR_result$males[which(vitalsens_ASR_result$males$variable == "Elasticity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity", alpha = 0.8) +
        geom_bar(data = vitalsens_ASR_result$females[which(vitalsens_ASR_result$females$variable == "Elasticity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity", alpha = 0.8) +
        coord_flip() +
        theme(text=element_text(family="Arial"),
              legend.position="none",
              legend.position = c(0, 1), 
              legend.justification = c(0, 1),
              legend.text=element_text(size=11),
              legend.title=element_blank(),
              legend.key.height=unit(0.8,"line"),
              legend.key.width=unit(0.8,"line"),
              legend.background = element_rect(fill=NA),
              axis.title.x = element_text(size=12, vjust=-0.1, margin = margin(10, 0, 0, 0)),
              axis.text.x  = element_text(size=10), 
              axis.title.y = element_text(size=12, margin = margin(0, 15, 0, 0)),
              axis.text.y  = element_text(size=10),#, angle=90, hjust = 0.5, vjust = 0.5),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_rect(linetype = "solid", colour = "grey"),
              panel.grid.minor = element_blank(),
              axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
              axis.ticks.length = unit(0.2, "cm"),
              axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
              plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) +
        ylab("Elasticity of adult sex ratio") +
        xlab("Apparent survival (\u03D5)          Mating function components") +
        scale_fill_manual(values = c(brewer.pal(8, "Set1")[c(1,2)], "#000000")) +
        scale_y_continuous(limits = c(-0.5, 0.5)) +
        annotate("text", x = 4.3, y = -0.45, vjust = 1.5, hjust = 0.5, label = "\u2640", 
                 size = 7, family="Arial") +
        annotate("text", x = 4.3, y = 0.45, vjust = 1.5, hjust = 0.5, label = "\u2642", 
                 size = 7, family="Arial") +
        annotate("rect", xmin = -Inf, xmax = 3.50, ymin = -Inf, ymax = Inf,
                 alpha = .2)
    Elasticity_bar_plot
    }

Sensitivity_bar_plot <- 
  function(vitalsens_ASR_result)
  {
      Sensitivity_bar_plot <- 
        ggplot() +  
        theme_bw() +
        geom_bar(data = vitalsens_ASR_result$males[which(vitalsens_ASR_result$males$variable == "Sensitivity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity", alpha = 0.8) +
        geom_bar(data = vitalsens_ASR_result$females[which(vitalsens_ASR_result$females$variable == "Sensitivity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity", alpha = 0.8) +
        coord_flip() +
        theme(text=element_text(family="Arial"),
              legend.position="none",
              legend.position = c(0, 1), 
              legend.justification = c(0, 1),
              legend.text=element_text(size=11),
              legend.title=element_blank(),
              legend.key.height=unit(0.8,"line"),
              legend.key.width=unit(0.8,"line"),
              legend.background = element_rect(fill=NA),
              axis.title.x = element_text(size=12, vjust=-0.1, margin = margin(10, 0, 0, 0)),
              axis.text.x  = element_text(size=10), 
              axis.title.y = element_text(size=12, margin = margin(0, 15, 0, 0)),
              axis.text.y  = element_text(size=10, hjust = 1),#, angle=90, hjust = 0.5, vjust = 0.5),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_rect(linetype = "solid", colour = "grey"),
              panel.grid.minor = element_blank(),
              axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
              axis.ticks.length = unit(0.2, "cm"),
              axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
              plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) +
        ylab("Sensitivity of adult sex ratio") +
        xlab("Apparent survival (\u03D5)          Mating function components") +
        scale_fill_manual(values = c(brewer.pal(8, "Set1")[c(1,2)], "#000000")) +
        scale_y_continuous(limits = c(-0.5, 0.5)) +
        annotate("text", x = 4.3, y = -0.45, vjust = 1.5, hjust = 0.5, label = "\u2640", 
                 size = 7, family="Arial") +
        annotate("text", x = 4.3, y = 0.45, vjust = 1.5, hjust = 0.5, label = "\u2642", 
                 size = 7, family="Arial") +
        annotate("rect", xmin = -Inf, xmax = 3.50, ymin = -Inf, ymax = Inf,
                 alpha = .2)
    Sensitivity_bar_plot
  }

LTRE_bar_plot <- 
  function(vitalsens_ASR_result)
  {
    LTRE_bar_plot <- 
      ggplot() +  
      theme_bw() +
      geom_bar(data = vitalsens_ASR_result$males[which(vitalsens_ASR_result$males$variable == "LTRE"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity", alpha = 0.8) +
      #geom_bar(data = vitalsens_ASR_result$females[which(vitalsens_ASR_result$females$variable == "LTRE"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity", alpha = 0.8) +
      coord_flip() +
      theme(text=element_text(family="Arial"),
            legend.position="none",
            legend.position = c(0, 1), 
            legend.justification = c(0, 1),
            legend.text=element_text(size=11),
            legend.title=element_blank(),
            legend.key.height=unit(0.8,"line"),
            legend.key.width=unit(0.8,"line"),
            legend.background = element_rect(fill=NA),
            axis.title.x = element_text(size=12, vjust=-0.1, margin = margin(10, 0, 0, 0)),
            axis.text.x  = element_text(size=10), 
            axis.title.y = element_text(size=12, margin = margin(0, 15, 0, 0)),
            axis.text.y  = element_text(size=10, hjust = 1),#, angle=90, hjust = 0.5, vjust = 0.5),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(linetype = "solid", colour = "grey"),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
            axis.ticks.length = unit(0.2, "cm"),
            axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
            plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) +
      ylab("Contribution to adult sex ratio") +
      xlab("Apparent survival (\u03D5)          Mating function components") +
      #scale_fill_brewer(palette = "Set1") +
      scale_fill_manual(values = c(brewer.pal(8, "Set1")[c(2)], "#000000")) +
      scale_y_continuous(limits = c(-0.5, 0.5)) +
      annotate("text", x = 4.3, y = -0.45, vjust = 1.5, hjust = 0.5, label = "\u2640", 
               size = 7, family="Arial") +
      annotate("text", x = 4.3, y = 0.45, vjust = 1.5, hjust = 0.5, label = "\u2642", 
               size = 7, family="Arial") +
      annotate("rect", xmin = -Inf, xmax = 3.50, ymin = -Inf, ymax = Inf,
               alpha = .2)
    LTRE_bar_plot
  }

LSA_bar_plot <- 
  function(LSA_result)
  {
    LSA_bar_plot <- 
      ggplot() +  
      theme_bw() +
      geom_bar(data = LSA_result[which(LSA_result$Sex == "Male"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity", alpha = 0.8) +
      geom_bar(data = LSA_result[which(LSA_result$Sex == "Female"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity", alpha = 0.8) +
      coord_flip() +
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
            axis.title.y = element_text(size=12, margin = margin(0, 61.5, 0, 0)),
            axis.text.y  = element_text(size=10, hjust = 1),#, angle=90, hjust = 0.5, vjust = 0.5),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(linetype = "solid", colour = "grey"),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
            axis.ticks.length = unit(0.2, "cm"),
            axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
            plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) +
      ylab(expression(paste("Adult sex ratio ", italic(r^"2")))) +
      xlab("Apparent survival (\u03D5)                                                      ") +
      scale_fill_manual(values = c(brewer.pal(8, "Set1")[c(1,2)], "#000000")) +
      scale_y_continuous(limits = c(-0.5, 0.5)) +
      annotate("text", x = 4.3, y = -0.45, vjust = 1.5, hjust = 0.5, label = "\u2640", 
               size = 7, family="Arial") +
      annotate("text", x = 4.3, y = 0.45, vjust = 1.5, hjust = 0.5, label = "\u2642", 
               size = 7, family="Arial") +
      annotate("text", x = 4, y = 0, vjust = 0.5, hjust = 0.5, label = "NA", 
               size = 5, family="Arial") +
      annotate("text", x = 5, y = 0, vjust = 0.5, hjust = 0.5, label = "NA", 
               size = 5, family="Arial") +
      annotate("text", x = 6, y = 0, vjust = 0.5, hjust = 0.5, label = "NA", 
               size = 5, family="Arial") +
      annotate("text", x = 7, y = 0, vjust = 0.5, hjust = 0.5, label = "NA", 
               size = 5, family="Arial") +
      annotate("text", x = 7.6, y = 0, vjust = 0.5, hjust = 0.5, label = "NA", 
               size = 5, family="Arial", color = "white") +
      annotate("rect", xmin = -Inf, xmax = 3.50, ymin = -Inf, ymax = Inf,
               alpha = .2)
    LSA_bar_plot
  }

###############################################################################
# implement the functions on the bootstrapped data output

Ceuta_surv_boot <- read.table(file = "output/Bootstrap/Ceuta_Survival_rates.txt", header = TRUE)
Ceuta_ASR_boot <- read.table(file = "output/Bootstrap/Ceuta_ASR.txt", header = TRUE)

Maio_surv_boot <- read.table(file = "output/Bootstrap/Maio_Survival_rates.txt", header = TRUE)
Maio_ASR_boot <- read.table(file = "output/Bootstrap/Maio_ASR.txt", header = TRUE)

Ceuta_surv_boot$iter <- as.factor(Ceuta_surv_boot$iter)
Maio_surv_boot$iter <- as.factor(Maio_surv_boot$iter)

Ceuta_avg_rates <- 
  Ceuta_surv_boot %>%
  group_by(Sex_Age) %>%
  summarise(Avg = mean(estimate))

Ceuta_avg_rates <- as.data.frame(Ceuta_avg_rates)

Maio_avg_rates <- 
  Maio_surv_boot %>%
  group_by(Sex_Age) %>%
  summarise(Avg = mean(estimate))

Maio_avg_rates <- as.data.frame(Maio_avg_rates)

# Define Ceuta vital rates estimated from mark-recapture analysis:
Ceuta_VR <- list(F_Chk_survl = Ceuta_avg_rates[2,2],
                 F_Juv_survl = Ceuta_avg_rates[3,2],
                 F_Adt_survl = Ceuta_avg_rates[1,2],
                 M_Chk_survl = Ceuta_avg_rates[5,2],
                 M_Juv_survl = Ceuta_avg_rates[6,2],
                 M_Adt_survl = Ceuta_avg_rates[4,2],
                 
# Define h (harem size, h = 1 is monogamy) and k (clutch size)
                 h = 1,
                 k = 3,
# Define primary sex ratio (assumed to be 0.5)
                 HSR = 0.5)

# Define Maio vital rates estimated from mark-recapture analysis:
Maio_VR <- list(F_Juv_survl = Maio_avg_rates[3,2],
                 F_Adt_survl = Maio_avg_rates[1,2],
                 M_Juv_survl = Maio_avg_rates[6,2],
                 M_Adt_survl = Maio_avg_rates[4,2],
                 
                 # Define h (harem size, h = 1 is monogamy) and k (clutch size)
                 h = 1,
                 k = 3,
                 # Define primary sex ratio (assumed to be 0.5)
                 HSR = 0.5)

# Ceuta matrix:
Ceuta_matrix <- plover_matrix(Ceuta_VR)
Maio_matrix <- plover_matrix(Maio_VR, chick_surv = FALSE)


# Determine the ASR at the stable stage distribution (assume h is 1, which 
# is monogamy...a conservative assumption)
Ceuta_ASR_h_1 <- freq_dep_SSD_ASR(A = Ceuta_matrix, h = 1, k = 3)
Maio_ASR_h_1 <- freq_dep_SSD_ASR(A = Maio_matrix, h = 1, k = 3)

# Lower-level vital rate sensitivity analysis
Ceuta_LLSA <- lower_level_sens_analysis(freq_dep_ASR = Ceuta_ASR_h_1, VR_list = Ceuta_VR)
Maio_LLSA <- lower_level_sens_analysis(freq_dep_ASR = Maio_ASR_h_1, VR_list = Maio_VR, chick_surv = FALSE)

# Calculate vital rate sensitivities and elasticities
Ceuta_VR_Sens_Elas <- vitalsens_ASR(elements = Ceuta_LLSA$Element_expression,
                                    VR_list = Ceuta_VR, freq_dep_ASR = Ceuta_ASR_h_1)

Maio_VR_Sens_Elas <- vitalsens_ASR(elements = Maio_LLSA$Element_expression,
                                    VR_list = Maio_VR, freq_dep_ASR = Maio_ASR_h_1)

Ceuta_LSA <- LSA_analysis(ASR = Ceuta_ASR_boot, survival_rates = Ceuta_surv_boot)

# define color palette
cbPalette <- c(brewer.pal(9,"Reds")[c(5,6,8)], brewer.pal(9,"Blues")[c(5,6,8)],
               brewer.pal(8, "Dark2")[c(1,8,5,2,7)])

# Plot linear sex- and stage-specific elasticities
Ceuta_Elasticity_line_plot <- Elasticity_line_plot(Ceuta_LLSA)

# Plot linear sex- and stage-specific sensitivities
Ceuta_Sensitivity_line_plot <- Sensitivity_line_plot(Ceuta_LLSA)

# Plot overall sex- and stage-specific elasticities
Ceuta_Elasticity_bar_plot <- Elasticity_bar_plot(Ceuta_VR_Sens_Elas)

# Plot overall sex- and stage-specific sensitivities
Ceuta_Sensitivity_bar_plot <- Sensitivity_bar_plot(Ceuta_VR_Sens_Elas)

# Plot overall sex- and stage-specific LTRE
Ceuta_LTRE_bar_plot <- LTRE_bar_plot(Ceuta_VR_Sens_Elas)

# Inspect LTRE values
Ceuta_VR_Sens_Elas$males[which(Ceuta_VR_Sens_Elas$males$variable == "LTRE"),]

# Plot overall sex- and stage-specific LTRE
Ceuta_LSA_bar_plot <- LSA_bar_plot(Ceuta_LSA)