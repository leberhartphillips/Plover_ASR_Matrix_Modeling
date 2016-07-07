# Plotting script of bootstrap results for comparative analysis
# Sex-differences in survival
# Luke J. Eberhart-Phillips
# July 7, 2016

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

# functions for perturbation analysis
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
        HSR = VR_list$HSR,
        M2 = unname(freq_dep_ASR$SSD_M2),
        F2 = unname(freq_dep_ASR$SSD_F2))
      # make a matrix of the elements
      el <- expression(0, ((k*M2)/(M2+(F2/h)))*(1-HSR), 0, ((k*F2)/(M2+(F2/h)))*(1-HSR),
                       F_Juv_survl, F_Adt_survl, 0, 0,
                       0, ((k*M2)/(M2+(F2/h)))*HSR, 0, ((k*F2)/(M2+(F2/h)))*HSR,
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
        HSR = VR_list$HSR,
        M2 = unname(freq_dep_ASR$SSD_M2),
        F2 = unname(freq_dep_ASR$SSD_F2))
      # make a matrix of the elements
      el <- expression(0, ((k*M2)/(M2+(F2/h)))*(1-HSR), 0, ((k*F2)/(M2+(F2/h)))*(1-HSR),
                       (F_Chk_survl*F_Juv_survl), F_Adt_survl, 0, 0,
                       0, ((k*M2)/(M2+(F2/h)))*HSR, 0, ((k*F2)/(M2+(F2/h)))*HSR,
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
        vrsen[i,h] <- eigen(A)$vectors[4,1]/(eigen(A)$vectors[2,1]+eigen(A)$vectors[4,1])
      }
    }
    # Next calculate rescaled elasticities
    vrelas <- matrix(numeric(n*length(vr_nums)), ncol=n, dimnames=list(vr_nums, names(vr)))
    for (h in 1:n)
    {
      for (i in 1:length(vr_nums))
      {
        vr2 <- vr
        vr2[[h]]<-vr_nums[i]*vr2[[h]]
        A<-matrix(sapply(el, eval, vr2 , NULL), nrow=sqrt(length(el)), byrow=TRUE)
        vrelas[i,h] <- (eigen(A)$vectors[4,1]/(eigen(A)$vectors[2,1]+eigen(A)$vectors[4,1]))/unname(freq_dep_ASR$ASR)
      }
    }
    if(!chick_surv)
    {
      # tidy up and label results
      colnames(vrsen) <- c("Female juvenile survival", 
                           "Female adult survival",
                           "Male fledgling survival", 
                           "Male adult survival",
                           "Mating system index (h)", 
                           "Clutch size",
                           "Hatching sex ratio",
                           "Breeding males",
                           "Breeding females")
      colnames(vrelas) <- c("Female juvenile survival", 
                            "Female adult survival", 
                            "Male fledgling survival", 
                            "Male adult survival",
                            "Mating system (h)", 
                            "Clutch size",
                            "Hatching sex ratio",
                            "Breeding males",
                            "Breeding females")
    }
    else
    {
      # tidy up and label results
      colnames(vrsen) <- c("Female chick survival", 
                           "Female fledgling survival", 
                           "Female adult survival", 
                           "Male chick survival", 
                           "Male fledgling survival", 
                           "Male adult survival",
                           "Mating system index (h)", 
                           "Clutch size",
                           "Hatching sex ratio",
                           "Breeding males",
                           "Breeding females")
      colnames(vrelas) <- c("Female chick survival", 
                            "Female fledgling survival", 
                            "Female adult survival", 
                            "Male chick survival", 
                            "Male fledgling survival", 
                            "Male adult survival",
                            "Mating system (h)", 
                            "Clutch size",
                            "Hatching sex ratio",
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
    ASR <- stable.stage[2]/stable.stage[4] # SSD ASR
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

vitalsens_ASR <- 
  function (elements, VR_list, freq_dep_ASR, chick_surv = TRUE, pop_name) 
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
        HSR = VR_list$HSR,
        M2 = unname(freq_dep_ASR$SSD_M2),
        F2 = unname(freq_dep_ASR$SSD_F2))
      
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
                      elasticity = 0)
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
      y <- x[which(x$Vital_rate =="F_Juv_survl"|
                     x$Vital_rate =="M_Juv_survl"|
                     x$Vital_rate =="F_Adt_survl"|
                     x$Vital_rate =="M_Adt_survl"|
                     x$Vital_rate =="h"|
                     x$Vital_rate == "k"|
                     x$Vital_rate == "M2"|
                     x$Vital_rate == "F2"|
                     x$Vital_rate =="HSR"),]
      colnames(y) <- c("Estimate", "Sensitivity", "Elasticity", "LTRE", "Vital_rate")
      y_melt <- melt(y[,c(2:5)])
      y_melt$VR <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"Juv"), "Juvenile",
                                    ifelse(str_detect(y_melt$Vital_rate,"Adt"), "Adult",
                                           ifelse(str_detect(y_melt$Vital_rate,"HSR"), "Hatching sex ratio", 
                                                  ifelse(str_detect(y_melt$Vital_rate,"h"), "Mating System",
                                                         ifelse(str_detect(y_melt$Vital_rate,"k"), "Clutch size", "No. breeding adults"))))))
      y_melt$VR <- factor(y_melt$VR, levels = c("Adult",
                                                "Juvenile",
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
                     x$Vital_rate =="h"|
                     x$Vital_rate == "k"|
                     x$Vital_rate == "M2"|
                     x$Vital_rate == "F2"|
                     x$Vital_rate =="HSR"),]
      colnames(y) <- c("Estimate", "Sensitivity", "Elasticity", "LTRE", "Vital_rate")
      y_melt <- melt(y[,c(2:5)])
      y_melt$VR <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"Chk"), "Chick",
                                    ifelse(str_detect(y_melt$Vital_rate,"Juv"), "Juvenile",
                                           ifelse(str_detect(y_melt$Vital_rate,"Adt"), "Adult",
                                                  ifelse(str_detect(y_melt$Vital_rate,"HSR"), "Hatching sex ratio", 
                                                         ifelse(str_detect(y_melt$Vital_rate,"h"), "Mating System",
                                                                ifelse(str_detect(y_melt$Vital_rate,"k"), "Clutch size", "No. breeding adults")))))))
      y_melt$VR <- factor(y_melt$VR, levels = c("Adult",
                                                "Juvenile",
                                                "Chick",
                                                " ",
                                                "No. breeding adults",
                                                "Mating System",
                                                "Clutch size",
                                                "Hatching sex ratio"))
    }
    y_melt$Sex <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"F_"), "Female", 
                                   ifelse(str_detect(y_melt$Vital_rate,"M_"),"Male", "Other")))
    y_melt$Sex <- 
      factor(y_melt$Sex,
             levels = c("Female","Male", "Other"))
    y_melt$population <- pop_name
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

###############################################################################
# CEUTA SNOWY PLOVER 
#--------------------
# implement the functions on the bootstrapped data output
Ceuta_surv_boot <- read.table(file = "output/Bootstrap/Ceuta_Survival_rates.txt", header = TRUE)
Ceuta_ASR_boot <- read.table(file = "output/Bootstrap/Ceuta_ASR.txt", header = TRUE)

Ceuta_surv_boot$iter <- as.factor(Ceuta_surv_boot$iter)

Ceuta_avg_rates <- 
  Ceuta_surv_boot %>%
  group_by(Sex_Age) %>%
  summarise(Avg = mean(estimate))

Ceuta_avg_rates <- as.data.frame(Ceuta_avg_rates)

# Define Ceuta vital rates estimated from mark-recapture analysis:
Ceuta_VR <- list(F_Chk_survl = Ceuta_avg_rates[2,2],
                 F_Juv_survl = Ceuta_avg_rates[3,2],
                 F_Adt_survl = Ceuta_avg_rates[1,2],
                 M_Chk_survl = Ceuta_avg_rates[5,2],
                 M_Juv_survl = Ceuta_avg_rates[6,2],
                 M_Adt_survl = Ceuta_avg_rates[4,2],
                 # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                 h = 1,
                 k = 3,
                 # Define primary sex ratio (assumed to be 0.5)
                 HSR = 0.5)

# Ceuta matrix:
Ceuta_matrix <- plover_matrix(Ceuta_VR)

# Determine the ASR at the stable stage distribution (assume h is 0.3, which 
# is polyandry with females having on average 3 males)
Ceuta_ASR_h_1 <- freq_dep_SSD_ASR(A = Ceuta_matrix, h = 1, k = 3)

# Lower-level vital rate sensitivity analysis
Ceuta_LLSA <- lower_level_sens_analysis(freq_dep_ASR = Ceuta_ASR_h_1, VR_list = Ceuta_VR)

# Calculate vital rate sensitivities and elasticities
Ceuta_VR_Sens_Elas <- vitalsens_ASR(elements = Ceuta_LLSA$Element_expression,
                                    VR_list = Ceuta_VR, freq_dep_ASR = Ceuta_ASR_h_1,
                                    pop_name = "Snowy")

###############################################################################
# TUZLA KENTISH PLOVER
# --------------------

Tuzla_surv_boot <- read.table(file = "output/Bootstrap/Tuzla_Survival_rates.txt", header = TRUE)
Tuzla_ASR_boot <- read.table(file = "output/Bootstrap/Tuzla_ASR.txt", header = TRUE)

Tuzla_surv_boot$iter <- as.factor(Tuzla_surv_boot$iter)

Tuzla_avg_rates <- 
  Tuzla_surv_boot %>%
  group_by(Sex_Age) %>%
  summarise(Avg = mean(estimate))

Tuzla_avg_rates <- as.data.frame(Tuzla_avg_rates)

# Define Tuzla vital rates estimated from mark-recapture analysis:
Tuzla_VR <- list(F_Juv_survl = Tuzla_avg_rates[3,2],
                 F_Adt_survl = Tuzla_avg_rates[1,2],
                 M_Juv_survl = Tuzla_avg_rates[6,2],
                 M_Adt_survl = Tuzla_avg_rates[4,2],
  # Define h (harem size, h < 1 is polyandry) and k (clutch size)
  h = 1,
  k = 3,
  # Define primary sex ratio (assumed to be 0.5)
  HSR = 0.5)

# Tuzla matrix:
Tuzla_matrix <- plover_matrix(Tuzla_VR, chick_surv = FALSE)

# Determine the ASR at the stable stage distribution
Tuzla_ASR_h_1 <- freq_dep_SSD_ASR(A = Tuzla_matrix, h = 1, k = 3)

# Lower-level vital rate sensitivity analysis
Tuzla_LLSA <- lower_level_sens_analysis(freq_dep_ASR = Tuzla_ASR_h_1,
                                        VR_list = Tuzla_VR, chick_surv = FALSE)

# Calculate vital rate sensitivities and elasticities
Tuzla_VR_Sens_Elas <- vitalsens_ASR(elements = Tuzla_LLSA$Element_expression,
                                    VR_list = Tuzla_VR, freq_dep_ASR = Tuzla_ASR_h_1,
                                    chick_surv = FALSE,
                                    pop_name = "Kentish-Tuzla")

###############################################################################
# MADAGASCAR PLOVER
# --------------------
# Define Madagascar plover vital rates estimated from mark-recapture analysis:
MP_surv_boot <- read.table(file = "output/Bootstrap/MP_Survival_rates.txt", header = TRUE)
MP_ASR_boot <- read.table(file = "output/Bootstrap/MP_ASR.txt", header = TRUE)

MP_surv_boot$iter <- as.factor(MP_surv_boot$iter)

MP_avg_rates <- 
  MP_surv_boot %>%
  group_by(Sex_Age) %>%
  summarise(Avg = mean(estimate))

MP_avg_rates <- as.data.frame(MP_avg_rates)

# Define MP vital rates estimated from mark-recapture analysis:
MP_VR <- list(F_Juv_survl = MP_avg_rates[3,2],
                 F_Adt_survl = MP_avg_rates[1,2],
                 M_Juv_survl = MP_avg_rates[6,2],
                 M_Adt_survl = MP_avg_rates[4,2],
                 # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                 h = 1,
                 k = 2,
                 # Define primary sex ratio (assumed to be 0.5)
                 HSR = 0.5)

# MP matrix:
MP_matrix <- plover_matrix(MP_VR, chick_surv = FALSE)

# Determine the ASR at the stable stage distribution
MP_ASR_h_1 <- freq_dep_SSD_ASR(A = MP_matrix, h = 1, k = 2)

# Lower-level vital rate sensitivity analysis
MP_LLSA <- lower_level_sens_analysis(freq_dep_ASR = MP_ASR_h_1,
                                        VR_list = MP_VR, chick_surv = FALSE)

# Calculate vital rate sensitivities and elasticities
MP_VR_Sens_Elas <- vitalsens_ASR(elements = MP_LLSA$Element_expression,
                                    VR_list = MP_VR, freq_dep_ASR = MP_ASR_h_1,
                                    chick_surv = FALSE,
                                 pop_name = "Madagascar")
###############################################################################
# WHITE-FRONTED PLOVER
# --------------------
# Define White_fronted plover vital rates estimated from mark-recapture analysis:
WfP_surv_boot <- read.table(file = "output/Bootstrap/WfP_Survival_rates.txt", header = TRUE)
WfP_ASR_boot <- read.table(file = "output/Bootstrap/WfP_ASR.txt", header = TRUE)

WfP_surv_boot$iter <- as.factor(WfP_surv_boot$iter)

WfP_avg_rates <- 
  WfP_surv_boot %>%
  group_by(Sex_Age) %>%
  summarise(Avg = mean(estimate))

WfP_avg_rates <- as.data.frame(WfP_avg_rates)

# Define WfP vital rates estimated from mark-recapture analysis:
WfP_VR <- list(F_Juv_survl = WfP_avg_rates[3,2],
                 F_Adt_survl = WfP_avg_rates[1,2],
                 M_Juv_survl = WfP_avg_rates[6,2],
                 M_Adt_survl = WfP_avg_rates[4,2],
                 # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                 h = 1,
                 k = 2,
                 # Define primary sex ratio (assumed to be 0.5)
                 HSR = 0.5)

# WfP matrix:
WfP_matrix <- plover_matrix(WfP_VR, chick_surv = FALSE)

# Determine the ASR at the stable stage distribution
WfP_ASR_h_1 <- freq_dep_SSD_ASR(A = WfP_matrix, h = 1, k = 2)

# Lower-level vital rate sensitivity analysis
WfP_LLSA <- lower_level_sens_analysis(freq_dep_ASR = WfP_ASR_h_1,
                                        VR_list = WfP_VR, chick_surv = FALSE)

# Calculate vital rate sensitivities and elasticities
WfP_VR_Sens_Elas <- vitalsens_ASR(elements = WfP_LLSA$Element_expression,
                                    VR_list = WfP_VR, freq_dep_ASR = WfP_ASR_h_1,
                                    chick_surv = FALSE,
                                  pop_name = "White-fronted")

###############################################################################
# MAIO KENTISH PLOVER
# --------------------
# Define Maio plover vital rates estimated from mark-recapture analysis:
Maio_surv_boot <- read.table(file = "output/Bootstrap/Maio_Survival_rates.txt", header = TRUE)
Maio_ASR_boot <- read.table(file = "output/Bootstrap/Maio_ASR.txt", header = TRUE)

Maio_surv_boot$iter <- as.factor(Maio_surv_boot$iter)

Maio_avg_rates <- 
  Maio_surv_boot %>%
  group_by(Sex_Age) %>%
  summarise(Avg = mean(estimate))

Maio_avg_rates <- as.data.frame(Maio_avg_rates)

# Define Maio vital rates estimated from mark-recapture analysis:
Maio_VR <- list(F_Juv_survl = Maio_avg_rates[3,2],
                 F_Adt_survl = Maio_avg_rates[1,2],
                 M_Juv_survl = Maio_avg_rates[6,2],
                 M_Adt_survl = Maio_avg_rates[4,2],
                 # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                 h = 1,
                 k = 3,
                 # Define primary sex ratio (assumed to be 0.5)
                 HSR = 0.5)

# Maio matrix:
Maio_matrix <- plover_matrix(Maio_VR, chick_surv = FALSE)

# Determine the ASR at the stable stage distribution (assume h is 0.5, which 
# is polyandry with females having on average 2 males)
Maio_ASR_h_1 <- freq_dep_SSD_ASR(A = Maio_matrix, h = 1, k = 3)

# Lower-level vital rate sensitivity analysis
Maio_LLSA <- lower_level_sens_analysis(freq_dep_ASR = Maio_ASR_h_1,
                                        VR_list = Maio_VR, chick_surv = FALSE)

# Calculate vital rate sensitivities and elasticities
Maio_VR_Sens_Elas <- vitalsens_ASR(elements = Maio_LLSA$Element_expression,
                                    VR_list = Maio_VR, freq_dep_ASR = Maio_ASR_h_1,
                                    chick_surv = FALSE,
                                   pop_name = "Kentish-Maio")

###############################################################################
# KITTLITZ'S PLOVER
# --------------------
# Define Kittlitz plover vital rates estimated from mark-recapture analysis:
KiP_surv_boot <- read.table(file = "output/Bootstrap/KiP_Survival_rates.txt", header = TRUE)
KiP_ASR_boot <- read.table(file = "output/Bootstrap/KiP_ASR.txt", header = TRUE)

KiP_surv_boot$iter <- as.factor(KiP_surv_boot$iter)

KiP_avg_rates <- 
  KiP_surv_boot %>%
  group_by(Sex_Age) %>%
  summarise(Avg = mean(estimate))

KiP_avg_rates <- as.data.frame(KiP_avg_rates)

# Define KiP vital rates estimated from mark-recapture analysis:
KiP_VR <- list(F_Juv_survl = KiP_avg_rates[3,2],
               F_Adt_survl = KiP_avg_rates[1,2],
               M_Juv_survl = KiP_avg_rates[6,2],
               M_Adt_survl = KiP_avg_rates[4,2],
               # Define h (harem size, h < 1 is polyandry) and k (clutch size)
               h = 1,
               k = 2,
               # Define primary sex ratio (assumed to be 0.5)
               HSR = 0.5)

# KiP matrix:
KiP_matrix <- plover_matrix(KiP_VR, chick_surv = FALSE)

# Determine the ASR at the stable stage distribution
KiP_ASR_h_1 <- freq_dep_SSD_ASR(A = KiP_matrix, h = 1, k = 2)

# Lower-level vital rate sensitivity analysis
KiP_LLSA <- lower_level_sens_analysis(freq_dep_ASR = KiP_ASR_h_1,
                                      VR_list = KiP_VR, chick_surv = FALSE)

# Calculate vital rate sensitivities and elasticities
KiP_VR_Sens_Elas <- vitalsens_ASR(elements = KiP_LLSA$Element_expression,
                                  VR_list = KiP_VR, freq_dep_ASR = KiP_ASR_h_1,
                                  chick_surv = FALSE,
                                  pop_name = "Kittlitz's")

# stack all the population specific sensitivity results
VR_Sens_Elas_males <- rbind(Ceuta_VR_Sens_Elas$males,
                            Tuzla_VR_Sens_Elas$males,
                            MP_VR_Sens_Elas$males,
                            WfP_VR_Sens_Elas$males,
                            Maio_VR_Sens_Elas$males,
                            KiP_VR_Sens_Elas$males)

# define correct levels of population and stage variables to arrange display
VR_Sens_Elas_males$population <- 
  factor(VR_Sens_Elas_males$population ,
         levels = c("Snowy",
                    "Kentish-Tuzla",
                    "Madagascar",
                    "Kentish-Maio",
                    "White-fronted",
                    "Kittlitz's"))

VR_Sens_Elas_males$VR <- 
  factor(VR_Sens_Elas_males$VR ,
         levels = c("Juvenile",
                    "Adult"))

# Custom color palette for the plotting of Juvenile and Adult stats
cbPalette <- c("#BDBDBD", "#737373")

# plot the comparative LTRE results
LTRE_bar_plot <- 
  ggplot() +  
  theme_bw() +
  geom_bar(data = VR_Sens_Elas_males[which(VR_Sens_Elas_males$variable == "LTRE" &
                                             (VR_Sens_Elas_males$VR == "Juvenile" |
                                                VR_Sens_Elas_males$VR == "Adult"| 
                                                VR_Sens_Elas_males$VR == "Fledgling")),], 
           aes(x=VR, y=value_trans, fill = VR), stat = "identity", alpha = 0.8) +
  facet_grid(. ~ population) +
  theme(text=element_text(family="Menlo"),
        legend.position="none",
        legend.position = c(0, 1), 
        legend.justification = c(0, 1),
        legend.text=element_text(size=11),
        legend.title=element_blank(),
        legend.key.height=unit(0.8,"line"),
        legend.key.width=unit(0.8,"line"),
        legend.background = element_rect(fill=NA),
        axis.title.x = element_text(size=12, vjust=-0.1),
        axis.text.x  = element_text(size=10, angle = 45, hjust = 1), 
        axis.title.y = element_text(size=12, margin = margin(0, 13, 0, 0)),
        axis.text.y  = element_text(size=10),#, hjust = 1),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "grey"),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        plot.margin = unit(c(0.5,1.75,0.5,0.5), "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size=11),
        panel.margin = unit(0.75, "lines")) +
  ylab("Contribution to adult sex ratio") +
  xlab("Apparent survival (\u03D5)") +
  scale_fill_manual(values = cbPalette) +
  scale_y_continuous(limits=c(0,0.6))
LTRE_bar_plot

LTRE_bar_plot_blank <- 
  ggplot(data = VR_Sens_Elas_males[which(VR_Sens_Elas_males$variable == "LTRE" &
                                           (VR_Sens_Elas_males$VR == "Juvenile" |
                                              VR_Sens_Elas_males$VR == "Adult"| 
                                              VR_Sens_Elas_males$VR == "Fledgling")),], 
         aes(x=VR, y=value_trans, fill = VR)) +  
  theme_bw() +
  geom_blank() +
  facet_grid(. ~ population) +
  theme(text=element_text(family="Menlo"),
        legend.position="none",
        legend.position = c(0, 1), 
        legend.justification = c(0, 1),
        legend.text=element_text(size=11),
        legend.title=element_blank(),
        legend.key.height=unit(0.8,"line"),
        legend.key.width=unit(0.8,"line"),
        legend.background = element_rect(fill=NA),
        axis.title.x = element_text(size=12, vjust=-0.1),
        axis.text.x  = element_text(size=10, angle = 45, hjust = 1), 
        axis.title.y = element_text(size=12, margin = margin(0, 13, 0, 0)),
        axis.text.y  = element_text(size=10),#, hjust = 1),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "grey"),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        plot.margin = unit(c(0.5,1.75,0.5,0.5), "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size=11),
        panel.margin = unit(0.75, "lines")) +
  ylab("Contribution to adult sex ratio") +
  xlab("Apparent survival (\u03D5)") +
  scale_fill_manual(values = cbPalette) +
  scale_y_continuous(limits=c(0,0.6))
LTRE_bar_plot_blank

# save LTRE plot
ggsave(LTRE_bar_plot, 
       filename = "LTRE_bar_plot.jpg", 
       path = "figs/",
       width = 10,
       height = 4.5, units = "in",
       dpi = 300)

ggsave(LTRE_bar_plot_blank, 
       filename = "LTRE_bar_plot_blank.jpg", 
       path = "figs/",
       width = 10,
       height = 4.5, units = "in",
       dpi = 300)

###############################################################################
# Comparative plot of sex- and stage specific elasticities vs. ASR
# ----------------------------------------------------------------
male_sens <- rbind(Ceuta_VR_Sens_Elas$males, Tuzla_VR_Sens_Elas$males, 
                   MP_VR_Sens_Elas$males, 
                   WfP_VR_Sens_Elas$males,
                   Maio_VR_Sens_Elas$males, KiP_VR_Sens_Elas$males)
female_sens <- rbind(Ceuta_VR_Sens_Elas$females, Tuzla_VR_Sens_Elas$females, 
                     MP_VR_Sens_Elas$females, 
                     WfP_VR_Sens_Elas$females,
                     Maio_VR_Sens_Elas$females, KiP_VR_Sens_Elas$females)

plover_sens <- rbind(male_sens, female_sens)
plover_sens$population <- as.factor(plover_sens$population)
plover_sens$population <- factor(plover_sens$population, levels = c("Snowy",
                                                        "Kentish-Tuzla",
                                                        "Madagascar",
                                                        "White-fronted",
                                                        "Kentish-Maio",
                                                        "Kittlitz's"))
plover_sens$VR <- factor(plover_sens$VR, levels = c("Adult",
                                                    "Juvenile",
                                                    " ",
                                                    "No. breeding adults",
                                                    "Mating System",
                                                    "Clutch size",
                                                    "Hatching sex ratio"))

data <- unname(c(Ceuta_ASR_h_1$ASR,
                 Tuzla_ASR_h_1$ASR,
                 MP_ASR_h_1$ASR,
                 WfP_ASR_h_1$ASR,
                 Maio_ASR_h_1$ASR,
                 KiP_ASR_h_1$ASR))
population <- c("Snowy",
          "Kentish-Tuzla",
          "Madagascar",
          "White-fronted",
          "Kentish-Maio",
          "Kittlitz's")
comp_ASR <- data.frame(population, data)
comp_ASR$trans <- comp_ASR$data-0.5
comp_ASR$population <- 
  factor(comp_ASR$population ,
         levels = c("Snowy",
                    "Kentish-Tuzla",
                    "Madagascar",
                    "White-fronted",
                    "Kentish-Maio",
                    "Kittlitz's"))

plover_sens$value <- as.numeric(plover_sens$value)
plover_sens <- full_join(plover_sens, comp_ASR, by = "population")
cbPalette <- c("#C0504D", "#558ED5")
plover_sens$abs_trans <- abs(plover_sens$trans)

mf <-lm(data ~ value, filter(plover_sens, variable == "LTRE" & VR == "Juvenile" & Sex == "Male"))
mfeq <- as.character(as.expression(substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                              list(a = format(coef(mf)[1], digits = 2), 
                                                   b = format(coef(mf)[2], digits = 2), 
                                                   r2 = format(summary(mf)$r.squared, digits = 3)))))
ma <-lm(data ~ value, filter(plover_sens, variable == "LTRE" & VR == "Juvenile" & Sex == "Male"))
maeq <- as.character(as.expression(substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                              list(a = format(coef(ma)[1], digits = 2), 
                                                   b = format(coef(ma)[2], digits = 2), 
                                                   r2 = format(summary(ma)$r.squared, digits = 3)))))
ff <-lm(data ~ value, filter(plover_sens, variable == "LTRE" & VR == "Juvenile" & Sex == "Female"))
ffeq <- as.character(as.expression(substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                              list(a = format(coef(ff)[1], digits = 2), 
                                                   b = format(coef(ff)[2], digits = 2), 
                                                   r2 = format(summary(ff)$r.squared, digits = 3)))))
fa <-lm(data ~ value, filter(plover_sens, variable == "LTRE" & VR == "Juvenile" & Sex == "Female"))
faeq <- as.character(as.expression(substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                              list(a = format(coef(fa)[1], digits = 2), 
                                                   b = format(coef(fa)[2], digits = 2), 
                                                   r2 = format(summary(fa)$r.squared, digits = 3)))))
lm_eqn <- data.frame(rbind(mfeq, maeq, ffeq, faeq))
lm_eqn$Sex <- c("Male", "Male", "Female", "Female")
lm_eqn$VR <- c("Juvenile", "Adult")
row.names(lm_eqn) <- NULL
colnames(lm_eqn) <- c("V1", "Sex", "VR")

Comp_Elasticity_plot <- 
  ggplot(filter(plover_sens, variable == "LTRE" & 
                  (VR == "Juvenile" | VR == "Adult")),
         aes(x = value_trans, y = data)) +
  geom_point(size = I(4), alpha = I(0.7)) + 
  geom_smooth(method = lm, aes(colour = Sex, fill = Sex), size=2) +
  theme_bw() +
  theme(text=element_text(family="Menlo"),
        legend.position="none",
        legend.position = c(0, 1), 
        legend.justification = c(0, 1),
        legend.text=element_text(size=11),
        legend.title=element_blank(),
        legend.key.height=unit(0.8,"line"),
        legend.key.width=unit(0.8,"line"),
        legend.background = element_rect(fill=NA),
        axis.title.x = element_text(size = 14),
        axis.text.x  = element_text(size = 13), 
        axis.title.y = element_text(size = 14),
        axis.text.y  = element_text(size = 13), 
        panel.grid.major = element_blank(),
        strip.text.x = element_text(size=14, face = "bold"),
        strip.text.y = element_text(size=14, face = "bold")) +
  facet_grid(VR ~ Sex, scales = "free") +
  ylab("Adult sex ratio (proportion \u2642 ± 95% CI)") +
  xlab("Vital rate elasticity") +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  geom_text(data = lm_eqn, aes(x = -Inf, y = 1,label = V1, family="Menlo"), 
            parse = TRUE, inherit.aes=FALSE, hjust = -0.05)

# Export 
ggsave(Comp_Elasticity_plot, 
       filename = "Comparative_Elasticity_ASR_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 12,
       height = 6, units = "in",
       dpi = 300)


LTRE_ASR_stage_regression <- 
  ggplot(filter(plover_sens, variable == "LTRE" & 
                  (VR == "Juvenile" | VR == "Adult") &
                  Sex == "Male"),
         aes(y = value_trans, x = data)) +
  geom_smooth(method = lm, formula = y ~ poly(x, 2), aes(colour = VR, fill = VR), size=2) +
  geom_point(size = I(4), alpha = I(0.7)) +
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
        axis.title.x = element_text(size = 14),
        axis.text.x  = element_text(size = 13), 
        axis.title.y = element_text(size = 14),
        axis.text.y  = element_text(size = 13), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=14, face = "bold"),
        strip.text.y = element_text(size=14, face = "bold")) +
  facet_grid(VR ~ .) +
  ylab("Contribution (± 95% CI)") +
  xlab("Adult sex ratio (proportion \u2642)")
Mating_system_ASR_sex_regression

# Export 
ggsave(Mating_system_ASR_sex_regression, 
       filename = "Mating_system_ASR_sex_regression.jpg", 
       path = "Figures/",
       width = 4,
       height = 6, units = "in",
       dpi = 300)