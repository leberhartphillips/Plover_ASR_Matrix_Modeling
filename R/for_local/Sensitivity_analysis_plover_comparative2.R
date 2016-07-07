###############################################################################
########## Deterministic Projection Matrix Model for Plover ASR  with #########
################## frequency dependent fecundity function #####################
###############################################################################
############################# Luke Eberhart-Phillips ##########################
###############################################################################

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
library(xkcd)

# Find fonts from computer that are candara or Candara
font_import(pattern="[C/c]andara|[X/x]kcd", prompt = FALSE) 
fonts()
fonttable()
loadfonts(device = "win") # load these into R

# Check that the most recent dev version of ggplot2 is loaded:
# install.packages("devtools")
# library(devtools)
# devtools::install_github("hadley/ggplot2")

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
# HSR ---> the primary sex ratio (default is 0.5)
# iterations ---> the number of iterations to simulate (default is 20)

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
  function (elements, VR_list, freq_dep_ASR, chick_surv = TRUE) 
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
      y_melt$VR <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"Juv"), "Fledgling",
                                    ifelse(str_detect(y_melt$Vital_rate,"Adt"), "Adult",
                                           ifelse(str_detect(y_melt$Vital_rate,"HSR"), "Hatching sex ratio", 
                                                  ifelse(str_detect(y_melt$Vital_rate,"h"), "Mating System",
                                                         ifelse(str_detect(y_melt$Vital_rate,"k"), "Clutch size", "No. breeding adults"))))))
      y_melt$VR <- factor(y_melt$VR, levels = c("Adult",
                                                "Fledgling",
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
    y_melt$Sex <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"F_"), "Female", 
                                   ifelse(str_detect(y_melt$Vital_rate,"M_"),"Male", "Other")))
    y_melt$Sex <- 
      factor(y_melt$Sex,
             levels = c("Female","Male", "Other"))
    # y_melt$value_trans <- ifelse(y_melt$Sex == "Female", abs(y_melt$value)*-1, y_melt$value)
    # males <- subset(y_melt, (Sex == "Male" | Sex == "Other"))
    # females <- subset(y_melt, Sex == "Female" )
    # results <- list(males = males,
    #                 females = females)
    # row.names(results$males) <- NULL
    # row.names(results$females) <- NULL
    # results$males[nrow(results$males)+1,] <- c("Dummy","Elasticity",0," ","Other",0)
    # results$females[nrow(results$females)+1,] <- c("Dummy","Elasticity",0," ","Other",0)
    # results$males[nrow(results$males)+1,] <- c("Dummy","Sensitivity",0," ","Other",0)
    # results$females[nrow(results$females)+1,] <- c("Dummy","Sensitivity",0," ","Other",0)
    # results$males$value_trans <- as.numeric(results$males$value_trans)
    # results$females$value_trans <- as.numeric(results$females$value_trans)
    # if(!chick_surv)
    # {
    #   results$males[nrow(results$males)+1,] <- c("Dummy","Elasticity",0,"Chick","Male",0)
    #   results$females[nrow(results$females)+1,] <- c("Dummy","Elasticity",0,"Chick","Female",0)
    #   results$males[nrow(results$males)+1,] <- c("Dummy","Sensitivity",0,"Chick","Male",0)
    #   results$females[nrow(results$females)+1,] <- c("Dummy","Sensitivity",0,"Chick","Female",0)
    #   results$males$value_trans <- as.numeric(results$males$value_trans)
    #   results$females$value_trans <- as.numeric(results$females$value_trans)
    # }
    # results
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

Elasticity_line_plot <- 
  function(lower_level_sens_analysis_result)
  {
    Elasticity_line_plot <- ggplot(lower_level_sens_analysis_result$Elasticities, aes(x = Perturbation, y = Elasticity, group = Vitalrate)) +  
      theme_bw() +
      geom_line(size = 1, aes(colour = Vitalrate)) +
      theme(text=element_text(family="Candara"),
            legend.position = c(0, 1), 
            legend.justification = c(0, 1),
            legend.text=element_text(size=11),
            legend.title=element_text(size=12),
            legend.key.height=unit(0.8,"line"),
            legend.key.width=unit(0.8,"line"),
            legend.background = element_rect(fill=NA),
            axis.title.x = element_text(size=12, vjust=-0.1),
            axis.text.x  = element_text(size=11), 
            axis.title.y = element_text(size=12, vjust=1.2),
            axis.text.y  = element_text(size=11), 
            panel.grid.major = element_blank()) +
      xlab("Proportion of current vital rate") +
      ylab("Proportion of current ASR") +
      scale_colour_manual(values=cbPalette, name="Vitalrate") +
      annotate("text", x = -Inf, y = Inf, vjust = 1, hjust = 0, label = "(b)", 
               size = 4, family="Candara")
    Elasticity_line_plot
  }

Sensitivity_line_plot <- 
  function(lower_level_sens_analysis_result)
  {
    Sensitivity_line_plot <- ggplot(lower_level_sens_analysis_result$Sensitivities, aes(x = Perturbation, y = Sensitivity, group = Vitalrate)) +  
      theme_bw() +
      geom_line(size = 1, aes(colour = Vitalrate)) +
      theme(text=element_text(family="Candara"),
            legend.position = c(0, 1), 
            legend.justification = c(0, 1),
            legend.text=element_text(size=11),
            legend.title=element_text(size=12),
            legend.key.height=unit(0.8,"line"),
            legend.key.width=unit(0.8,"line"),
            legend.background = element_rect(fill=NA),
            axis.title.x = element_text(size=12, vjust=-0.1),
            axis.text.x  = element_text(size=11), 
            axis.title.y = element_text(size=12, vjust=1.2),
            axis.text.y  = element_text(size=11), 
            panel.grid.major = element_blank()) +
      xlab("Value of vital rate") +
      ylab(expression(paste("Adult sex ratio"))) +
      scale_colour_manual(values=cbPalette, name="Vitalrate") +
      annotate("text", x = -Inf, y = Inf, vjust = 1, hjust = 0, label = "(a)", 
               size = 4, family="Candara")
    Sensitivity_line_plot
  }

Elasticity_bar_plot <- 
  function(vitalsens_ASR_result, chick_surv = TRUE)
  {
    if(!chick_surv)
    {
      Elasticity_bar_plot <- 
        ggplot() +  
        theme_bw() +
        geom_bar(data = vitalsens_ASR_result$males[which(vitalsens_ASR_result$males$variable == "Elasticity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity") +
        geom_bar(data = vitalsens_ASR_result$females[which(vitalsens_ASR_result$females$variable == "Elasticity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity") +
        coord_flip() +
        theme(text=element_text(family="Candara"),
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
              axis.title.y = element_blank(),
              axis.text.y  = element_text(size=10),#, angle=90, hjust = 0.5, vjust = 0.5),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_rect(linetype = "solid", colour = "grey")) +
        ylab("Elasticity of adult sex ratio") +
        scale_fill_brewer(palette = "Set1") +
        scale_y_continuous(limits = c(-1.8, 1.8)) +
        annotate("text", x = 4, y = -1, vjust = 1.5, hjust = 0.5, label = "\u2640", 
                 size = 7, family="Candara") +
        annotate("text", x = 4, y = 1, vjust = 1.5, hjust = 0.5, label = "\u2642", 
                 size = 7, family="Candara") +
        annotate("text", x = 3, y = -1, vjust = 0.5, hjust = 0.5, label = 'italic("NA")', 
                 size = 3, family="Candara", parse = TRUE) +
        annotate("text", x = 3, y = 1, vjust = 0.5, hjust = 0.5, label = 'italic("NA")', 
                 size = 3, family="Candara", parse = TRUE) +
        annotate("rect", xmin = -Inf, xmax = 4.25, ymin = -Inf, ymax = Inf,
                 alpha = .1)
    }
    else
    {
      Elasticity_bar_plot <- 
        ggplot() +  
        theme_bw() +
        geom_bar(data = vitalsens_ASR_result$males[which(vitalsens_ASR_result$males$variable == "Elasticity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity") +
        geom_bar(data = vitalsens_ASR_result$females[which(vitalsens_ASR_result$females$variable == "Elasticity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity") +
        coord_flip() +
        theme_xkcd() +
        theme(text=element_text(family="Candara"),
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
              axis.title.y = element_blank(),
              axis.text.y  = element_text(size=10),#, angle=90, hjust = 0.5, vjust = 0.5),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_rect(linetype = "solid", colour = "grey")) +
        ylab("Elasticity of adult sex ratio") +
        scale_fill_brewer(palette = "Set1") +
        scale_y_continuous(limits = c(-1.8, 1.8)) +
        annotate("text", x = 4, y = -1, vjust = 1.5, hjust = 0.5, label = "\u2640", 
                 size = 7, family="Candara") +
        annotate("text", x = 4, y = 1, vjust = 1.5, hjust = 0.5, label = "\u2642", 
                 size = 7, family="Candara") +
        annotate("rect", xmin = -Inf, xmax = 4.25, ymin = -Inf, ymax = Inf,
                 alpha = .1)
    }
    Elasticity_bar_plot
  }

Sensitivity_bar_plot <- 
  function(vitalsens_ASR_result, chick_surv = TRUE)
  {
    if(!chick_surv)
    {
      Sensitivity_bar_plot <- 
        ggplot() +  
        theme_bw() +
        geom_bar(data = vitalsens_ASR_result$males[which(vitalsens_ASR_result$males$variable == "Sensitivity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity") +
        geom_bar(data = vitalsens_ASR_result$females[which(vitalsens_ASR_result$females$variable == "Sensitivity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity") +
        coord_flip() +
        theme(text=element_text(family="Candara"),
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
              axis.title.y = element_blank(),
              axis.text.y  = element_text(size=10, hjust = 1),#, angle=90, hjust = 0.5, vjust = 0.5),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_rect(linetype = "solid", colour = "grey")) +
        ylab("Sensitivity of adult sex ratio") +
        scale_fill_brewer(palette = "Set1") +
        scale_y_continuous(limits = c(-1, 1)) +
        annotate("text", x = 4, y = -0.65, vjust = 1.5, hjust = 0.5, label = "\u2640", 
                 size = 7, family="Candara") +
        annotate("text", x = 4, y = 0.65, vjust = 1.5, hjust = 0.5, label = "\u2642", 
                 size = 7, family="Candara") +
        annotate("text", x = 3, y = -0.65, vjust = 0.5, hjust = 0.5, label = 'italic("NA")', 
                 size = 3, family="Candara", parse = TRUE) +
        annotate("text", x = 3, y = 0.65, vjust = 0.5, hjust = 0.5, label = 'italic("NA")', 
                 size = 3, family="Candara", parse = TRUE) +
        annotate("rect", xmin = -Inf, xmax = 4.25, ymin = -Inf, ymax = Inf,
                 alpha = .1)
    }
    else
    {
      Sensitivity_bar_plot <- 
        ggplot() +  
        theme_bw() +
        geom_bar(data = vitalsens_ASR_result$males[which(vitalsens_ASR_result$males$variable == "Sensitivity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity") +
        geom_bar(data = vitalsens_ASR_result$females[which(vitalsens_ASR_result$females$variable == "Sensitivity"),], aes(x=VR, y=value_trans, fill=Sex),stat = "identity") +
        coord_flip() +
        theme(text=element_text(family="Candara"),
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
              axis.title.y = element_blank(),
              axis.text.y  = element_text(size=10, hjust = 1),#, angle=90, hjust = 0.5, vjust = 0.5),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_rect(linetype = "solid", colour = "grey")) +
        ylab("Sensitivity of adult sex ratio") +
        scale_fill_brewer(palette = "Set1") +
        scale_y_continuous(limits = c(-1, 1)) +
        annotate("text", x = 4, y = -0.65, vjust = 1.5, hjust = 0.5, label = "\u2640", 
                 size = 7, family="Candara") +
        annotate("text", x = 4, y = 0.65, vjust = 1.5, hjust = 0.5, label = "\u2642", 
                 size = 7, family="Candara") +
        annotate("rect", xmin = -Inf, xmax = 4.25, ymin = -Inf, ymax = Inf,
                 alpha = .1)
    }
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

###############################################################################
# CEUTA SNOWY PLOVER 
#--------------------

# Define Ceuta vital rates estimated from mark-recapture analysis:
Ceuta_VR <- list(F_Chk_survl = 0.4101527,
                 F_Juv_survl = 0.1541992,
                 F_Adt_survl = 0.6607083,
                 M_Chk_survl = 0.4640493,
                 M_Juv_survl = 0.2282648,
                 M_Adt_survl = 0.7111575,
                 # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                 h = 0.3,
                 k = 3,
                 # Define primary sex ratio (assumed to be 0.5)
                 HSR = 0.5)

# Ceuta matrix:
Ceuta_matrix <- plover_matrix(Ceuta_VR)

# Determine the ASR at the stable stage distribution (assume h is 0.3, which 
# is polyandry with females having on average 3 males)
Ceuta_ASR_h_0.3 <- freq_dep_SSD_ASR(A = Ceuta_matrix, h = 0.3, k = 3)

# Lower-level vital rate sensitivity analysis
Ceuta_LLSA <- lower_level_sens_analysis(freq_dep_ASR = Ceuta_ASR_h_0.3, VR_list = Ceuta_VR)

# Calculate vital rate sensitivities and elasticities
Ceuta_VR_Sens_Elas <- vitalsens_ASR(elements = Ceuta_LLSA$Element_expression,
                                        VR_list = Ceuta_VR, freq_dep_ASR = Ceuta_ASR_h_0.3)

# define color palette
cbPalette <- c(brewer.pal(9,"Reds")[c(5,7,9)], brewer.pal(9,"Blues")[c(5,7,9)],
               brewer.pal(8, "Dark2")[c(1,8,5,2,7)])

# Plot linear sex- and stage-specific elasticities
Ceuta_Elasticity_line_plot <- Elasticity_line_plot(Ceuta_LLSA)

# Plot linear sex- and stage-specific sensitivities
Ceuta_Sensitivity_line_plot <- Sensitivity_line_plot(Ceuta_LLSA)

# Plot overall sex- and stage-specific elasticities
Ceuta_Elasticity_bar_plot <- Elasticity_bar_plot(Ceuta_VR_Sens_Elas, chick_surv = TRUE)

# Export elasticity bar plot
ggsave(Ceuta_Elasticity_bar_plot, 
       filename = "Ceuta_Elasticity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

# Plot overall sex- and stage-specific sensitivities
Ceuta_Sensitivity_bar_plot <- Sensitivity_bar_plot(Ceuta_VR_Sens_Elas, chick_surv = TRUE)

# Export sensitivity bar plot
ggsave(Ceuta_Sensitivity_bar_plot, 
       filename = "Ceuta_Sensitivity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

# # Save arrangement in the plotting device as a grob
# Ceuta_Sen_Ela_combo_plot <- arrangeGrob(Ceuta_Sensitivity_bar_plot,
#                                         Ceuta_Elasticity_bar_plot,
#                                         ncol=1)
# 
# # export the grob to the directory
# ggsave(Sen_Ela_combo_plot, 
#        filename = "Ceuta_Elasticities_and_Sensitivities.jpg", 
#        path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
#        width = 4,
#        height = 8, units = "in",
#        dpi = 300)

###############################################################################
# TUZLA KENTISH PLOVER
# --------------------
# Define Tuzla vital rates estimated from mark-recapture analysis:
Tuzla_VR <- list(
  F_Juv_survl = 0.3796332,
  F_Adt_survl = 0.7583189,
  M_Juv_survl = 0.5778207,
  M_Adt_survl = 0.7431910,
  # Define h (harem size, h < 1 is polyandry) and k (clutch size)
  h = 0.5,
  k = 3,
  # Define primary sex ratio (assumed to be 0.5)
  HSR = 0.5)

# Tuzla matrix:
Tuzla_matrix <- plover_matrix(Tuzla_VR, chick_surv = FALSE)

# Determine the ASR at the stable stage distribution (assume h is 0.5, which 
# is polyandry with females having on average 2 males)
Tuzla_ASR_h_0.5 <- freq_dep_SSD_ASR(A = Tuzla_matrix, h = 0.5, k = 3)
Tuzla_ASR_h_0.3 <- freq_dep_SSD_ASR(A = Tuzla_matrix, h = 0.3, k = 3)

# Lower-level vital rate sensitivity analysis
Tuzla_LLSA <- lower_level_sens_analysis(freq_dep_ASR = Tuzla_ASR_h_0.3,
                                        VR_list = Tuzla_VR, chick_surv = FALSE)

# Calculate vital rate sensitivities and elasticities
Tuzla_VR_Sens_Elas <- vitalsens_ASR(elements = Tuzla_LLSA$Element_expression,
                                    VR_list = Tuzla_VR, freq_dep_ASR = Tuzla_ASR_h_0.3,
                                    chick_surv = FALSE)

# define color palette
cbPalette <- c(brewer.pal(9,"Reds")[c(7,9)], brewer.pal(9,"Blues")[c(7,9)],
               brewer.pal(8, "Dark2")[c(1,8,5,2,7)])

# Plot linear sex- and stage-specific elasticities
Tuzla_Elasticity_line_plot <- Elasticity_line_plot(Tuzla_LLSA)

# Plot linear sex- and stage-specific sensitivities
Tuzla_Sensitivity_line_plot <- Sensitivity_line_plot(Tuzla_LLSA)

# Plot overall sex- and stage-specific elasticities
Tuzla_Elasticity_bar_plot <- Elasticity_bar_plot(Tuzla_VR_Sens_Elas, chick_surv = FALSE)

# Export elasticity bar plot
ggsave(Tuzla_Elasticity_bar_plot, 
       filename = "Tuzla_Elasticity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

# Plot overall sex- and stage-specific sensitivities
Tuzla_Sensitivity_bar_plot <- Sensitivity_bar_plot(Tuzla_VR_Sens_Elas, chick_surv = FALSE)

# Export sensitivity bar plot
ggsave(Tuzla_Sensitivity_bar_plot, 
       filename = "Tuzla_Sensitivity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

# Plot overall sex- and stage-specific LTRE
Tuzla_LTRE_bar_plot <- LTRE_bar_plot(Tuzla_VR_Sens_Elas)

###############################################################################
# MADAGASCAR PLOVER
# --------------------
# Define Madagascar plover vital rates estimated from mark-recapture analysis:
Madagascar_VR <- list(
  F_Juv_survl = 0.0643717,
  F_Adt_survl = 0.6490020,
  M_Juv_survl = 0.1391483,
  M_Adt_survl = 0.5962389,
  # Define h (harem size, h<1 is polyandry) and k (clutch size)
  h = 1,
  k = 2,
  # Define primary sex ratio (assumed to be 0.5)
  HSR = 0.5)

# Madagascar matrix:
Madagascar_matrix <- plover_matrix(Madagascar_VR, chick_surv = FALSE)

# Determine the ASR at the stable stage distribution (assume h is 1, which 
# is monogamy with females having on average 1 males, and k is 2 (i.e., 2eggs))
Madagascar_ASR_h_1 <- freq_dep_SSD_ASR(A = Madagascar_matrix, h = 1, k = 2)

# Lower-level vital rate sensitivity analysis
Madagascar_LLSA <- lower_level_sens_analysis(freq_dep_ASR = Madagascar_ASR_h_1,
                                             VR_list = Madagascar_VR, 
                                             chick_surv = FALSE)

# Calculate vital rate sensitivities and elasticities
Madagascar_VR_Sens_Elas <- vitalsens_ASR(elements = Madagascar_LLSA$Element_expression,
                                         VR_list = Madagascar_VR, freq_dep_ASR = Madagascar_ASR_h_1,
                                         chick_surv = FALSE)

# define color palette
cbPalette <- c(brewer.pal(9,"Reds")[c(7,9)], brewer.pal(9,"Blues")[c(7,9)],
               brewer.pal(8, "Dark2")[c(1,8,5,2,7)])

# Plot linear sex- and stage-specific elasticities
Madagascar_Elasticity_line_plot <- Elasticity_line_plot(Madagascar_LLSA)

# Plot linear sex- and stage-specific sensitivities
Madagascar_Sensitivity_line_plot <- Sensitivity_line_plot(Madagascar_LLSA)

# Plot overall sex- and stage-specific elasticities
Madagascar_Elasticity_bar_plot <- Elasticity_bar_plot(Madagascar_VR_Sens_Elas, chick_surv = FALSE)

# Export elasticity bar plot
ggsave(Madagascar_Elasticity_bar_plot, 
       filename = "Madagascar_Elasticity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

# Plot overall sex- and stage-specific sensitivities
Madagascar_Sensitivity_bar_plot <- Sensitivity_bar_plot(Madagascar_VR_Sens_Elas, chick_surv = FALSE)

# Export sensitivity bar plot
ggsave(Madagascar_Sensitivity_bar_plot, 
       filename = "Madagascar_Sensitivity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

###############################################################################
# WHITE-FRONTED PLOVER
# --------------------
# Define White_fronted plover vital rates estimated from mark-recapture analysis:
White_fronted_VR <- list(
  F_Juv_survl = 0.5699494,
  F_Adt_survl = 0.8106450,
  M_Juv_survl = 0.4559940,
  M_Adt_survl = 0.8784390,
  # Define h (harem size, h<1 is polyandry) and k (clutch size)
  h = 1,
  k = 2,
  # Define primary sex ratio (assumed to be 0.5)
  HSR = 0.5)

# White_fronted matrix:
White_fronted_matrix <- plover_matrix(White_fronted_VR, chick_surv = FALSE)

# Determine the ASR at the stable stage distribution (assume h is 1, which 
# is monogamy with females having on average 1 males, and k is 2 (i.e., 2eggs))
White_fronted_ASR_h_1 <- freq_dep_SSD_ASR(A = White_fronted_matrix, h = 1, k = 2)

# Lower-level vital rate sensitivity analysis
White_fronted_LLSA <- lower_level_sens_analysis(freq_dep_ASR = White_fronted_ASR_h_1,
                                                VR_list = White_fronted_VR, 
                                                chick_surv = FALSE)

# Calculate vital rate sensitivities and elasticities
White_fronted_VR_Sens_Elas <- vitalsens_ASR(elements = White_fronted_LLSA$Element_expression,
                                            VR_list = White_fronted_VR, freq_dep_ASR = White_fronted_ASR_h_1,
                                            chick_surv = FALSE)

# define color palette
cbPalette <- c(brewer.pal(9,"Reds")[c(7,9)], brewer.pal(9,"Blues")[c(7,9)],
               brewer.pal(8, "Dark2")[c(1,8,5,2,7)])

# Plot linear sex- and stage-specific elasticities
White_fronted_Elasticity_line_plot <- Elasticity_line_plot(White_fronted_LLSA)

# Plot linear sex- and stage-specific sensitivities
White_fronted_Sensitivity_line_plot <- Sensitivity_line_plot(White_fronted_LLSA)

# Plot overall sex- and stage-specific elasticities
White_fronted_Elasticity_bar_plot <- Elasticity_bar_plot(White_fronted_VR_Sens_Elas, chick_surv = FALSE)

# Export elasticity bar plot
ggsave(White_fronted_Elasticity_bar_plot, 
       filename = "White_fronted_Elasticity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

# Plot overall sex- and stage-specific sensitivities
White_fronted_Sensitivity_bar_plot <- Sensitivity_bar_plot(White_fronted_VR_Sens_Elas, chick_surv = FALSE)

# Export sensitivity bar plot
ggsave(White_fronted_Sensitivity_bar_plot, 
       filename = "White_fronted_Sensitivity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

###############################################################################
# MAIO KENTISH PLOVER
# --------------------
# Define Maio plover vital rates estimated from mark-recapture analysis:
Maio_VR <- list(
  F_Juv_survl = 0.2172390,
  F_Adt_survl = 0.8315268,
  M_Juv_survl = 0.1949260,
  M_Adt_survl = 0.8175885,
  # Define h (harem size, h<1 is polyandry) and k (clutch size)
  h = 1,
  k = 3,
  # Define primary sex ratio (assumed to be 0.5)
  HSR = 0.5)

# Maio matrix:
Maio_matrix <- plover_matrix(Maio_VR, chick_surv = FALSE)

# Determine the ASR at the stable stage distribution (assume h is 1, which 
# is monogamy with females having on average 1 males, and k is 2 (i.e., 2eggs))
Maio_ASR_h_1 <- freq_dep_SSD_ASR(A = Maio_matrix, h = 1, k = 3)

# Lower-level vital rate sensitivity analysis
Maio_LLSA <- lower_level_sens_analysis(freq_dep_ASR = Maio_ASR_h_1,
                                       VR_list = Maio_VR, 
                                       chick_surv = FALSE)

# Calculate vital rate sensitivities and elasticities
Maio_VR_Sens_Elas <- vitalsens_ASR(elements = Maio_LLSA$Element_expression,
                                   VR_list = Maio_VR, freq_dep_ASR = Maio_ASR_h_1,
                                   chick_surv = FALSE)

# define color palette
cbPalette <- c(brewer.pal(9,"Reds")[c(7,9)], brewer.pal(9,"Blues")[c(7,9)],
               brewer.pal(8, "Dark2")[c(1,8,5,2,7)])

# Plot linear sex- and stage-specific elasticities
Maio_Elasticity_line_plot <- Elasticity_line_plot(Maio_LLSA)

# Plot linear sex- and stage-specific sensitivities
Maio_Sensitivity_line_plot <- Sensitivity_line_plot(Maio_LLSA)

# Plot overall sex- and stage-specific elasticities
Maio_Elasticity_bar_plot <- Elasticity_bar_plot(Maio_VR_Sens_Elas, chick_surv = FALSE)

# Export elasticity bar plot
ggsave(Maio_Elasticity_bar_plot, 
       filename = "Maio_Elasticity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

# Plot overall sex- and stage-specific sensitivities
Maio_Sensitivity_bar_plot <- Sensitivity_bar_plot(Maio_VR_Sens_Elas, chick_surv = FALSE)

# Export sensitivity bar plot
ggsave(Maio_Sensitivity_bar_plot, 
       filename = "Maio_Sensitivity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

###############################################################################
# KITTLITZ'S PLOVER
# --------------------
# Define Kittlitz plover vital rates estimated from mark-recapture analysis:
Kittlitz_VR <- list(
  F_Juv_survl = 0.5600676,
  F_Adt_survl = 0.6148649,
  M_Juv_survl = 0.1310105,
  M_Adt_survl = 0.7064647,
  # Define h (harem size, h<1 is polyandry) and k (clutch size)
  h = 1,
  k = 2,
  # Define primary sex ratio (assumed to be 0.5)
  HSR = 0.5)

# Kittlitz matrix:
Kittlitz_matrix <- plover_matrix(Kittlitz_VR, chick_surv = FALSE)

# Determine the ASR at the stable stage distribution (assume h is 1, which 
# is monogamy with females having on average 1 males, and k is 2 (i.e., 2eggs))
Kittlitz_ASR_h_1 <- freq_dep_SSD_ASR(A = Kittlitz_matrix, h = 1, k = 3)

# Lower-level vital rate sensitivity analysis
Kittlitz_LLSA <- lower_level_sens_analysis(freq_dep_ASR = Kittlitz_ASR_h_1,
                                           VR_list = Kittlitz_VR, 
                                           chick_surv = FALSE)

# Calculate vital rate sensitivities and elasticities
Kittlitz_VR_Sens_Elas <- vitalsens_ASR(elements = Kittlitz_LLSA$Element_expression,
                                       VR_list = Kittlitz_VR, freq_dep_ASR = Kittlitz_ASR_h_1,
                                       chick_surv = FALSE)

# define color palette
cbPalette <- c(brewer.pal(9,"Reds")[c(7,9)], brewer.pal(9,"Blues")[c(7,9)],
               brewer.pal(8, "Dark2")[c(1,8,5,2,7)])

# Plot linear sex- and stage-specific elasticities
Kittlitz_Elasticity_line_plot <- Elasticity_line_plot(Kittlitz_LLSA)

# Plot linear sex- and stage-specific sensitivities
Kittlitz_Sensitivity_line_plot <- Sensitivity_line_plot(Kittlitz_LLSA)

# Plot overall sex- and stage-specific elasticities
Kittlitz_Elasticity_bar_plot <- Elasticity_bar_plot(Kittlitz_VR_Sens_Elas, chick_surv = FALSE)

# Export elasticity bar plot
ggsave(Kittlitz_Elasticity_bar_plot, 
       filename = "Kittlitz_Elasticity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

# Plot overall sex- and stage-specific sensitivities
Kittlitz_Sensitivity_bar_plot <- Sensitivity_bar_plot(Kittlitz_VR_Sens_Elas, chick_surv = FALSE)

# Export sensitivity bar plot
ggsave(Kittlitz_Sensitivity_bar_plot, 
       filename = "Kittlitz_Sensitivity_bar_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 3,
       height = 4.7, units = "in",
       dpi = 300)

###############################################################################
# Comparative plot of ASR
# ------------------------
data <- unname(c(Ceuta_ASR_h_0.3$ASR_f,
                 Tuzla_ASR_h_0.3$ASR_f,
                 Madagascar_ASR_h_1$ASR_f,
                 White_fronted_ASR_h_1$ASR_f,
                 Maio_ASR_h_1$ASR_f,
                 Kittlitz_ASR_h_1$ASR_f))
pops <- c("Snowy",
          "Kentish (Tuzla)",
          "Madagascar",
          "White-fronted",
          "Kentish (Maio)",
          "Kittlitz's")
comp_ASR <- data.frame(pops, data)
comp_ASR$trans <- comp_ASR$data-0.5
comp_ASR$pops <- 
  factor(comp_ASR$pops ,
         levels = c("Snowy",
                    "Kentish (Tuzla)",
                    "Madagascar",
                    "White-fronted",
                    "Kentish (Maio)",
                    "Kittlitz's"))

comp_ASR_plot <- 
  ggplot(comp_ASR, aes(x = pops, y = trans)) +  
  theme_bw() +
  #theme_xkcd() +
  geom_bar(stat = "identity", position=position_dodge(), alpha = 0.1) +
  theme(text=element_text(family="Candara"),
        legend.position = c(0, 1), 
        legend.justification = c(0, 1),
        legend.text=element_text(size=11),
        legend.title=element_blank(),
        legend.key.height=unit(0.8,"line"),
        legend.key.width=unit(0.8,"line"),
        legend.background = element_rect(fill=NA),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size=11, angle = 45, hjust = 1), 
        axis.title.y = element_text(size=12, vjust=1.2),
        axis.text.y  = element_text(size=11), 
        panel.grid.major = element_blank()) +
  ylab("Adult sex ratio (proportion \u2642)") +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = NULL) +
  scale_x_discrete(breaks = NULL) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=Inf, alpha=0.5,
           fill=brewer.pal(8, "Set1")[c(2)]) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=-Inf, alpha=0.5,
           fill=brewer.pal(8, "Set1")[c(1)]) +
  annotate("text", x = 1:3, y = c(0.24267471, 0.09240868, 0.08616934),
           label = c("0.74", "0.59", "0.58"), vjust = -1, size = 4,
           family="Candara") +
  annotate("text", x = 4:6, y = c(-0.01414658, -0.03949824, -0.26598766),
           label = c("0.48", "0.46", "0.23"), vjust = 1.5, size = 4,
           family="Candara") +
  annotate("text", x = c(-Inf,Inf), y = c(-Inf, Inf),
           label = c("\u2640", "\u2642"), size = 7,
           family="Candara", vjust = c(-1,2), hjust = c(-0.5,1.5))
comp_ASR_plot

# Export comparative ASR_f plot
ggsave(comp_ASR_plot, 
       filename = "Comparative_ASR_plot_empty.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 8,
       height = 3.5, units = "in",
       dpi = 300)
