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

KiP_VR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/KiP_Survival_rates.txt",
                      colClasses = c("factor", "numeric", "factor", "factor"),
                      header = TRUE)

WfP_VR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/WfP_Survival_rates.txt",
                     colClasses = c("factor", "numeric", "factor", "factor"),
                     header = TRUE)

MP_VR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/MP_Survival_rates.txt",
                    colClasses = c("factor", "numeric", "factor", "factor"),
                    header = TRUE) 

Tuzla_VR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/Tuzla_Survival_rates.txt",
                     colClasses = c("factor", "numeric", "factor", "factor"),
                     header = TRUE)

Maio_VR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/Maio_Survival_rates.txt",
                       colClasses = c("factor", "numeric", "factor", "factor"),
                       header = TRUE)

Ceuta_VR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/Ceuta_Survival_rates.txt",
                      colClasses = c("factor", "numeric", "factor", "factor"),
                      header = TRUE)

###### ASR estimation #######################################################
# Define WfP vital rates estimated from mark-recapture analysis:
Ceuta_ASR_output <- numeric(1000)
for(i in 1:1000){
VR <- list(#F_Chk_survl = Ceuta_VR[which(Ceuta_VRCeuta_VR$iter == i),2][5],
           F_Juv_survl = Ceuta_VR[which(Ceuta_VR$iter == i),2][3],
           F_Adt_survl = Ceuta_VR[which(Ceuta_VR$iter == i),2][1],
           #M_Chk_survl = Ceuta_VR[which(Ceuta_VR$iter == i),2][6],
           M_Juv_survl = Ceuta_VR[which(Ceuta_VR$iter == i),2][4],
           M_Adt_survl = Ceuta_VR[which(Ceuta_VR$iter == i),2][2],
           # Define h (harem size, h < 1 is polyandry) and k (clutch size)
           h = 1,
           k = 3,
           # Define primary sex ratio (assumed to be 0.5)
           PSR = 0.5)

# WfP matrix:
matrix <- plover_matrix(VR, chick_surv = FALSE)

# Determine the ASR at the stable stage distribution (assume h is 1, which 
# is monogamy)
ASR_h_1 <- freq_dep_SSD_ASR(A = matrix, h = 1, k = 3)

# Extract ASR
Ceuta_ASR_output[i] <- ASR_h_1$ASR
}