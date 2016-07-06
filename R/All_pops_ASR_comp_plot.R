library(RColorBrewer)
library(dplyr)

# import bootstrap results
Ceuta_ASR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/Ceuta_ASR.txt",
                        colClasses = c("numeric", "factor"),
                        header = TRUE)
Tuzla_ASR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/Tuzla_ASR.txt",
                        colClasses = c("numeric", "factor"),
                        header = TRUE)
Maio_ASR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/Maio_ASR.txt",
                       colClasses = c("numeric", "factor"),
                       header = TRUE)
WfP_ASR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/WfP_ASR.txt",
                      colClasses = c("numeric", "factor"),
                      header = TRUE)
KiP_ASR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/KiP_ASR.txt",
                      colClasses = c("numeric", "factor"),
                      header = TRUE)
MP_ASR <- read.table("/home/luke/comparative_ASR/Bootstrap/Total_output/MP_ASR.txt",
                     colClasses = c("numeric", "factor"),
                     header = TRUE)

# KiP_ASR$ASR <- KiP_ASR_output
# WfP_ASR$ASR <- WfP_ASR_output
# Tuzla_ASR$ASR <- Tuzla_ASR_output
# MP_ASR$ASR <- MP_ASR_output
# Maio_ASR$ASR <- Maio_ASR_output
Ceuta_ASR$ASR <- Ceuta_ASR_output
CI <- 0.95

Ceuta_ASR_95CI <- c(arrange(Ceuta_ASR, ASR)[25,1], arrange(Ceuta_ASR, ASR)[975,1])
Maio_ASR_95CI <- c(arrange(Maio_ASR, ASR)[25,1], arrange(Maio_ASR, ASR)[975,1])
Tuzla_ASR_95CI <- c(arrange(Tuzla_ASR, ASR)[25,1], arrange(Tuzla_ASR, ASR)[975,1])
WfP_ASR_95CI <- c(arrange(WfP_ASR, ASR)[25,1], arrange(WfP_ASR, ASR)[975,1])
KiP_ASR_95CI <- c(arrange(KiP_ASR, ASR)[25,1], arrange(KiP_ASR, ASR)[975,1])
MP_ASR_95CI <- c(arrange(MP_ASR, ASR)[25,1], arrange(MP_ASR, ASR)[975,1])

Ceuta_ASR_mean <- mean(Ceuta_ASR$ASR)
Maio_ASR_mean <- mean(Maio_ASR$ASR)
Tuzla_ASR_mean <- mean(Tuzla_ASR$ASR)
WfP_ASR_mean <- mean(WfP_ASR$ASR)
KiP_ASR_mean <- mean(KiP_ASR$ASR)
MP_ASR_mean <- mean(MP_ASR$ASR)

Ceuta_ASR_median <- median(Ceuta_ASR$ASR)
Maio_ASR_median <- median(Maio_ASR$ASR)
Tuzla_ASR_median <- median(Tuzla_ASR$ASR)
WfP_ASR_median <- median(WfP_ASR$ASR)
KiP_ASR_median <- median(KiP_ASR$ASR)
MP_ASR_median <- median(MP_ASR$ASR)

KiP_ASR_95CI_quan <- stats::quantile(KiP_ASR$ASR, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
KiP_ASR_95CI_rank <- c(arrange(KiP_ASR, ASR)[25,1], arrange(KiP_ASR, ASR)[975,1])
KiP_ASR_95CI_equa <- c(mean(KiP_ASR$ASR)-((sd(KiP_ASR$ASR)/sqrt(length(KiP_ASR$ASR)))*1.96),
                       mean(KiP_ASR$ASR)+((sd(KiP_ASR$ASR)/sqrt(length(KiP_ASR$ASR)))*1.96)) 

All_pops_CI <- rbind(Ceuta_ASR_95CI,
                     Tuzla_ASR_95CI,
                     MP_ASR_95CI,
                     Maio_ASR_95CI,
                     WfP_ASR_95CI,
                     KiP_ASR_95CI)

All_pops_mean <- rbind(Ceuta_ASR_mean,
                     Tuzla_ASR_mean,
                     MP_ASR_mean,
                     Maio_ASR_mean,
                     WfP_ASR_mean,
                     KiP_ASR_mean)

population <- c("Snowy",
                "Kentish (Tuzla)",
                "Madagascar",
                "Kentish (Maio)",
                "White-fronted",
                "Kittlitz's")

All_pops_CI <- as.data.frame(cbind(population, All_pops_CI, All_pops_mean))

rownames(All_pops_CI) <- NULL
colnames(All_pops_CI) <- c("population", "lcl", "ucl", "mean")
All_pops_CI$lcl <- as.numeric(as.character(All_pops_CI$lcl))
All_pops_CI$ucl <- as.numeric(as.character(All_pops_CI$ucl))
All_pops_CI$mean <- as.numeric(as.character(All_pops_CI$mean))

Ceuta_ASR$population <- "Snowy"
Tuzla_ASR$population <- "Kentish (Tuzla)"
Maio_ASR$population <- "Kentish (Maio)"
WfP_ASR$population <- "White-fronted"
KiP_ASR$population <- "Kittlitz's"
MP_ASR$population <- "Madagascar"

All_ASR <- rbind(Ceuta_ASR,
                 Tuzla_ASR,
                 Maio_ASR,
                 WfP_ASR,
                 KiP_ASR,
                 MP_ASR)

All_ASR$population <- 
  factor(All_ASR$population ,
         levels = c("Snowy",
                    "Kentish (Tuzla)",
                    "Madagascar",
                    "Kentish (Maio)",
                    "White-fronted",
                    "Kittlitz's"))

write.table(All_ASR, file = "Data_files/ASR_bootstrap_raw_data.txt", sep = "\t")

Bootstrap_histogram_all <- 
  ggplot() +
  annotate("rect", xmin=-Inf, xmax=0.5, ymin=-Inf, ymax=Inf, alpha=0.5,
           fill=brewer.pal(8, "Set1")[c(1)]) +
  annotate("rect", xmin=0.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.5,
           fill=brewer.pal(8, "Set1")[c(2)]) +
  annotate("text", x = c(-Inf,Inf), y = c(120, 120),
           label = c("\u2640", "\u2642"), size = 7,
           family="Candara", vjust = c(1.5,1.5), hjust = c(-0.5,1.5)) +
  geom_histogram(binwidth = 0.01, data = All_ASR, aes(x = ASR)) +
  geom_errorbarh(data = All_pops_CI, aes(y = 160, x = lcl, xmin = lcl, xmax = ucl), color = "black", size = 0.8, linetype = "solid") +
  facet_grid(population ~ .) +
  theme_bw() +
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
        axis.title.y = element_text(size=12),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        #strip.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "grey"),
        plot.margin = unit(c(2,1,1,1), "cm"),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.margin = unit(0.75, "lines")) +
  ylab("Frequency") +
  xlab("Adult sex ratio (proportion \u2642)") +
  scale_x_continuous(limits = c(0.2, 0.8)) +
  scale_y_continuous(limits = c(0, 170))
Bootstrap_histogram_all

ggsave(Bootstrap_histogram_all, 
       filename = "ASR_bootstrap_histogram_All_pops.jpg", 
       path = "/home/luke/comparative_ASR/Bootstrap/Total_output/Figures",
       width = 5,
       height = 10, units = "in",
       dpi = 300)

comp_ASR_plot <- 
  ggplot(comp_ASR, aes(x = pops, y = trans)) +  
  theme_bw() +
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