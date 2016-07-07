# Plotting script of bootstrap results for comparative analysis
# Population ASR histograms
# Luke J. Eberhart-Phillips
# July 7, 2016

# Load the required libraries
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(extrafont)

# Find fonts from computer that are candara or Candara
font_import(pattern="[H/h]eiti", prompt = FALSE) 
fonts()
fonttable()
loadfonts() # load these into R

# import bootstrap results
Ceuta_ASR <- read.table("output/Bootstrap/Ceuta_ASR.txt",
                        colClasses = c("numeric", "factor"),
                        header = TRUE)
Tuzla_ASR <- read.table("output/Bootstrap/Tuzla_ASR.txt",
                        colClasses = c("numeric", "factor"),
                        header = TRUE)
Maio_ASR <- read.table("output/Bootstrap/Maio_ASR.txt",
                       colClasses = c("numeric", "factor"),
                       header = TRUE)
WfP_ASR <- read.table("output/Bootstrap/WfP_ASR.txt",
                      colClasses = c("numeric", "factor"),
                      header = TRUE)
KiP_ASR <- read.table("output/Bootstrap/KiP_ASR.txt",
                      colClasses = c("numeric", "factor"),
                      header = TRUE)
MP_ASR <- read.table("output/Bootstrap/MP_ASR.txt",
                     colClasses = c("numeric", "factor"),
                     header = TRUE)

# Set the confidence interval to 95%
CI <- 0.95

# arrange the bootstrap results and choose the middle 95% to find the 95% CI
Ceuta_ASR_95CI <- c(arrange(Ceuta_ASR, ASR)[25,1], arrange(Ceuta_ASR, ASR)[975,1])
Maio_ASR_95CI <- c(arrange(Maio_ASR, ASR)[25,1], arrange(Maio_ASR, ASR)[975,1])
Tuzla_ASR_95CI <- c(arrange(Tuzla_ASR, ASR)[25,1], arrange(Tuzla_ASR, ASR)[975,1])
WfP_ASR_95CI <- c(arrange(WfP_ASR, ASR)[25,1], arrange(WfP_ASR, ASR)[975,1])
KiP_ASR_95CI <- c(arrange(KiP_ASR, ASR)[25,1], arrange(KiP_ASR, ASR)[975,1])
MP_ASR_95CI <- c(arrange(MP_ASR, ASR)[25,1], arrange(MP_ASR, ASR)[975,1])

# Calculate the mean ASR from the bootstrap results
Ceuta_ASR_mean <- mean(Ceuta_ASR$ASR)
Maio_ASR_mean <- mean(Maio_ASR$ASR)
Tuzla_ASR_mean <- mean(Tuzla_ASR$ASR)
WfP_ASR_mean <- mean(WfP_ASR$ASR)
KiP_ASR_mean <- mean(KiP_ASR$ASR)
MP_ASR_mean <- mean(MP_ASR$ASR)

# Calculate the median ASR from the bootstrap results
Ceuta_ASR_median <- median(Ceuta_ASR$ASR)
Maio_ASR_median <- median(Maio_ASR$ASR)
Tuzla_ASR_median <- median(Tuzla_ASR$ASR)
WfP_ASR_median <- median(WfP_ASR$ASR)
KiP_ASR_median <- median(KiP_ASR$ASR)
MP_ASR_median <- median(MP_ASR$ASR)

# Other ways to calculate the 95% confidence interval
KiP_ASR_95CI_quan <- stats::quantile(KiP_ASR$ASR, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
KiP_ASR_95CI_rank <- c(arrange(KiP_ASR, ASR)[25,1], arrange(KiP_ASR, ASR)[975,1])
KiP_ASR_95CI_equa <- c(mean(KiP_ASR$ASR)-((sd(KiP_ASR$ASR)/sqrt(length(KiP_ASR$ASR)))*1.96),
                       mean(KiP_ASR$ASR)+((sd(KiP_ASR$ASR)/sqrt(length(KiP_ASR$ASR)))*1.96)) 

# stack all summaries
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

# Add population names
population <- c("Snowy",
                "Kentish-Tuzla",
                "Madagascar",
                "Kentish-Maio",
                "White-fronted",
                "Kittlitz's")

# create dataframe of summary
All_pops_CI <- as.data.frame(cbind(population, All_pops_CI, All_pops_mean))

# Tidy up dataframe
rownames(All_pops_CI) <- NULL
colnames(All_pops_CI) <- c("population", "lcl", "ucl", "mean")
All_pops_CI$lcl <- as.numeric(as.character(All_pops_CI$lcl))
All_pops_CI$ucl <- as.numeric(as.character(All_pops_CI$ucl))
All_pops_CI$mean <- as.numeric(as.character(All_pops_CI$mean))

# Add population names to ASR bootstrap data
Ceuta_ASR$population <- "Snowy"
Tuzla_ASR$population <- "Kentish-Tuzla"
Maio_ASR$population <- "Kentish-Maio"
WfP_ASR$population <- "White-fronted"
KiP_ASR$population <- "Kittlitz's"
MP_ASR$population <- "Madagascar"

# stack bootstrap results from all populations
All_ASR <- rbind(Ceuta_ASR,
                 Tuzla_ASR,
                 Maio_ASR,
                 WfP_ASR,
                 KiP_ASR,
                 MP_ASR)

# define correct levels of population variable to sort ASR from male to female
All_ASR$population <- 
  factor(All_ASR$population ,
         levels = c("Snowy",
                    "Kentish-Tuzla",
                    "Madagascar",
                    "Kentish-Maio",
                    "White-fronted",
                    "Kittlitz's"))

# export colated ASR bootstrap results
#write.table(All_ASR, file = "Data_files/ASR_bootstrap_raw_data.txt", sep = "\t")

# Comparative plot of bootstrap results of ASR
Bootstrap_histogram_all <- 
  ggplot() +
  coord_flip() +
  annotate("rect", xmin=0.15, xmax=0.5, ymin=0, ymax=170, alpha=0.7,
           fill= brewer.pal(8, "Dark2")[c(1)]) +
  annotate("rect", xmin=0.5, xmax=0.85, ymin=0, ymax=170, alpha=0.7,
           fill= brewer.pal(8, "Dark2")[c(2)]) +
  annotate("text", x = c(0.2), y = c(85),
           label = c("\u2640"), size = 7,
           family="Menlo", vjust = c(0.5), hjust = c(0.5)) +
  annotate("text", x = c(0.8), y = c(85),
           label = c("\u2642"), size = 7,
           family="Menlo", vjust = c(0.5), hjust = c(0.5)) +
  geom_histogram(binwidth = 0.01, data = All_ASR, aes(x = ASR), fill = "grey30") +
  geom_errorbarh(data = All_pops_CI, aes(y = 160, x = lcl, xmin = lcl, xmax = ucl), color = "black", size = 0.8, linetype = "solid") +
  facet_grid(. ~ population) +
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
        axis.title.y = element_text(size=12, vjust=-0.1),
        axis.text.y  = element_text(size=10), 
        axis.title.x = element_text(size=12),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0.5,1.75,1.95,0.7), "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size=11),
        panel.margin = unit(0.75, "lines")) +
  ylab("Frequency") +
  xlab("Adult sex ratio (proportion \u2642)") +
  scale_x_continuous(limits = c(0.15, 0.85), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 170), expand = c(0, 0))
Bootstrap_histogram_all

# Save plot
ggsave(Bootstrap_histogram_all, 
       filename = "ASR_bootstrap_histogram_All_pops.jpg", 
       path = "figs/",
       width = 10,
       height = 4.5, units = "in",
       dpi = 300)

# blank plot for presentations
Bootstrap_histogram_all_blank <- 
  ggplot() +
  coord_flip() +
  annotate("rect", xmin=0.15, xmax=0.5, ymin=0, ymax=170, alpha=0.7,
           fill= brewer.pal(8, "Dark2")[c(1)]) +
  annotate("rect", xmin=0.5, xmax=0.85, ymin=0, ymax=170, alpha=0.7,
           fill= brewer.pal(8, "Dark2")[c(2)]) +
  annotate("text", x = c(0.2), y = c(85),
           label = c("\u2640"), size = 7,
           family="Menlo", vjust = c(0.5), hjust = c(0.5)) +
  annotate("text", x = c(0.8), y = c(85),
           label = c("\u2642"), size = 7,
           family="Menlo", vjust = c(0.5), hjust = c(0.5)) +
  geom_blank(data = All_ASR, aes(x = ASR)) +
  facet_grid(. ~ population) +
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
        axis.title.y = element_text(size=12, vjust=-0.1),
        axis.text.y  = element_text(size=10), 
        axis.title.x = element_text(size=12),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0.5,1.75,1.95,0.7), "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size=11),
        panel.margin = unit(0.75, "lines")) +
  ylab("Frequency") +
  xlab("Adult sex ratio (proportion \u2642)") +
  scale_x_continuous(limits = c(0.15, 0.85), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 170), expand = c(0, 0))
Bootstrap_histogram_all_blank

# save blank plot
ggsave(Bootstrap_histogram_all_blank, 
       filename = "ASR_bootstrap_histogram_All_pops_blank.jpg", 
       path = "figs/",
       width = 10,
       height = 4.5, units = "in",
       dpi = 300)