library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)

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

Vital_rates <- rbind(KiP_VR,
                     WfP_VR,
                     MP_VR,
                     Tuzla_VR,
                     Maio_VR,
                     Ceuta_VR)

Population_Stage_Survival_Summary <- Rmisc::summarySE(Vital_rates, measurevar = "estimate", groupvars = c("species", "Sex_Age"))
Population_Stage_Survival_Summary <- Population_Stage_Survival_Summary[,-c(3,6,7)]
colnames(Population_Stage_Survival_Summary) <- c("Species", "Sex_Age", "Survival", "Std_dev")
write.csv(x = Population_Stage_Survival_Summary, file = "/home/luke/comparative_ASR/Bootstrap/Total_output/Population_Stage_Survival_Summary.txt", row.names = FALSE)
###### ASR estimation #######################################################
# Define WfP vital rates estimated from mark-recapture analysis:
sex_diff_surv <- function(VR) {
  sex_diff_surv_output <- data.frame(Adult = numeric(1000),
                                     Juvenile = numeric(1000),
                                     Chick = numeric(1000))
  for(i in 1:1000){
    Adult <- VR[which(VR$iter == i), 2][2] - VR[which(VR$iter == i), 2][1]
    Juvenile <- VR[which(VR$iter == i), 2][4] - VR[which(VR$iter == i), 2][3]
    Chick <- VR[which(VR$iter == i), 2][6] - VR[which(VR$iter == i), 2][5]
    
    sex_diff_surv_output[i, 1] <- Adult
    sex_diff_surv_output[i, 2] <- Juvenile
    sex_diff_surv_output[i, 3] <- Chick
  }
  sex_diff_surv_output <- melt(sex_diff_surv_output)
  colnames(sex_diff_surv_output) <- c("Stage", "Difference", "Population")
  sex_diff_surv_output
}

Ceuta_sex_diff <- sex_diff_surv(Ceuta_VR)
Ceuta_sex_diff$Population <- c("Snowy")

Tuzla_sex_diff <- sex_diff_surv(Tuzla_VR)
Tuzla_sex_diff$Population <- c("Kentish (Tuzla)")

MP_sex_diff <- sex_diff_surv(MP_VR)
MP_sex_diff$Population <- c("Madagascar")

Maio_sex_diff <- sex_diff_surv(Maio_VR)
Maio_sex_diff$Population <- c("Kentish (Maio)")

WfP_sex_diff <- sex_diff_surv(WfP_VR)
WfP_sex_diff$Population <- c("White-fronted")

KiP_sex_diff <- sex_diff_surv(KiP_VR)
KiP_sex_diff$Population <- c("Kittlitz's")

All_pops_sex_diff <- rbind(Ceuta_sex_diff,
                           Tuzla_sex_diff,
                           MP_sex_diff,
                           Maio_sex_diff,
                           WfP_sex_diff,
                           KiP_sex_diff)

colnames(All_pops_sex_diff) <- c("Stage", "Difference", "Population")

All_pops_sex_diff$Population <- 
  factor(All_pops_sex_diff$Population ,
         levels = c("Snowy",
                    "Kentish (Tuzla)",
                    "Madagascar",
                    "Kentish (Maio)",
                    "White-fronted",
                    "Kittlitz's"))

sex_diff_summary <- 
  All_pops_sex_diff %>%
  group_by(Population, Stage) %>%
  summarise(avg = mean(Difference),
            median = median(Difference),
            var = var(Difference))

cbPalette <- brewer.pal(3, "Dark2")
cbPalette <- c("#737373", "#BDBDBD")


All_pops_sex_diff_A_F <- filter(All_pops_sex_diff, Stage != "Chick")
All_pops_sex_diff_A_F$Stage <- 
  factor(All_pops_sex_diff_A_F$Stage, levels = c("Adult", "Juvenile"))

Background <- 
  ggplot(aes(y = Difference, x = Stage, fill = Stage), data = All_pops_sex_diff_A_F) + 
  coord_flip() +
  theme_bw() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, alpha=0.5,
           fill=brewer.pal(8, "Set1")[c(1)]) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=Inf, alpha=0.5,
           fill=brewer.pal(8, "Set1")[c(2)]) +
  annotate("text", x = c(2,2), y = c(-Inf, Inf),
           label = c("\u2640", "\u2642"), size = 7,
           family="Candara", vjust = c(0.5,0.5), hjust = c(-0.3,1.3)) +
  facet_grid(Population ~ .) +
  theme(text = element_text(family="Candara", color = "white"), # set the font as Candara
        legend.position = "none",
        axis.title.x = element_text(size=12, vjust=-0.1),
        axis.text.x  = element_text(size=10), 
        axis.title.y = element_text(size=12, hjust=0.5, vjust = 3.5),
        axis.text.y  = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.5, colour = "white"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_blank(),
        plot.margin = unit(c(2,1,1,1), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  scale_x_continuous(limits=c(0,4),breaks=c(0,1), labels=c("Fledgling", "Adult")) +
  scale_y_continuous(limits=c(-0.4,0.4)) +
  xlab("Life-stage") + 
  ylab("Sex-bias in survival")
#Background

Bootstrap_sex_diff_VR_plot <- 
  ggplot(aes(y = Difference, x = Stage, fill = Stage), data = All_pops_sex_diff_A_F) + 
  coord_flip() +
  theme_bw() +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  #geom_dotplot(binaxis = "y", binwidth = 0.003, stackdir = "center") +
  facet_grid(Population ~ .) +
  theme(text = element_text(family="Candara"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.x = element_text(size=12, vjust=-0.1),
        axis.text.x  = element_text(size=10), 
        axis.title.y = element_text(size=12, hjust=0.5, vjust = 3.5),
        axis.text.y  = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        panel.border = element_rect(linetype = "solid", colour = "grey"),
        plot.margin = unit(c(2,1,1,1), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  scale_fill_manual(values = cbPalette) +
  scale_y_continuous(limits=c(-0.4,0.4)) +
  xlab("Life-stage") + 
  ylab("Sex-bias in survival")
Bootstrap_sex_diff_VR_plot

jpeg(filename = "/home/luke/comparative_ASR/Bootstrap/Total_output/Figures/Sex-differences_in_survival.jpg",
     quality = 100,
     width = 5,
     height = 10, 
     units = "in",
     res = 300) 

grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) ) 
print( Background + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( Bootstrap_sex_diff_VR_plot + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
dev.off()


ggsave(#SP_fecund_plot, 
       filename = "Sex-differences_in_survival.jpg", 
       path = "/home/luke/comparative_ASR/Bootstrap/Total_output/Figures",
       width = 4,
       height = 4, units = "in",
       dpi = 300)