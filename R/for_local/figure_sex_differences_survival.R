# Plotting script of bootstrap results for comparative analysis
# Sex-differences in survival
# Luke J. Eberhart-Phillips
# July 6, 2016

# Load the required libraries
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(extrafont)

# Find fonts from computer that are candara or Candara
font_import(pattern="[M/m]enlo", prompt = FALSE) 
fonts()
fonttable()
loadfonts() # load these into R

# Import the bootstrapped survival analysis results
# each row is the sex- and stage-specific survival estimate for a given iteration
KiP_VR <- read.table("output/Bootstrap/KiP_Survival_rates.txt",
                     colClasses = c("factor", "numeric", "factor", "factor"),
                     header = TRUE)

WfP_VR <- read.table("output/Bootstrap/WfP_Survival_rates.txt",
                     colClasses = c("factor", "numeric", "factor", "factor"),
                     header = TRUE)

MP_VR <- read.table("output/Bootstrap/MP_Survival_rates.txt",
                    colClasses = c("factor", "numeric", "factor", "factor"),
                    header = TRUE) 

Tuzla_VR <- read.table("output/Bootstrap/Tuzla_Survival_rates.txt",
                       colClasses = c("factor", "numeric", "factor", "factor"),
                       header = TRUE)

Maio_VR <- read.table("output/Bootstrap/Maio_Survival_rates.txt",
                      colClasses = c("factor", "numeric", "factor", "factor"),
                      header = TRUE)

Ceuta_VR <- read.table("output/Bootstrap/Ceuta_Survival_rates.txt",
                       colClasses = c("factor", "numeric", "factor", "factor"),
                       header = TRUE)

# stack all estimates into one dataframe
Vital_rates <- rbind(KiP_VR,
                     WfP_VR,
                     MP_VR,
                     Tuzla_VR,
                     Maio_VR,
                     Ceuta_VR)

# summarize survival by popualtion, sex, and stage
Population_Stage_Survival_Summary <- 
  Rmisc::summarySE(Vital_rates, measurevar = "estimate", 
                   groupvars = c("species", "Sex_Age"))

# Remove N, se, and ci (arbitray since they are biased by bootstrap N)
Population_Stage_Survival_Summary <- 
  Population_Stage_Survival_Summary[,-c(3,6,7)]

# Tidy up column names
colnames(Population_Stage_Survival_Summary) <- 
  c("Species", "Sex_Age", "Survival", "Std_dev")

# Export survival summary statistics
write.csv(x = Population_Stage_Survival_Summary, 
          file = "output/Bootstrap/Population_Stage_Survival_Summary.txt", row.names = FALSE)


# Function to calculate the sex differences in stage-specific survival by population
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
  colnames(sex_diff_surv_output) <- c("Stage", "Difference")
  sex_diff_surv_output
}

# apply the function to each population
Ceuta_sex_diff <- sex_diff_surv(Ceuta_VR)
Ceuta_sex_diff$Population <- c("Snowy")

Tuzla_sex_diff <- sex_diff_surv(Tuzla_VR)
Tuzla_sex_diff$Population <- c("Kentish-Tuzla")

MP_sex_diff <- sex_diff_surv(MP_VR)
MP_sex_diff$Population <- c("Madagascar")

Maio_sex_diff <- sex_diff_surv(Maio_VR)
Maio_sex_diff$Population <- c("Kentish-Maio")

WfP_sex_diff <- sex_diff_surv(WfP_VR)
WfP_sex_diff$Population <- c("White-fronted")

KiP_sex_diff <- sex_diff_surv(KiP_VR)
KiP_sex_diff$Population <- c("Kittlitz's")

# stack the results into one dataframe and tidy column names
All_pops_sex_diff <- rbind(Ceuta_sex_diff,
                           Tuzla_sex_diff,
                           MP_sex_diff,
                           Maio_sex_diff,
                           WfP_sex_diff,
                           KiP_sex_diff)

colnames(All_pops_sex_diff) <- c("Stage", "Difference", "Population")

# define the factor levels of the population variable so that the populations
# are in an order that reflects the ASR (male biased to female biased)
All_pops_sex_diff$Population <- 
  factor(All_pops_sex_diff$Population ,
         levels = c("Snowy",
                    "Kentish-Tuzla",
                    "Madagascar",
                    "Kentish-Maio",
                    "White-fronted",
                    "Kittlitz's"))

# summarize the sex differences according to population
sex_diff_summary <- 
  All_pops_sex_diff %>%
  group_by(Population, Stage) %>%
  summarise(avg = mean(Difference),
            median = median(Difference),
            var = var(Difference))

# Custom color palette for the plotting of Juvenile and Adult stats
cbPalette <- c("#BDBDBD", "#737373")

# Drop all chick estiamtes, since these are not used in the comparative
# analysis due to the lack of precision in chick survival from most of 
# the populations
All_pops_sex_diff_A_F <- filter(All_pops_sex_diff, Stage != "Chick")

# redefine the factor levels to acknowledge that "chick" does not exist
All_pops_sex_diff_A_F$Stage <- 
  factor(All_pops_sex_diff_A_F$Stage, levels = c("Juvenile", "Adult"))

# draw the background plot (i.e. the two shades of color for male or
# female bias)
Background <- 
  ggplot(aes(y = Difference, x = Stage, fill = Stage), data = All_pops_sex_diff_A_F) + 
  theme_bw() +
  annotate("rect", xmin=0, xmax=4, ymin=-0.43, ymax=0, alpha=0.7,
           fill=brewer.pal(8, "Dark2")[c(1)]) +
  annotate("rect", xmin=0, xmax=4, ymin=0, ymax=0.43, alpha=0.7,
           fill=brewer.pal(8, "Dark2")[c(2)]) +
  annotate("text", x = c(2), y = c(-0.37),
           label = c("\u2640"), size = 7,
           family="Menlo", vjust = c(0.5), hjust = c(0.5)) +
  annotate("text", x = c(2), y = c(0.37),
           label = c("\u2642"), size = 7,
           family="Menlo", vjust = c(0.5), hjust = c(0.5)) +
  facet_grid(. ~ Population) +
  theme(text = element_text(family="Menlo", colour = "white"),
        legend.position = "none",
        axis.title.x = element_text(size=12, vjust=-0.1),
        axis.text.x  = element_text(size=10, angle = 45, hjust = 1), 
        axis.title.y = element_text(size=12, hjust=0.5, vjust = 3.5),
        axis.text.y  = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.5, colour = "white"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_blank(),
        plot.margin = unit(c(0.5,1.75,0.5,0.5), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_text(size=11)) +
  scale_x_continuous(limits=c(0,4),breaks=c(0,1), labels=c("Juvenile", "Adult"), expand = c(0, 0)) +
  scale_y_continuous(limits=c(-0.43,0.43), expand = c(0, 0)) +
  xlab("Life-stage") + 
  ylab("Sex-bias in survival")

# draw the plot with the same parameters as the background plot but this time add data
Bootstrap_sex_diff_VR_plot <- 
  ggplot(aes(y = Difference, x = Stage, fill = Stage), data = All_pops_sex_diff_A_F) + 
  theme_bw() +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_boxplot(width = 0.2) +
  facet_grid(. ~ Population) +
  theme(text = element_text(family="Menlo"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.x = element_text(size=12, vjust=-0.1),
        axis.text.x  = element_text(size=10, angle = 45, hjust = 1), 
        axis.title.y = element_text(size=12, hjust=0.5, vjust = 3.5),
        axis.text.y  = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        panel.border = element_blank(),
        plot.margin = unit(c(0.5,1.75,0.5,0.5), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_text(size=11)) +
  scale_fill_manual(values = cbPalette) +
  scale_y_continuous(limits=c(-0.43,0.43), expand = c(0, 0)) +
  xlab("Life-stage") + 
  ylab("Sex-bias in survival")
Bootstrap_sex_diff_VR_plot

# blank plot for presentations
Bootstrap_sex_diff_VR_plot_blank <- 
  ggplot(aes(y = Difference, x = Stage, fill = Stage), data = All_pops_sex_diff_A_F) + 
  theme_bw() +
  geom_blank() +
  facet_grid(. ~ Population) +
  theme(text = element_text(family="Menlo"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.x = element_text(size=12, vjust=-0.1),
        axis.text.x  = element_text(size=10, angle = 45, hjust = 1), 
        axis.title.y = element_text(size=12, hjust=0.5, vjust = 3.5),
        axis.text.y  = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        panel.border = element_blank(),
        plot.margin = unit(c(0.5,1.75,0.5,0.5), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_text(size=11)) +
  scale_fill_manual(values = cbPalette) +
  scale_y_continuous(limits=c(-0.43,0.43), expand = c(0, 0)) +
  xlab("Life-stage") + 
  ylab("Sex-bias in survival")
Bootstrap_sex_diff_VR_plot_blank

# prepare export of plot
jpeg(filename = "figs/Sex-differences_in_survival.jpg",
     quality = 100,
     width = 10,
     height = 4.5, 
     units = "in",
     res = 500) 

# overlay the background and data plots on top of eachother and export to disk
grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) ) 
print( Background + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( Bootstrap_sex_diff_VR_plot + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
dev.off()

# prepare export of blank plot
jpeg(filename = "figs/Sex-differences_in_survival_blank.jpg",
     quality = 100,
     width = 10,
     height = 4.5, 
     units = "in",
     res = 500) 

# overlay the background and data plots on top of eachother and export to disk
grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) ) 
print( Background + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( Bootstrap_sex_diff_VR_plot_blank + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
dev.off()