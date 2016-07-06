# Plotting script of bootstrap results for comparative analysis
# Luke J. Eberhart-Phillips
# July 6, 2016

# Load the required libraries
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(extrafont)

# Violin plot of sex-differences in survival for 1st years and adults
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