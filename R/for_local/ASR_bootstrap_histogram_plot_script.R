# Comparative ASR bootstrap histogram plot

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