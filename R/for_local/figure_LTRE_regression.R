# Comparative elasticity regression plot

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
lm_eqn$VR <- c("Fledgling", "Adult")
row.names(lm_eqn) <- NULL
colnames(lm_eqn) <- c("V1", "Sex", "VR")

Comp_Elasticity_plot <- 
  ggplot(filter(plover_sens, variable == "Elasticity" & VR != "Chick" & 
                  Vital_rate != "PSR" & Vital_rate != "h"),
         aes(x = value_trans, y = data)) +
  geom_point(size = I(4), alpha = I(0.7)) + 
  geom_smooth(method = lm, aes(colour = Sex, fill = Sex), size=2) +
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
  geom_text(data = lm_eqn, aes(x = -Inf, y = 1,label = V1, family="Candara"), 
            parse = TRUE, inherit.aes=FALSE, hjust = -0.05)

# Export 
ggsave(Comp_Elasticity_plot, 
       filename = "Comparative_Elasticity_ASR_plot.jpg", 
       path = "L:/PhD/Plover_Matrix_Modelling/Plover_Matrix_Modelling/Figures",
       width = 12,
       height = 6, units = "in",
       dpi = 300)