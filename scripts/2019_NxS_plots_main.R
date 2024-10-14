##########################################################################
## Figure prep
#########################################################################
## Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)
library(lme4)
library(emmeans)
library(car)

## Central figure theme
pubtheme <- theme_bw(base_size = 16) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(size = 1.5, fill = NA),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 16),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor.y = element_blank(),
        legend.text.align = 0)

## Load data and marginal mean + SE summary sheet
data <- read.csv("../data/2019_NxS_datasheet.csv",
                 stringsAsFactors = FALSE,
                 na.strings = "NA") %>%
  mutate(anet.mass = a400 / marea,
         vcmax.mass = vcmax25 / marea,
         jmax.mass = jmax25 / marea,
         treatment = factor(treatment, levels = c("C", "NO3", "AS", "S")))

plot_data <- subset(data, nrcs.code == "ACSA3")

## Remove outliers based on statistical models
data$a400[data$a400 < 0.2] <- NA
data$anet.mass[data$a400 < 0.2] <- NA
data$pnue[data$pnue < 0] <- NA
data$vcmax.chi[c(85)] <- NA


##########################################################################
## Nmass - soil N
##########################################################################
nmass <- lmer(leaf.n ~ soil.n.norm + mineral.pH + 
                (1 | site), data = subset(data, nrcs.code == "ACSA3"))
Anova(nmass)
nmass.trend <- data.frame(emmeans(nmass, ~1, "soil.n.norm", 
                                  at = list(soil.n.norm = seq(0, 30, 0.1))))

nmass.plot <- ggplot(data = plot_data, 
                     aes(x = soil.n.norm, y = leaf.n/100)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = nmass.trend, 
              aes(y = emmean/100, ymin = lower.CL/100, ymax = upper.CL/100),
              alpha = 0.3) +
  geom_smooth(data = nmass.trend, aes(y = emmean/100), size = 2, se = FALSE,
              color = "black", method = 'lm', linetype = "dashed") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0.015, 0.035), 
                     breaks = seq(0.015, 0.035, 0.005)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil N (μg N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"["mass"]*" (gN g"["dry_mass"]*""^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
nmass.plot

##########################################################################
## Nmass - soil pH
##########################################################################

nmass.ph.plot <- ggplot(data = plot_data,
                        aes(x = mineral.pH, y = leaf.n/100)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(3.5, 5.5), 
                     breaks = seq(3.5, 5.5, 0.5)) +
  scale_y_continuous(limits = c(0.015, 0.035), 
                     breaks = seq(0.015, 0.035, 0.005)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil pH")),
       y = expression(bold("N"["mass"]*" (gN g"["dry_mass"]*""^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
nmass.ph.plot

##########################################################################
## Marea - soil N
##########################################################################
marea <-  lmer(marea ~ soil.n.norm + mineral.pH + 
                 (1 | site), data = subset(data, nrcs.code == "ACSA3"))
Anova(marea)

marea.trend <- data.frame(emmeans(marea, ~1, "soil.n.norm", 
                                  at = list(soil.n.norm = seq(0, 30, 0.1))))

marea.plot <- ggplot(data = plot_data, 
                     aes(x = soil.n.norm, y = marea)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = marea.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.3) +
  geom_smooth(data = marea.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black", method = 'lm', linetype = "dashed") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil N (μg N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("M"[area]*" (g"["dry_mass"]*" m"["leaf"]*""^"-2"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
marea.plot

##########################################################################
## Marea - soil pH
##########################################################################
marea.ph.plot <- ggplot(data = plot_data,
                        aes(x = mineral.pH, y = marea)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(3.5, 5.5),
                     breaks = seq(3.5, 5.5, 0.5)) +
  scale_y_continuous(limits = c(0, 120),
                     breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil pH")),
       y = expression(bold("M"[area]*" (g"["dry_mass"]*" m"["leaf"]*""^"-2"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
marea.ph.plot

##########################################################################
## Narea - soil N
##########################################################################
narea <-  lmer(narea ~ soil.n.norm + mineral.pH + 
                 (1 | site), data = subset(data, nrcs.code == "ACSA3"))
Anova(narea)
narea.trend <- data.frame(emmeans(narea, ~1, "soil.n.norm", 
                                  at = list(soil.n.norm = seq(0, 30, 0.1))))


narea.plot <- ggplot(data = plot_data, aes(x = soil.n.norm, y = narea)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = narea.trend, 
              aes(y = emmean,
                  ymin = lower.CL, ymax = upper.CL),
              alpha = 0.3) +
  geom_smooth(data = narea.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black", method = "lm") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0.5, 2.5), breaks = seq(0.5, 2.5, 0.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil N (μg N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
narea.plot


##########################################################################
## Narea - soil pH
##########################################################################
narea.ph.plot <- ggplot(data = plot_data,
                        aes(x = mineral.pH, y = narea)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(3.5, 5.5), breaks = seq(3.5, 5.5, 0.5)) +
  scale_y_continuous(limits = c(0.5, 2.5), breaks = seq(0.5, 2.5, 0.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil pH")),
       y = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
narea.ph.plot

##########################################################################
## Anet,area - soil N
##########################################################################
a400 <- lmer(sqrt(a400) ~ soil.n.norm + mineral.pH + 
               (1 | site), data = subset(data, nrcs.code == "ACSA3"))
Anova(a400)

a400.plot <- ggplot(data = plot_data, aes(x = soil.n.norm, y = a400)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil N (μg N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("A"["net"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
a400.plot

a400.ph.plot <- ggplot(data = plot_data,
                        aes(x = mineral.pH, y = a400)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(3.5, 5.5), breaks = seq(3.5, 5.5, 0.5)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil pH")),
       y = expression(bold("A"["net"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
a400.ph.plot

##########################################################################
## Anet,area - leaf N
##########################################################################
a400.narea <- ggplot(data = plot_data, aes(x = narea, y = a400)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0.5, 2.5), breaks = seq(0.5, 2.5, 0.5)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")), 
       y = expression(bold("A"["net"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
a400.narea

##########################################################################
## Vcmax - soil N
##########################################################################
vcmax <- lmer(vcmax25 ~ soil.n.norm + mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))
Anova(vcmax)

vcmax.trend <- data.frame(emmeans(vcmax, ~1, "soil.n.norm", 
                                  at = list(soil.n.norm = seq(0, 30, 0.1))))

vcmax.plot <- ggplot(data = plot_data, aes(x = soil.n.norm, y = vcmax25)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = vcmax.trend,
              aes(y = emmean,
                  ymin = lower.CL, ymax = upper.CL),
              alpha = 0.3) +
  geom_smooth(data = vcmax.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil N (μg N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("V"["cmax25"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.plot

##########################################################################
## Vcmax - soil pH
##########################################################################
Anova(vcmax)

vcmax.ph.plot <- ggplot(data = plot_data,
                       aes(x = mineral.pH, y = vcmax25)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(3.5, 5.5), breaks = seq(3.5, 5.5, 0.5)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil pH")),
       y = expression(bold("V"["cmax25"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.ph.plot

##########################################################################
## Vcmax - leaf N
##########################################################################
vcmax.leaf <- lmer(vcmax25 ~ narea + (1 | site), 
                   subset(data, nrcs.code == "ACSA3"))
Anova(vcmax.leaf)

vcmax.narea.trend <- data.frame(emmeans(vcmax.leaf, ~1, "narea",
                                        at = list(narea = seq(0.7, 2.25, 0.01))))

vcmax.narea <- ggplot(data = plot_data, aes(x = narea, y = vcmax25)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = vcmax.narea.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_smooth(data = vcmax.narea.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black", method = 'lm') +
  scale_x_continuous(limits = c(0.5, 2.5), breaks = seq(0.5, 2.5, 0.5)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")),
       y = expression(bold("V"["cmax25"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.narea

##########################################################################
## Jmax - soil N
##########################################################################
jmax <- lmer(jmax25 ~ soil.n.norm + mineral.pH + 
               (1 | site), data = subset(data, nrcs.code == "ACSA3"))
Anova(jmax)
jmax.trend <- data.frame(emmeans(jmax, ~1, "soil.n.norm", 
                                  at = list(soil.n.norm = seq(0, 30, 0.1))))

jmax.plot <- ggplot(data = plot_data, aes(x = soil.n.norm, y = jmax25)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = jmax.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.3) +
  geom_smooth(data = jmax.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black", method = 'lm') +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil N (μg N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("J"["max25"]*" (μmol m"^"-2" ~ "s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title.x = element_text(size = 16))
jmax.plot

##########################################################################
## Jmax - soil pH
##########################################################################
Anova(jmax)

jmax.ph.plot <- ggplot(data = plot_data,
                        aes(x = mineral.pH, y = jmax25)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(3.5, 5.5), breaks = seq(3.5, 5.5, 0.5)) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil pH")),
       y = expression(bold("J"["max25"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
jmax.ph.plot

##########################################################################
## Jmax - leaf N
##########################################################################
jmax.leaf <- lmer(jmax25 ~ narea + (1 | site), 
                  data = subset(data, nrcs.code == "ACSA3"))
Anova(jmax.leaf)

jmax.narea.trend <- data.frame(emmeans(jmax.leaf, ~1, "narea",
                                       at = list(narea = seq(0.7, 2.25, 0.01))))

jmax.narea <- ggplot(data = plot_data, aes(x = narea, y = jmax25)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = jmax.narea.trend,
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_smooth(data = jmax.narea.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black", method = 'lm') +
  scale_x_continuous(limits = c(0.5, 2.5), breaks = seq(0.5, 2.5, 0.5)) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")),
       y = expression(bold("J"["max25"]*" (μmol m"^"-2" ~ "s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title.x = element_text(size = 16))
jmax.narea

##########################################################################
## chi - soil N
##########################################################################
data$chi[66] <- NA

chi <- lmer(chi ~ soil.n.norm + mineral.pH + 
               (1 | site), data = subset(data, nrcs.code == "ACSA3"))
Anova(chi)

chi.trend <- data.frame(emmeans(chi, ~1, "soil.n.norm", 
                                 at = list(soil.n.norm = seq(0, 30, 0.1))))

chi <- ggplot(data = plot_data, aes(x = soil.n.norm, y = chi)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = chi.trend, 
              aes(y = emmean,
                  ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_smooth(data = chi.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0.57, 1.0), breaks = seq(0.6, 1.0, 0.1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil N (μg N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold(chi*" (unitless)")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
chi

##########################################################################
## chi - soil pH
##########################################################################
Anova(chi)

chi.ph.plot <- ggplot(data = plot_data,
                       aes(x = mineral.pH, y = chi)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(3.5, 5.5), breaks = seq(3.5, 5.5, 0.5)) +
  scale_y_continuous(limits = c(0.57, 1.0), breaks = seq(0.6, 1.0, 0.1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil pH")),
       y = expression(bold(chi*" (unitless)")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
chi.ph.plot

##########################################################################
## chi - leaf N
##########################################################################
chi.nleaf <- lmer(chi ~ narea + (1 | site), data = subset(data, nrcs.code == "ACSA3"))
Anova(chi.nleaf)

chi.narea.trend <- data.frame(emmeans(chi.nleaf, ~1, "narea",
                                      at = list(narea = seq(0.7, 2.25, 0.01))))

chi.narea <- ggplot(data = plot_data, aes(x = narea, y = chi)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = chi.narea.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_smooth(data = chi.narea.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black") +
  scale_x_continuous(limits = c(0.5, 2.5), breaks = seq(0.5, 2.5, 0.5)) +
  scale_y_continuous(limits = c(0.57, 0.9), breaks = seq(0.6, 0.9, 0.1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")),
       y = expression(bold(chi*" (unitless)")), shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
chi.narea

##########################################################################
## PNUE - soil N
##########################################################################
pnue <- lmer(sqrt(pnue) ~ soil.n.norm + mineral.pH + 
               (1 | site), data = subset(data, nrcs.code == "ACSA3"))
Anova(pnue)
pnue.trend <- data.frame(emmeans(pnue, ~1, "soil.n.norm",
                                 at = list(soil.n.norm = seq(0, 30, 0.1)),
                                 type = "response"))

pnue.plot <- ggplot(data = plot_data, aes(x = soil.n.norm, y = pnue)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = pnue.trend, aes(y = response, ymin = lower.CL, 
                  ymax = upper.CL), alpha = 0.3) +
  geom_smooth(data = pnue.trend, aes(y = response), size = 2, se = FALSE,
              color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil N (μg N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("PNUE (μmol CO"["2"]*" mol"^"-1"*"N s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title.y = element_text(size = 13))
pnue.plot

##########################################################################
## PNUE - soil pH
##########################################################################
Anova(pnue)

pnue.ph.plot <- ggplot(data = plot_data,
                       aes(x = mineral.pH, y = pnue)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(3.5, 5.5), breaks = seq(3.5, 5.5, 0.5)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil pH")),
       y = expression(bold("PNUE (μmol CO"["2"]*" mol"^"-1"*"N s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title.y = element_text(size = 13))
pnue.ph.plot

##########################################################################
## Narea-chi - soil N
##########################################################################
narea.chi <- lmer(narea.chi ~ soil.n.norm + mineral.pH +
                    (1 | site), data = subset(data, nrcs.code == "ACSA3"))
Anova(narea.chi)

narea.chi.trend <- data.frame(emmeans(narea.chi, ~1, "soil.n.norm",
                                      at = list(soil.n.norm = seq(0, 30,0.1)),
                                      type = "response"))

narea.chi.plot <- ggplot(data = plot_data, aes(x = soil.n.norm, y = narea.chi)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = narea.chi.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_smooth(data = narea.chi.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil N (μg N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"["area"]*" : "*chi*" (gN m"["leaf"]*""^"-2"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
narea.chi.plot

##########################################################################
## Narea.chi - soil pH
##########################################################################
Anova(narea.chi)

narea.chi.ph.plot <- ggplot(data = plot_data,
                       aes(x = mineral.pH, y = narea.chi)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(3.5, 5.5), breaks = seq(3.5, 5.5, 0.5)) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil pH")),
       y = expression(bold("N"["area"]*" : "*chi*" (gN m"["leaf"]*""^"-2"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
narea.chi.ph.plot

##########################################################################
## Vcmax-chi - soil N
##########################################################################
vcmax.chi <- lmer(vcmax.chi ~ soil.n.norm + mineral.pH + 
                     (1 | site), data = subset(data, nrcs.code == "ACSA3"))
Anova(vcmax.chi)

vcmax.chi.trend <- data.frame(emmeans(vcmax.chi, ~1, "soil.n.norm",
                                      at = list(soil.n.norm = seq(0, 30,0.1)),
                                      type = "response"))

vcmax.chi.plot <- ggplot(data = plot_data, aes(x = soil.n.norm, y = vcmax.chi)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = vcmax.chi.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_smooth(data = vcmax.chi.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 40)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil N (μg N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("V"["cmax25"]*" : "*chi*" (μmol m"^"-2"*"s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.chi.plot

##########################################################################
## Narea.chi - soil pH
##########################################################################
Anova(vcmax.chi)

vcmax.chi.ph.plot <- ggplot(data = plot_data,
                            aes(x = mineral.pH, y = vcmax.chi)) +
  geom_point(data = plot_data, aes(shape = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(3.5, 5.5), breaks = seq(3.5, 5.5, 0.5)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 40)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil pH")),
       y = expression(bold("V"["cmax25"]*" : "*chi*" (μmol m"^"-2"*"s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.chi.ph.plot

##########################################################################
## Vcmax-chi - leaf N
##########################################################################
vcmax.chi.nleaf <- lmer(vcmax.chi ~ narea + (1 | site), 
                      data = subset(data, nrcs.code == "ACSA3"))
Anova(vcmax.chi.nleaf)

vcmax.chi.nleaf.trend <- data.frame(emmeans(vcmax.chi.nleaf, ~1, "narea",
                                         at = list(narea = seq(0.7, 2.25, 0.01))))

vcmax.chi.narea <- ggplot(data = plot_data, aes(x = narea, y = vcmax.chi)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = vcmax.chi.nleaf.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_smooth(data = vcmax.chi.nleaf.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black") +
  scale_x_continuous(limits = c(0.5, 2.5), breaks = seq(0.5, 2.5, 0.5)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 40)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")),
       y = expression(bold("V"["cmax25"]*" : "*chi*" (μmol m"^"-2"*"s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.chi.narea

##########################################################################
## Figure 1: leaf N
##########################################################################
# png("[insert path here]",
#     width = 10, height = 12, units = 'in', res = 600)
ggarrange(narea.plot, narea.ph.plot,
          nmass.plot, nmass.ph.plot,
          marea.plot, marea.ph.plot,
          ncol = 2, nrow = 3, align = "hv",
          common.legend = TRUE, legend = "right", 
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
          font.label = list(size = 18, face = "bold"))
# dev.off()

##########################################################################
## Figure 2: Leaf biochemistry
##########################################################################
# png("[insert path here]",
#     width = 13, height = 11, units = 'in', res = 600)
ggarrange(a400.plot, a400.ph.plot, a400.narea,
          vcmax.plot, vcmax.ph.plot, vcmax.narea,
          jmax.plot, jmax.ph.plot, jmax.narea,
          ncol = 3, nrow = 3, align = "hv",
          common.legend = TRUE, legend = "right", 
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)",
                     "(f)", "(g)", "(h)", "(i)"),
          font.label = list(size = 18, face = "bold"))
# dev.off()

##########################################################################
## Figure 3: PNUE/iWUE 
##########################################################################
# png("[insert path here]",
#     width = 10, height = 15, units = 'in', res = 600)
ggarrange(chi, chi.ph.plot,
          pnue.plot, pnue.ph.plot, 
          narea.chi.plot, narea.chi.ph.plot,
          vcmax.chi.plot, vcmax.chi.ph.plot,
          ncol = 2, nrow = 4, align = "hv",
          common.legend = TRUE, legend = "right", 
          labels = c("(a)", "(b)", "(c)", "(d)", 
                     "(e)", "(f)", "(g)", "(h)"),
          font.label = list(size = 18, face = "bold"))
# dev.off()

##########################################################################
## Figure 3: Relationships between Narea and chi/vcmax:chi
##########################################################################
# png("[insert path here]",
#     width = 6, height = 8, units = 'in', res = 600)
ggarrange(chi.narea, vcmax.chi.narea,
          ncol = 1, nrow = 2, align = "hv",
          common.legend = TRUE, legend = "right", 
          labels = c("(a)", "(b)"),
          font.label = list(size = 18, face = "bold"))
# dev.off()
