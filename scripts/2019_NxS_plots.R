## NOTE: Script assumes that the`scripts` folder is the root directory

## Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(lme4)
library(sjPlot)
library(patchwork)
library(gtable)
library(grid)

## Central figure theme
pubtheme <- theme_bw(base_size = 18) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(size = 1.5, fill = NA),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.title = element_text(face = "bold"),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor.y = element_blank(),
        legend.text.align = 0)

## Add colorblind friendly palette
cbbPalette <- c("#DDAA33", "#BB5566", "#004488", "#BBBBBB", "#FFFFFF")

## Load data and marginal mean + SE summary sheet
data <- read.csv("../data/2019_NxS_datasheet.csv",
                 stringsAsFactors = FALSE,
                 na.strings = "NA")
spp.data <- read.csv("../data/2019_NxS_figs_emmeanOutputs.csv")
soil.data <- read.csv("../data/2019_NxS_figs_soilemmeanOutputs.csv")
soil.data$.group <- tolower(soil.data$.group)

## Remove outliers based on statistical models
data$a400[data$a400 < 0.2] <- NA
data$vcmax25[11] <- NA
data$jmax25[11] <- NA
data$stom.lim[c(75)] <- NA
data$pnue[data$pnue < 0] <- NA
data$pnue[c(16, 102)] <- NA
data$iwue[data$iwue < 0] <- NA
data$iwue[c(16, 102)] <- NA
data$vcmax.gs[11] <- NA
data$ba.2011.2019[data$ba.2011.2019 <= 0] <- NA
data$growth.2011.2019[data$growth.2011.2019 <= 0] <- NA

## Subset by five primary species
data <- subset(data, nrcs.code == "ACRU" | nrcs.code == "ACSA3" |
                 nrcs.code == "QURU" | nrcs.code == "FAGR" | 
                 nrcs.code == "FRAM2")

# Add s.trt facet labels
facet.labels <- c("Sulfur added", "No sulfur added")
names(facet.labels) <- c("high.sulf", "low.sulf")

##########################################################################
## Categorical treatment combination effects on soil nitrogen availability
##########################################################################
soiln.plot.nsa <- ggplot(data = subset(soil.data, 
                                       variable == "soil.n" & s.trt == "low.sulf"), 
                         aes(x = factor(n.trt, levels = c("low.nit", "high.nit")),
                             y = emmean,
                             fill = factor(n.trt, levels = c("low.nit", "high.nit")))) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                position = position_dodge(0.9), width = 0.5, size = 1.5) +
  geom_point(aes(fill = n.trt), shape = 21, size = 6) +
  facet_grid(.~s.trt, labeller = labeller(s.trt = facet.labels)) +
  geom_text(aes(y = 30, label = .group), size = 7) +
  scale_y_continuous(limits = c(-0.5, 32), breaks = seq(0, 32, 8)) +
  scale_x_discrete(labels = c("low.nit" = "-N", "high.nit" = "+N")) +
  scale_fill_manual(values = c("#BB5566", "#004488")) +
  labs(x = NULL,
       y = expression(bold("Soil N (μg N g"^"-1"~"resin day"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
soiln.plot.nsa

soiln.plot.sa <- ggplot(data = subset(soil.data, 
                                      variable == "soil.n" & s.trt == "high.sulf"), 
                         aes(x = factor(n.trt, levels = c("low.nit", "high.nit")),
                             y = emmean, fill = factor(n.trt, levels = c("low.nit", 
                                                                         "high.nit")))) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                position = position_dodge(0.9), width = 0.5, size = 1.5) +
  geom_point(aes(fill = n.trt), shape = 21, size = 6) +
  facet_grid(.~s.trt, labeller = labeller(s.trt = facet.labels)) +
  geom_text(aes(y = 30, label = .group), size = 7) +
  scale_y_continuous(limits = c(-0.5, 32), breaks = seq(0, 32, 8)) +
  scale_x_discrete(labels = c("low.nit" = "-N", "high.nit" = "+N")) +
  scale_fill_manual(values = c("#BB5566", "#004488")) +
  labs(x = NULL, y = NULL) +
  guides(fill = "none") +
  pubtheme
soiln.plot.sa

##########################################################################
## Categorical treatment combination effects on soil pH
##########################################################################
ph.plot.nsa <- ggplot(data = subset(soil.data, 
                                   variable == "soil.pH" & s.trt == "low.sulf"),
                  aes(x = factor(n.trt, levels = c("low.nit", "high.nit")),
                      y = emmean, fill = factor(n.trt, levels = c("low.nit", "high.nit")))) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                position = position_dodge(0.9), width = 0.5, size = 1.5) +
  geom_point(shape = 21, size = 6) +
  facet_grid(.~s.trt, labeller = labeller(s.trt = facet.labels)) +
  geom_text(aes(y = 5.1, label = .group), size = 7) +
  scale_y_continuous(limits = c(3.9, 5.15), breaks = seq(3.9, 5.1, 0.3)) +
  scale_x_discrete(labels = c("low.nit" = "-N", "high.nit" = "+N")) +
  scale_fill_manual(values = c("#BB5566", "#004488")) +
  labs(x = NULL,
       y = expression(bold("Soil pH"))) +
  guides(fill = "none") +
  pubtheme +
  theme(strip.text = element_blank(),
        panel.grid.minor = element_blank())
ph.plot.nsa

ph.plot.sa <- ggplot(data = subset(soil.data, 
                                   variable == "soil.pH" & s.trt == "high.sulf"),
                     aes(x = factor(n.trt, levels = c("low.nit", "high.nit")),
                         y = emmean, 
                         fill = factor(n.trt, levels = c("low.nit", "high.nit")))) +
  #geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                position = position_dodge(0.9), width = 0.5, size = 1.5) +
  geom_point(shape = 21, size = 6) +
  facet_grid(.~s.trt, labeller = labeller(s.trt = facet.labels)) +
  geom_text(aes(y = 5.1, label = .group), size = 7) +
  scale_y_continuous(limits = c(3.9, 5.15), breaks = seq(3.9, 5.1, 0.3)) +
  scale_x_discrete(labels = c("low.nit" = "-N", "high.nit" = "+N")) +
  scale_fill_manual(values = c("#BB5566", "#004488")) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme +
  theme(strip.text = element_blank())
ph.plot.sa

##########################################################################
## Nleaf
##########################################################################
leaf.n <- lmer(leaf.n ~ soil.n.total.day + mineral.pH + nrcs.code + 
                 (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                             nrcs.code == "ACSA3" |
                                             nrcs.code == "QURU" |
                                             nrcs.code == "FAGR" | 
                                             nrcs.code == "FRAM2"))
leaf.n.fit <- data.frame(get_model_data(leaf.n, type = "pred", 
                                        terms = "soil.n.total.day"))


nleaf <- ggplot(data = data, aes(x = soil.n.total.day, y = leaf.n)) +
  geom_jitter(data = data, aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, 10)) +
  geom_ribbon(data = leaf.n.fit, aes(x = x, y = predicted,
                                     ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line(data = leaf.n.fit, aes(x = x, y = predicted), size = 1) +
  scale_y_continuous(limits = c(1, 4), 
                     breaks = seq(1, 4, 1)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("N"["mass"]~"(g g"^"-1"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
nleaf

##########################################################################
## SLA
##########################################################################
sla <- lmer(log(sla) ~ soil.n.total.day + mineral.pH * nrcs.code + 
              (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                          nrcs.code == "ACSA3" |
                                          nrcs.code == "QURU" |
                                          nrcs.code == "FAGR" | 
                                          nrcs.code == "FRAM2"))
sla.fit <- data.frame(get_model_data(sla, type = "pred", 
                                     terms = "soil.n.total.day"))

sla <- ggplot(data = data, aes(x = soil.n.total.day, y = sla)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  #geom_ribbon(data = sla.fit, aes(x = x, y = predicted, ymin = conf.low,
  #                                ymax = conf.high), alpha = 0.3) +
  #geom_line(data = slice(sla.fit, c(1,12)), 
  #          aes(x = x, y = predicted), size = 1) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, 150)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL, 
       y = expression(bold("SLA (cm"^"-2"~"g"^"-1"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
sla

##########################################################################
## Narea
##########################################################################
narea <- lmer(narea ~ soil.n.total.day + mineral.pH + nrcs.code +
                (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                            nrcs.code == "ACSA3" |
                                            nrcs.code == "QURU" |
                                            nrcs.code == "FAGR" | 
                                            nrcs.code == "FRAM2"))
narea.fit <- data.frame(get_model_data(narea, type = "pred", 
                                       terms = "soil.n.total.day"))

narea <- ggplot(data = data, aes(x = soil.n.total.day, y = narea)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  geom_ribbon(data = narea.fit, aes(x = x, y = predicted,
                                    ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  geom_line(data = narea.fit, aes(x = x, y = predicted), size = 1) +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("N"["area"]~"(g m"^"-2"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
narea

##########################################################################
## Net photosynthesis (area basis)
##########################################################################
a.area <- ggplot(data = data, aes(x = soil.n.total.day, y = a400)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 3)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("A"["net"]~ "(μmol m"^"-2"~"s"^"-1"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
a.area

##########################################################################
## Vcmax.area standardized to 25degC 
##########################################################################
vcmax <- ggplot(data = data, aes(x = soil.n.total.day, y = vcmax25)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("V"["cmax25"]~"(μmol m"^"-2" ~ "s"^"-1" ~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax

##########################################################################
## Jmax.area standardized to 25degC 
##########################################################################
jmax <- ggplot(data = data, aes(x = soil.n.total.day, y = jmax25)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 200), 
                     breaks = seq(0, 200, 50)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("J"["max25"]~"(μmol m"^"-2" ~ "s"^"-1" ~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
jmax

##########################################################################
## Jmax:Vcmax standardized to 25degC 
##########################################################################
vjmax <- lmer(log(jmax.vcmax) ~ soil.n.total.day + mineral.pH + nrcs.code + 
                (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                            nrcs.code == "ACSA3" |
                                            nrcs.code == "QURU" |
                                            nrcs.code == "FAGR" | 
                                            nrcs.code == "FRAM2"))

jmax.vcmax <- ggplot(data = data, aes(x = soil.n.total.day, y = jmax.vcmax)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(1, 3), breaks = seq(1, 3, 0.5)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("J"["max25"]~":V"["cmax25"])),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
jmax.vcmax

##########################################################################
## Stomatal conductance
##########################################################################
gs <- lmer(log(gsw) ~ soil.n.total.day + mineral.pH + nrcs.code + 
             (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                         nrcs.code == "ACSA3" |
                                         nrcs.code == "QURU" |
                                         nrcs.code == "FAGR" | 
                                         nrcs.code == "FRAM2"))

gs <- ggplot(data = data, aes(x = soil.n.total.day, y = gsw)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 0.18), breaks = seq(0, 0.18, 0.06)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("g"["s400"]~"(mol m"^"-2"~"s"^"-1"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
gs

##########################################################################
## chi
##########################################################################
chi <- lmer(chi ~ soil.n.total.day + mineral.pH + nrcs.code +
              (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                          nrcs.code == "ACSA3" |
                                          nrcs.code == "QURU" |
                                          nrcs.code == "FAGR" |
                                          nrcs.code == "FRAM2"))

chi <- ggplot(data = data, aes(x = soil.n.total.day, y = chi)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.2)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("χ (Pa Pa"^"-1"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title.y = element_text(size = 15))  
chi

##########################################################################
## Stomatal limitation
##########################################################################
l <- lmer(stom.lim ~ soil.n.total.day + mineral.pH + nrcs.code +
            (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                        nrcs.code == "ACSA3" |
                                        nrcs.code == "QURU" |
                                        nrcs.code == "FAGR" |
                                        nrcs.code == "FRAM2"))
l.fit <- data.frame(get_model_data(l, type = "pred", 
                                   terms = "soil.n.total.day"))

l <- ggplot(data = data, aes(x = soil.n.total.day, y = stom.lim)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  #geom_ribbon(data = l.fit, aes(x = x, y = predicted,
  #                              ymin = conf.low, ymax = conf.high), 
  #            alpha = 0.3) +
  #geom_line(data = l.fit, aes(x = x, y = predicted), 
  #          size = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0.7, 1), breaks = seq(0.7, 1, 0.1)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("Stomatal limitation")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
l

##########################################################################
## PNUE
##########################################################################
pnue <- lmer(sqrt(pnue) ~ soil.n.total.day + mineral.pH + nrcs.code +
               (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                           nrcs.code == "ACSA3" |
                                           nrcs.code == "QURU" |
                                           nrcs.code == "FAGR" |
                                           nrcs.code == "FRAM2"))
pnue.fit <- data.frame(get_model_data(pnue, type = "pred",
                                      terms = "soil.n.total.day"))
pnue <- ggplot(data = data, aes(x = soil.n.total.day, y = pnue)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  geom_ribbon(data = pnue.fit, aes(x = x, y = predicted,
                                   ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  geom_line(data = slice(pnue.fit, c(1,12)), 
            aes(x = x, y = predicted), size = 1) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, 3)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("PNUE (μmol CO"["2"]~" gN"^"-1"~"s"^"-1"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title.y = element_text(size = 15))  
pnue

##########################################################################
## iWUE
##########################################################################
iwue <- lmer(iwue ~ soil.n.total.day + mineral.pH + nrcs.code +
               (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                           nrcs.code == "ACSA3" |
                                           nrcs.code == "QURU" |
                                           nrcs.code == "FAGR" |
                                           nrcs.code == "FRAM2"))

iwue <- ggplot(data = data, aes(x = soil.n.total.day, y = iwue)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, 60)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("A"["400"]~": g"["s400"]~"(μmol CO"["2"]~" mol"^"-1"~"H"["2"]~"O)")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
iwue

##########################################################################
## narea-gs
##########################################################################
narea.gs <- lmer(log(narea.gs) ~ soil.n.total.day + mineral.pH + nrcs.code +
                   (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                               nrcs.code == "ACSA3" |
                                               nrcs.code == "QURU" |
                                               nrcs.code == "FAGR" |
                                               nrcs.code == "FRAM2"))
narea.gs.fit <- data.frame(get_model_data(narea.gs, type = "pred",
                                          terms = "soil.n.total.day"))

narea.gs <- ggplot(data = data, aes(x = soil.n.total.day, y = narea.gs)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  geom_ribbon(data = narea.gs.fit, 
              aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  geom_line(data = slice(narea.gs.fit, c(1,12)), 
            aes(x = x, y = predicted), size = 1) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("N"["area"]~":g"["s"]~"(gN s mol"^"-1"~"H"["2"]~"O)")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme + 
  theme(axis.title.y = element_text(size = 15))  
narea.gs

##########################################################################
## vcmax-gs
##########################################################################
vcmax.gs <- lmer(log(vcmax.gs) ~ soil.n.total.day * mineral.pH + nrcs.code +
                   (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                               nrcs.code == "ACSA3" |
                                               nrcs.code == "QURU" |
                                               nrcs.code == "FAGR" |
                                               nrcs.code == "FRAM2"))


vcmax.gs.fit <- data.frame(get_model_data(vcmax.gs, type = "pred",
                                          terms = "soil.n.total.day"))
vcmax.gs <- ggplot(data = data, aes(x = soil.n.total.day, y = vcmax.gs)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  geom_ribbon(data = vcmax.gs.fit, 
              aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  geom_line(data = slice(vcmax.gs.fit, c(1,12)), 
            aes(x = x, y = predicted), size = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 8000), breaks = seq(0, 8000, 2000)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("V"["cmax"]~":g"["s"]~"(μmol CO"["2"]~"mol"^"-1"~"H"["2"]~"O)")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title.y = element_text(size = 15))
vcmax.gs

##########################################################################
## Change in basal area
##########################################################################
ba <- lmer(sqrt(ba.2011.2019) ~ soil.n.total.day + mineral.pH + nrcs.code + 
             (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                         nrcs.code == "ACSA3" |
                                         nrcs.code == "QURU" |
                                         nrcs.code == "FAGR" |
                                         nrcs.code == "FRAM2"))

ba <- ggplot(data = data, aes(x = soil.n.total.day, y = ba.2011.2019)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 140), breaks = seq(0, 140, 35)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("BA"["11_19"]~"(cm"^"2"~"yr"^"-1"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
ba

##########################################################################
## Relative growth rate
##########################################################################
growth <- lmer(sqrt(growth.2011.2019) ~ soil.n.total.day + mineral.pH + nrcs.code +
                 (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                             nrcs.code == "ACSA3" |
                                             nrcs.code == "QURU" |
                                             nrcs.code == "FAGR" |
                                             nrcs.code == "FRAM2"))

rgr <- ggplot(data = data, aes(x = soil.n.total.day, y = growth.2011.2019)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 3, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("AS" = "+N; +S",
                                                            "NO3" = "+N; -S",
                                                            "S" = "-N; +S",
                                                            "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("RGR"["11_19"]~"(kg kg"^"-1"~"yr"^"-1"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
rgr

##########################################################################
## Figure 1: Categorical nutrient treatment effects on soil N and pH
##########################################################################
png("../working_drafts/figs/NxS_fig1_categorical.png",
    width = 8, height = 8, units = 'in', res = 600)
ggarrange(soiln.plot.nsa, soiln.plot.sa, ph.plot.nsa, ph.plot.sa,
          ncol = 2, nrow = 2, common.legend = TRUE, align = "v",
          legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold")) %>%
  annotate_figure(bottom = text_grob(
    expression(bold("Soil nitrogen addition")),
    size = 18))
dev.off()

##########################################################################
## Figure 2: leaf N allocation
##########################################################################
png("../working_drafts/figs/NxS_fig2_leafn.png",
    width = 12, height = 6, units = 'in', res = 600)
ggarrange(narea, ggarrange(nleaf, sla, ncol = 1, nrow = 2,
                           legend = "none", align = "hv", labels = c("B", "C"),
                           font.label = list(size = 18, face = "bold")),
          common.legend = TRUE, legend = "right", labels = c("A"),
          widths = c(1.5, 1),
          font.label = list(size = 18, face = "bold")) %>%
  annotate_figure(bottom = text_grob(
            expression(bold("Soil nitrogen (μg N g"^"-1"~"resin day"^"-1"~")")),
            size = 18))
dev.off()

##########################################################################
## Figure 2: Leaf biochemistry
##########################################################################
png("../working_drafts/figs/NxS_fig3_leafbiochem.png",
    width = 10, height = 8, units = 'in', res = 600)
ggarrange(a.area, vcmax, jmax, jmax.vcmax,
          ncol = 2, nrow = 2, common.legend = TRUE,
          align = "v", legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold")) %>%
  annotate_figure(bottom = text_grob(
    expression(bold("Soil nitrogen (μg N g"^"-1"~"resin day"^"-1"~")")),
    size = 18))
dev.off()

##########################################################################
## Figure 3: PNUE/iWUE tradeoffs
##########################################################################
png("../working_drafts/figs/NxS_fig4_pnueiwue.png",
    width = 10, height = 8, units = 'in', res = 600)
ggarrange(pnue, chi, narea.gs, vcmax.gs, ncol = 2, nrow = 2, 
          common.legend = TRUE, align = "v", legend = "right", 
          labels = "AUTO", font.label = list(size = 18, face = "bold")) %>%
  annotate_figure(bottom = text_grob(
    expression(bold("Soil nitrogen (μg N g"^"-1"~"resin day"^"-1"~")")),
    size = 18))
dev.off()

##########################################################################
## Figure 4: Whole plant measurements
##########################################################################
png("../working_drafts/figs/NxS_fig5_wp.png",
    width = 10, height = 4.5, units = 'in', res = 600)
ggarrange(ba, rgr, ncol = 2, nrow = 1, common.legend = TRUE,
          align = "v", legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold")) %>%
  annotate_figure(bottom = text_grob(
                    expression(
                      bold("Soil nitrogen (μg N g"^"-1"~"resin day"^"-1"~")")),
                    size = 18))
dev.off()


##########################################################################
## Presentation Fig. 1: Leaf N and net photosynthesis
##########################################################################
ggarrange(narea, a.area, ncol = 2, nrow = 1, common.legend = TRUE,
          align = "hv", legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold")) %>%
  annotate_figure(bottom = text_grob(
    expression(
      bold("Soil nitrogen (μg N g"^"-1"~"resin day"^"-1"~")")),
    size = 15)) %>%
  ggexport(filename = "../working_drafts/figs/NxS_talkfig1_leafNphoto.jpg", 
           width = 6000, height = 2250, res = 600)


##########################################################################
## Presentation Fig. 2: PNUE, iWUE
##########################################################################
ggarrange(pnue, chi, ncol = 2, nrow = 1, common.legend = TRUE,
          align = "hv", legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold")) %>%
  annotate_figure(bottom = text_grob(
    expression(
      bold("Soil nitrogen (μg N g"^"-1"~"resin day"^"-1"~")")),
    size = 15)) %>%
  ggexport(filename = "../working_drafts/figs/NxS_talkfig2_pnueiwue.jpg", 
           width = 6000, height = 2250, res = 600)

##########################################################################
## Presentation Fig. 3: Narea:gs
##########################################################################
ggarrange(narea.gs, vcmax.gs, ncol = 2, nrow = 1, common.legend = TRUE,
          align = "hv", legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold")) %>%
  annotate_figure(bottom = text_grob(
    expression(
      bold("Soil nitrogen (μg N g"^"-1"~"resin day"^"-1"~")")),
    size = 15)) %>%
  ggexport(filename = "../working_drafts/figs/NxS_talkfig3_nareags.jpg", 
           width = 6500, height = 2550, res = 600)



