##########################################################################
## Figure prep
#########################################################################
## Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(lme4)
library(emmeans)

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
        legend.text = element_text(size = 18),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor.y = element_blank(),
        legend.text.align = 0)

## Add colorblind friendly palette
cbbPalette <- c("#FFFFFF", "#DDAA33", "#BB5566", "#004488", "#000000")

## Load data and marginal mean + SE summary sheet
data <- read.csv("../data/2019_NxS_datasheet.csv",
                 stringsAsFactors = FALSE,
                 na.strings = "NA")
  
spp.data <- read.csv("../data/2019_NxS_figs_emmeanOutputs.csv") %>%
  mutate(nrcs.code = ifelse(nrcs.code == "ACSA3",
                            "ACSA",
                            ifelse(nrcs.code == "FRAM2", 
                                   "FRAM",
                                   nrcs.code)),
         nrcs.code = factor(nrcs.code, 
                            levels = c("ACRU", "ACSA", "FRAM", 
                                       "FAGR", "QURU")))

## Remove outliers based on statistical models
data$narea[12] <- NA
data$a400[data$a400 < 0.2] <- NA
data$vcmax25[11] <- NA
data$jmax25[11] <- NA
data$chi[66] <- NA
data$stom.lim[c(68, 75)] <- NA
data$pnue[data$pnue < 0] <- NA
data$iwue[data$iwue < 0] <- NA
data$iwue[c(16, 102)] <- NA
data$p.bioe[11] <- NA
data$p.photo[11] <- NA
data$p.structure[49] <- NA
data$narea.chi[12] <- NA
data$vcmax.chi[c(11, 37)] <- NA

## Subset by five primary species
data <- subset(data, nrcs.code == "ACRU" | nrcs.code == "ACSA3" |
                 nrcs.code == "QURU" | nrcs.code == "FAGR" | 
                 nrcs.code == "FRAM2") %>%
  mutate(nrcs.code = ifelse(nrcs.code == "ACSA3",
                            "ACSA",
                            ifelse(nrcs.code == "FRAM2", 
                                   "FRAM",
                                   nrcs.code)),
         nrcs.code = factor(nrcs.code, 
                            levels = c("ACRU", "ACSA", "FRAM", 
                                       "FAGR", "QURU")))

## Create blank plot as spacer plot
blank.plot <- ggplot() + 
  theme_bw() +
  theme(panel.background = element_rect(color = "white",
                                        fill = "white"),
        panel.border = element_rect(color = "white"))

##########################################################################
## Nmass - soil N
##########################################################################
nmass <- lmer(leaf.n ~ soil.n.norm + mineral.pH + nrcs.code + 
                 (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                             nrcs.code == "ACSA" |
                                             nrcs.code == "QURU" |
                                             nrcs.code == "FAGR" | 
                                             nrcs.code == "FRAM"))
car::Anova(nmass)
nmass.trend <- data.frame(emmeans(nmass, ~1, "soil.n.norm", 
                                  at = list(soil.n.norm = c(0, 4, seq(20,30,0.1)))))

nleaf <- ggplot(data = data, aes(x = soil.n.norm, y = leaf.n/100)) +
  geom_jitter(data = data, aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  geom_ribbon(data = nmass.trend, 
              aes(y = emmean/100,
                  ymin = lower.CL/100, ymax = upper.CL/100),
              alpha = 0.3) +
  geom_smooth(data = nmass.trend, aes(y = emmean/100), size = 1, se = FALSE,
              color = "black", method = 'lm') +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0.015, 0.035), 
                     breaks = seq(0.015, 0.035, 0.005)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("N"["mass"]*" (gN g"["dry biomass"]*""^"-1"*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
nleaf

##########################################################################
## Nmass - SPECIES
##########################################################################
nleaf.spp <- ggplot(data = subset(spp.data, variable == "leaf.n"),
                    aes(x = nrcs.code, y = emmean/100, 
                        fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = leaf.n/100,
                               fill = nrcs.code, shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL/100, ymax =upper.CL/100),
                size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 0.035, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0.015, 0.035), breaks = seq(0.015, 0.035, 0.005)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = NULL,
       y = NULL,
       shape = "Treatment",
       fill = "Species") +
  pubtheme
nleaf.spp

##########################################################################
## Marea - soil N
##########################################################################
marea <- lmer(marea ~ soil.n.norm + mineral.pH + nrcs.code + 
                (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                            nrcs.code == "ACSA" |
                                            nrcs.code == "QURU" |
                                            nrcs.code == "FAGR" | 
                                            nrcs.code == "FRAM"))

marea.trend <- data.frame(emmeans(marea, ~1, "soil.n.norm", 
                                  at = list(soil.n.norm = c(0, 4, seq(20,30,0.1)))))



marea <- ggplot(data = data, aes(x = soil.n.norm, y = marea)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  geom_ribbon(data = marea.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_line(data = marea.trend, 
            aes(y = emmean), size = 1) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(20, 110), breaks = seq(20, 110, 30)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = expression(bold("Soil N (μg N g"^"-1"*" resin d"^"-1"*")")), 
       y = expression(bold("M"[area]*" (g"["dry biomass"]*" m"["leaf"]*""^"-2"*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
marea

##########################################################################
## Marea - SPECIES
##########################################################################
marea.spp <- ggplot(data = subset(spp.data, variable == "marea"),
                  aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = marea, fill = nrcs.code,
                               shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 110, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(20, 110), breaks = seq(20, 110, 30)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), 
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = "Species",
       y = NULL,
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
marea.spp

##########################################################################
## Narea - soil N
##########################################################################
narea <- lmer(narea ~ soil.n.norm + mineral.pH + nrcs.code +
                (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                            nrcs.code == "ACSA" |
                                            nrcs.code == "QURU" |
                                            nrcs.code == "FAGR" | 
                                            nrcs.code == "FRAM"))
narea.trend <- data.frame(emmeans(narea, ~1, "soil.n.norm", 
                                  at = list(soil.n.norm = c(0, 4, seq(20,30,0.1)))))


narea <- ggplot(data = data, aes(x = soil.n.norm, y = narea)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  geom_ribbon(data = narea.trend,
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_line(data = narea.trend, aes(y = emmean), size = 1) +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
narea

##########################################################################
## Narea - SPECIES
##########################################################################
narea.spp <- ggplot(data = subset(spp.data, variable == "narea"),
                    aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = narea, fill = nrcs.code,
                               shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 3, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = NULL,
       y = NULL,
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
narea.spp

##########################################################################
## A400 - soil N
##########################################################################
a400 <- ggplot(data = data, aes(x = soil.n.norm, y = a400)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, 4)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("A"["net"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
a400

##########################################################################
## A400 - SPECIES
##########################################################################
a400.spp <- ggplot(data = subset(spp.data, variable == "a400"),
                     aes(x = nrcs.code, y = emmean^2, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = a400, fill = nrcs.code,
                               shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL^2, ymax = upper.CL^2), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 16, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, 4)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), 
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = NULL,
       y = NULL,
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
a400.spp

##########################################################################
## A400 - leaf N
##########################################################################
a400.leaf <- lmer(sqrt(a400) ~ narea + (1 | nrcs.code) + (1 | site), data = data)

a400.nleaf.trend <- data.frame(emmeans(a400.leaf, ~1, "narea", 
                                       at = list(narea = seq(0.72, 2.67, 0.01)),
                                       type = "response"))

a400.narea <- ggplot(data = data, aes(x = narea, y = a400)) +
  geom_point(aes(fill = nrcs.code, shape = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = a400.nleaf.trend, 
              aes(y = response, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_line(data = a400.nleaf.trend, aes(y = response), size = 1) +
  scale_x_continuous(limits = c(0.5, 3), breaks = seq(0.5, 3, 0.5)) +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, 4)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = NULL,
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
a400.narea

##########################################################################
## Vcmax - soil N
##########################################################################
vcmax <- ggplot(data = data, aes(x = soil.n.norm, y = vcmax25)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("V"["cmax25"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax

##########################################################################
## Vcmax - SPECIES
##########################################################################
vcmax.spp <- ggplot(data = subset(spp.data, variable == "vcmax25"),
                    aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, 
              aes(x = nrcs.code, y = vcmax25, fill = nrcs.code, shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 120, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), 
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = NULL,
       y = NULL,
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.spp

##########################################################################
## Vcmax - leaf N
##########################################################################
vcmax.leaf <- lmer(vcmax25 ~ narea + (1 | nrcs.code) + (1 | site), data = data)

vcmax.narea.trend <- data.frame(emmeans(vcmax.leaf, ~1, "narea",
                                        at = list(narea = seq(0.7, 2.67, 0.01))))

vcmax.narea <- ggplot(data = data, aes(x = narea, y = vcmax25)) +
  geom_point(aes(fill = nrcs.code, shape = treatment), 
              size = 4, alpha = 0.75) +
  geom_ribbon(data = vcmax.narea.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_line(data = vcmax.narea.trend, aes(y = emmean), size = 1) +
  scale_x_continuous(limits = c(0.5, 3), breaks = seq(0.5, 3, 0.5)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = NULL,
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.narea

##########################################################################
## Jmax - soil N
##########################################################################
jmax <- ggplot(data = data, aes(x = soil.n.norm, y = jmax25)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = expression(bold("Soil N (μg N g"^"-1"*" resin d"^"-1"*")")),
       y = expression(bold("J"["max25"]*" (μmol m"^"-2" ~ "s"^"-1"*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title.x = element_text(size = 16))
jmax

##########################################################################
## Jmax - SPECIES
##########################################################################
jmax.spp <- ggplot(data = subset(spp.data, variable == "jmax25"),
                   aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = jmax25, fill = nrcs.code,
                               shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 200, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), 
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = "Species",
       y = NULL,
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + 
  pubtheme +
  theme(axis.title.x = element_text(size = 16))
jmax.spp

##########################################################################
## Jmax - leaf N
##########################################################################
jmax.leaf <- lmer(jmax25 ~ narea +  (1 | nrcs.code) + (1 | site), data = data)


jmax.narea.trend <- data.frame(emmeans(jmax.leaf, ~1, "narea",
                                       at = list(narea = c(0.7, 2.67, 0.01))))


jmax.narea <- ggplot(data = data, aes(x = narea, y = jmax25)) +
  geom_point(aes(fill = nrcs.code, shape = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = jmax.narea.trend, 
              aes(y = emmean,
                  ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_line(data = jmax.narea.trend, aes(y = emmean), 
            size = 1) +
  scale_x_continuous(limits = c(0.5, 3), breaks = seq(0.5, 3, 0.5)) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")),
       y = NULL,
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title.x = element_text(size = 16))
jmax.narea

##########################################################################
## chi - soil N
##########################################################################
chi <- ggplot(data = data, aes(x = soil.n.norm, y = chi)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold(chi*" (unitless)")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
chi

##########################################################################
## chi - SPECIES
##########################################################################
chi.spp <- ggplot(data = subset(spp.data, variable == "chi"),
                  aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data,
              aes(x = nrcs.code, y = chi, fill = nrcs.code, shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 1, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = NULL,
       y = NULL,
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme 
chi.spp

##########################################################################
## chi - leaf N
##########################################################################
chi.leaf <- lmer(chi ~ narea + (1 | nrcs.code) + (1 | site), data = data)

chi.narea.trend <- data.frame(emmeans(chi.leaf, ~1, "narea",
                                      at = list(narea = seq(0.7, 2.67, 0.01))))

chi.narea <- ggplot(data = data, aes(x = narea, y = chi)) +
  geom_point(aes(fill = nrcs.code, shape = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = chi.narea.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_line(data = chi.narea.trend, aes(y = emmean), size = 1) +
  scale_x_continuous(limits = c(0.5, 3), breaks = seq(0.5, 3, 0.5)) +
  scale_y_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = NULL,
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
chi.narea

##########################################################################
## PNUE - soil N
##########################################################################
pnue <- lmer(sqrt(pnue) ~ soil.n.norm + mineral.pH + nrcs.code +
               (1 | site), data = data)

pnue.trend <- data.frame(emmeans(pnue, ~1, "soil.n.norm",
                                 at = list(soil.n.norm = c(0, 4, seq(20,30,0.1))),
                                 type = "response"))

pnue.plot <- ggplot(data = data, aes(x = soil.n.norm, y = pnue)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 2, size = 4, alpha = 0.75) +
  geom_ribbon(data = pnue.trend,
              aes(y = response,
                  ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_line(data = pnue.trend, 
            aes(y = response), size = 1) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("PNUE (μmol CO"["2"]*" mol"^"-1"*"N s"^"-1"*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title = element_text(size = 15))  
pnue.plot

##########################################################################
## PNUE - SPECIES
##########################################################################
pnue.spp <- ggplot(data = subset(spp.data, variable == "pnue"),
                   aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = pnue, fill = nrcs.code,
                               shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 120, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = NULL,
       y = NULL,
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title = element_text(size = 15))
pnue.spp

##########################################################################
## Narea-chi - soil N
##########################################################################
narea.chi <- lmer(narea.chi ~ soil.n.norm + mineral.pH + nrcs.code +
                   (1 | site), data = data)

narea.chi.trend <- data.frame(emmeans(narea.chi, ~1, "soil.n.norm",
                                      at = list(soil.n.norm = c(0, 4, seq(20,30,0.1)))))

narea.chi.plot <- ggplot(data = data, aes(x = soil.n.norm, y = narea.chi)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  geom_ribbon(data = narea.chi.trend , 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_line(data = narea.chi.trend , 
            aes(y = emmean), size = 1) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 4.25), breaks = seq(0, 4, 1)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold("N"["area"]*" : "*chi*" (gN m"["leaf"]*""^"-2"*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
narea.chi.plot

##########################################################################
## Narea-chi - SPECIES
##########################################################################
narea.chi.spp <- ggplot(data = subset(spp.data, variable == "narea.chi"),
                       aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data,
              aes(x = nrcs.code, y = narea.chi, fill = nrcs.code, shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 4.25, label = .group), 
            fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 4.25), breaks = seq(0, 4, 1)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = NULL,
       y = NULL,
       fill = "Species",
       shape = "Treatment") +
  pubtheme
narea.chi.spp

##########################################################################
## Vcmax-chi - soil N
##########################################################################
vcmax.chi <- ggplot(data = data, aes(x = soil.n.norm, y = vcmax.chi)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, 60)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = expression(bold("Soil N (μg N g"^"-1"*" resin d"^"-1"*")")),
       y = expression(bold("V"["cmax"]*" : "*chi*" (μmol m"^"-2"*"s"^"-1"*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme #+
 # theme(axis.title = element_text(size = 15))
vcmax.chi

##########################################################################
## Vcmax-chi - SPECIES
##########################################################################
vcmax.chi.spp <- ggplot(data = subset(spp.data, variable == "vcmax.chi"),
                       aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = vcmax.chi, fill = nrcs.code,
                               shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 180, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, 60)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = "Species",
       y = NULL,
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.chi.spp

##########################################################################
## Vcmax-gs - leaf N
##########################################################################
vcmax.chi.leaf <- lmer(vcmax.chi ~ narea + (1 | nrcs.code) + (1 | site), 
                      data = data)
vcmax.chi.leaf.trend <- data.frame(emmeans(vcmax.chi.leaf, ~1, "narea",
                                         at = list(narea = seq(0.7, 2.67, 0.01))))


vcmax.chi.narea <- ggplot(data = data, aes(x = narea, y = vcmax.chi)) +
  geom_point(aes(fill = nrcs.code, shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = vcmax.chi.leaf.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3) +
  geom_line(data = vcmax.chi.leaf.trend, aes(y = emmean), size = 1) +
  scale_x_continuous(limits = c(0.5, 3), breaks = seq(0.5, 3, 0.5)) +
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, 60)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")),
       y = NULL,
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.chi.narea

##########################################################################
## N in photosynthesis - soil N
##########################################################################
p.photo <- ggplot(data = data, aes(x = soil.n.norm, y = p.photo)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold(rho["photo"]*" (gN"["photo"]* " gN"["leaf"]*""^-1*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
p.photo

##########################################################################
## N in photosynthesis - SPECIES
##########################################################################
p.photo.spp <- ggplot(data = subset(spp.data, variable == "p.photo"),
                       aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, 
              aes(x = nrcs.code, y = p.photo, fill = nrcs.code,
                  shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 0.8, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = NULL,
       y = NULL,
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
p.photo.spp

##########################################################################
## N in rubisco - soil N
##########################################################################
p.rubisco <- ggplot(data = data, aes(x = soil.n.norm, y = p.rubisco)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 0.65), breaks = seq(0, 0.6, 0.2)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold(rho["rub"]*" (gN"["rub"]* " gN"["leaf"]*""^-1*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
p.rubisco

##########################################################################
## Prop N in rubisco - SPECIES
##########################################################################
p.rubisco.spp <- ggplot(data = subset(spp.data, variable == "p.rubisco"),
                          aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, 
              aes(x = nrcs.code, y = p.rubisco, fill = nrcs.code,
                  shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 0.65, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 0.65), breaks = seq(0, 0.6, 0.2)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = NULL,
       y = NULL,
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
p.rubisco.spp

##########################################################################
## Prop N in bioenergetics - soil N
##########################################################################
p.bioenergetics <- ggplot(data = data, aes(x = soil.n.norm, y = p.bioe)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.025)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = NULL,
       y = expression(bold(rho["bioe"]*" (gN"["bioe"]* " gN"["leaf"]*""^-1*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
p.bioenergetics

##########################################################################
## Prop N in bioenergetics - SPECIES
##########################################################################
p.bio.spp <- ggplot(data = subset(spp.data, variable == "p.bioe"),
                        aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, 
              aes(x = nrcs.code, y = p.bioe, fill = nrcs.code,
                  shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 0.1, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.025)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = NULL,
       y = NULL,
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
p.bio.spp

##########################################################################
## Prop N in structure - soil N
##########################################################################
p.structure <- ggplot(data = data, aes(x = soil.n.norm, y = p.structure)) +
  geom_jitter(aes(fill = nrcs.code, shape = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, 0.03)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24), 
                     labels = c("AS" = "+N; +S",
                                "NO3" = "+N; -S",
                                "S" = "-N; +S",
                                "C" = "-N; -S")) +
  labs(x = expression(bold("Soil N (μg N g"^"-1"*" resin d"^"-1"*")")),
       y = expression(bold(rho["str"]*" (gN"["str"]* " gN"["leaf"]*""^-1*")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
p.structure

##########################################################################
## Prop N in structure - SPECIES
##########################################################################
p.str.spp <- ggplot(data = subset(spp.data, variable == "p.structure"),
                    aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, 
              aes(x = nrcs.code, y = p.structure, fill = nrcs.code,
                  shape = treatment),
              size = 4, width = 0.2, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 0.12, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, 0.03)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA" = expression(italic("A. saccharum")),
                               "FRAM" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  labs(x = "Species",
       y = NULL,
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
p.str.spp

##########################################################################
## Figure 1: leaf N
##########################################################################
png("../working_drafts/figs/NxS_fig1_leafn.png",
    width = 12, height = 12, units = 'in', res = 600)
ggarrange(narea, narea.spp,
          nleaf, nleaf.spp,
          marea, marea.spp, ncol = 2, nrow = 3, align = "hv",
          common.legend = TRUE, legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold"))
dev.off()

##########################################################################
## Figure 2: Leaf biochemistry
##########################################################################
png("../working_drafts/figs/NxS_fig2_leafbiochem.png",
    width = 16, height = 12, units = 'in', res = 600)
ggarrange(a400, a400.spp, a400.narea,
          vcmax, vcmax.spp, vcmax.narea,
          jmax, jmax.spp, jmax.narea,
          ncol = 3, nrow = 3, align = "hv",
          common.legend = TRUE, legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold"))
dev.off()

##########################################################################
## Figure S2: leaf N allocation
##########################################################################
png("../working_drafts/figs/NxS_figS2_leafn_allocation.png",
    width = 12, height = 16, units = 'in', res = 600)
ggarrange(p.rubisco, p.rubisco.spp,
          p.bioenergetics, p.bio.spp,
          p.photo, p.photo.spp,
          p.structure, p.str.spp,
          ncol = 2, nrow = 4, align = "hv",
          common.legend = TRUE, legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold"))
dev.off()

##########################################################################
## Figure 3: PNUE/iWUE 
##########################################################################
png("../working_drafts/figs/NxS_fig3_pnueiwue.png",
    width = 16, height = 16, units = 'in', res = 600)
ggarrange(chi, chi.spp, chi.narea, 
          pnue.plot, pnue.spp, blank.plot,
          narea.chi.plot, narea.chi.spp, blank.plot,
          vcmax.chi, vcmax.chi.spp, vcmax.chi.narea,
          ncol = 3, nrow = 4, align = "hv",
          common.legend = TRUE, legend = "right", 
          labels = c("A", "B", "C", "D", "E", NA, "F", "G", NA,
                     "H", "I", "J"),
          font.label = list(size = 18, face = "bold"))
dev.off()
