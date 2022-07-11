## Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)

## Central figure theme
pubtheme <- theme_bw(base_size = 18) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
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
data <- read.csv("../data_sheets/NxS_datasheet.csv",
                 stringsAsFactors = FALSE,
                 na.strings = "NA")
spp.data <- read.csv("../data_sheets/NxS_figs_emmeanOutputs.csv",
                     strip.white = TRUE)
spp.data$.group <- tolower(spp.data$.group)

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

## Reorder treatments
data$treatment <- factor(data$treatment, levels = c("NO3", "AS", "S", "C"))

##########################################################################
## Nleaf
##########################################################################
nleaf.spp <- ggplot(data = subset(spp.data, 
                                  variable == "leaf.n"),
                    aes(x = nrcs.code, 
                        y = emmean, 
                        fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = leaf.n, fill = nrcs.code,
                               shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax =upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 4, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(1, 4), breaks = seq(1, 4, 1)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = NULL,
       y = expression(bold("N"["mass"]~"(g g"^"-1"~")")),
       shape = "Treatment",
       fill = "Species") +
  pubtheme
nleaf.spp

##########################################################################
## SLA
##########################################################################
sla.spp <- ggplot(data = subset(spp.data, variable == "sla"),
                  aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = sla, fill = nrcs.code,
                               shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 500, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 125)) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  labs(x = "Species",
       y = expression(bold("SLA (cm"^"2"~"g"^"-1"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
sla.spp

##########################################################################
## Narea
##########################################################################
narea.spp <- ggplot(data = subset(spp.data, variable == "narea"),
                    aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = narea, fill = nrcs.code,
                               shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 3, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  labs(x = "Species",
       y = expression(bold("N"["area"]~"(g m"^"-2"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
narea.spp

##########################################################################
## Net photosynthesis (area basis)
##########################################################################
a.area.spp <- ggplot(data = subset(spp.data, variable == "a400"),
                     aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = a400, fill = nrcs.code,
                               shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 16, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, 4)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  labs(x = NULL,
       y = expression(bold("A"["net"]~ "(μmol m"^"-2"~"s"^"-1"~")")),
       shape = "Treatment",
       fill = "Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
a.area.spp

##########################################################################
## Vcmax.area standardized to 25degC 
##########################################################################
vcmax.spp <- ggplot(data = subset(spp.data, variable == "vcmax25"),
                    aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, 
              aes(x = nrcs.code, y = vcmax25, fill = nrcs.code, shape = treatment),
              size = 3, width = 0.1, alpha = 0.9) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 120, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  labs(x = NULL,
       y = expression(bold("V"["cmax25"]~"(μmol m"^"-2" ~ "s"^"-1" ~")")),
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.spp

##########################################################################
## Jmax.area standardized to 25degC 
##########################################################################
jmax.spp <- ggplot(data = subset(spp.data, variable == "jmax25"),
                   aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = jmax25, fill = nrcs.code,
                               shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 200, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 200), 
                     breaks = seq(0, 200, 50)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  labs(x = NULL,
       y = expression(bold("J"["max25"]~"(μmol m"^"-2" ~ "s"^"-1" ~")")),
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + 
  pubtheme
jmax.spp

##########################################################################
## Jmax:Vcmax standardized to 25degC 
##########################################################################
jmax.vcmax.spp <- ggplot(data = subset(spp.data, variable == "jmax:vcmax"),
                         aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = jmax.vcmax,
                               fill = nrcs.code,
                               shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 3, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(1, 3), breaks = seq(1, 3, 0.5)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  labs(x = NULL,
       y = expression(bold("J"["max25"]~":V"["cmax25"])),
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
jmax.vcmax.spp

##########################################################################
## Stomatal conductance
##########################################################################
gs.spp <- ggplot(data = subset(spp.data, variable == "gs"),
                 aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = gsw, fill = nrcs.code,
                               shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 0.20, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 0.20), breaks = seq(0, 0.20, 0.05)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  labs(x = NULL,
       y = expression(bold("g"["s"]~"(mol m"^"-2"~"s"^"-1"~")")),
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
gs.spp

##########################################################################
## chi
##########################################################################
chi.spp <- ggplot(data = subset(spp.data, variable == "chi"),
                  aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data,
              aes(x = nrcs.code, y = chi, fill = nrcs.code, shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 1, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.2)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  labs(x = NULL,
       y = expression(bold("χ (Pa Pa"^"-1"~")")),
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
chi.spp

##########################################################################
## Stomatal limitation
##########################################################################
stomlim.spp <- ggplot(data = subset(spp.data, variable == "stom.lim"),
                      aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, 
              aes(x = nrcs.code, y = stom.lim, fill = nrcs.code, shape = treatment),
              size = 3, width = 0.1, alpha = 0.9) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 1.05, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0.75, 1.05), breaks = seq(0.75, 1.05, 0.1)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  labs(x = NULL,
       y = expression(bold("Stomatal limitation")),
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
stomlim.spp

##########################################################################
## PNUE
##########################################################################
pnue.spp <- ggplot(data = subset(spp.data, variable == "pnue"),
                   aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = pnue, fill = nrcs.code,
                               shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 9, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, 3)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  labs(x = NULL,
       y = expression(bold("PNUE (μmol CO"["2"]~" gN"^"-1"~"s"^"-1"~")")),
       fill = "Species",
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme +
  theme(axis.title.y = element_text(size = 15))
pnue.spp

##########################################################################
## Narea-gs
##########################################################################
narea.gs.spp <- ggplot(data = subset(spp.data, variable == "narea.gs"),
                       aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data,
              aes(x = nrcs.code, y = narea.gs, fill = nrcs.code, shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 200, label = .group), 
            fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = NULL,
       y = expression(bold("N"["area"]~":g"["s"]~"(gN s mol"^"-1"~"H"["2"]~"O)")),
       fill = "Species",
       shape = "Treatment") +
  pubtheme  +
  theme(axis.title.y = element_text(size = 15))
narea.gs.spp

##########################################################################
## vcmax-gs
##########################################################################
vcmax.gs.spp <- ggplot(data = subset(spp.data, variable == "vcmax.gs"),
                       aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = vcmax.gs, fill = nrcs.code,
                               shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 8000, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 8000), breaks = seq(0, 8000, 2000)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = NULL,
       y = expression(bold("V"["cmax"]~":g"["s"]~"(μmol CO"["2"]~"mol"^"-1"~"H"["2"]~"O)")),
       fill = "Species",
       shape = "Treatment") +
  pubtheme +
  theme(axis.title.y = element_text(size = 15))
vcmax.gs.spp

##########################################################################
## Change in basal area
##########################################################################
ba.spp <- ggplot(data = subset(spp.data, variable == "basal.area"),
                 aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = ba.2011.2019,
                               fill = nrcs.code, shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 140, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 140), breaks = seq(0, 140, 35)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = NULL,
       y = expression(bold("Basal area (cm"^"2"~"yr"^"-1"~")")),
       fill = "Species",
       shape = "Treatment") +
  pubtheme
ba.spp

##########################################################################
## Relative growth rate
##########################################################################
rgr.spp <- ggplot(data = subset(spp.data, variable == "rgr"),
                  aes(x = nrcs.code, y = emmean, fill = nrcs.code)) +
  geom_jitter(data = data, aes(x = nrcs.code, y = growth.2011.2019,
                               fill = nrcs.code, shape = treatment),
              size = 3, width = 0.1, alpha = 0.75) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 1, width = 0.5) +
  geom_point(size = 5, shape = 21, fill = "black") +
  geom_text(aes(y = 200, label = .group), fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  scale_shape_manual(values = c(21, 22, 23, 25), labels = c("AS" = "+N, +S",
                                                            "NO3" = "+N, -S",
                                                            "S" = "-N, +S",
                                                            "C" = "-N, -S")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = NULL,
       y = expression(bold("RGR (g kg"^"-1"~"yr"^"-1"~")")),
       fill = "Species",
       shape = "Treatment") +
  pubtheme
rgr.spp

##########################################################################
## Figure S2: leaf N allocation
##########################################################################
png("../working_drafts/figs/NxS_figS2_leafn.png",
    width = 12, height = 6, units = 'in', res = 600)
ggarrange(narea.spp, ggarrange(nleaf.spp, sla.spp, ncol = 1, nrow = 2,
                           legend = "none", align = "hv", labels = c("B", "C"),
                           font.label = list(size = 18, face = "bold")),
          common.legend = TRUE, legend = "right", labels = c("A"),
          widths = c(1.5, 1),
          font.label = list(size = 18, face = "bold"))
dev.off()

##########################################################################
## Figure S3: Leaf biochemistry
##########################################################################
png("../working_drafts/figs/NxS_figS3_leafbiochem.png",
    width = 10, height = 8, units = 'in', res = 600)
ggarrange(a.area.spp, vcmax.spp, jmax.spp, jmax.vcmax.spp,
          ncol = 2, nrow = 2, common.legend = TRUE,
          align = "hv", legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold"))
dev.off()

##########################################################################
## Figure S4: PNUE/iWUE tradeoffs
##########################################################################
png("../working_drafts/figs/NxS_figS4_pnueiwue.png",
    width = 10, height = 8, units = 'in', res = 600)
ggarrange(pnue.spp, chi.spp, narea.gs.spp, vcmax.gs.spp, ncol = 2, nrow = 2, 
          common.legend = TRUE, align = "hv", legend = "right", 
          labels = "AUTO", font.label = list(size = 18, face = "bold"))
dev.off()

##########################################################################
## Figure S5: Whole plant measurements
##########################################################################
png("../working_drafts/figs/NxS_figS5_wp.png",
    width = 10, height = 4.5, units = 'in', res = 600)
ggarrange(ba.spp, rgr.spp, ncol = 2, nrow = 1, common.legend = TRUE,
          align = "hv", legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold"))
dev.off()

