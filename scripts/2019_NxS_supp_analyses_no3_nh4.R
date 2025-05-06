##########################################################################
## Load libraries and import data
##########################################################################
# Load libraries
library(lme4)
library(emmeans)
library(car)
library(tidyverse)
library(MuMIn)
library(ggpubr)

# Remove emmeans digit limiter
emm_options(opt.digits = FALSE)

# Import datasheet
data <- read.csv("../data/2019_NxS_datasheet.csv", 
                 stringsAsFactors = FALSE,
                 na.strings = "NA") %>%
  mutate(treatment = factor(treatment, 
                            levels = c("C", "NO3", "AS", "S")))
head(data)

# Subset dataset to only include Acer saccharum
plot_data <- subset(data, nrcs.code == "ACSA3")

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

##########################################################################
##########################################################################
## Nitrate-specific analyses
##########################################################################
##########################################################################

##########################################################################
## Nleaf - soil NO3-N
##########################################################################
nmass_no3n <- lmer(leaf.n ~ soil.no3n.norm + mineral.pH + (1 | site), 
               data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(nmass_no3n)
qqnorm(residuals(nmass_no3n))
qqline(residuals(nmass_no3n))
hist(residuals(nmass_no3n))
shapiro.test(residuals(nmass_no3n))
outlierTest(nmass_no3n)

# Model output
summary(nmass_no3n)
Anova(nmass_no3n)
r.squaredGLMM(nmass_no3n)

# Post-hoc tests
test(emtrends(nmass_no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
nmass.no3n.trend <- data.frame(
  emmeans(nmass_no3n, ~1, "soil.no3n.norm", 
          at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
nmass.no3n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.no3n.norm, y = leaf.n/100)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = nmass.no3n.trend, 
              aes(y = emmean/100, ymin = lower.CL/100, ymax = upper.CL/100),
              alpha = 0.2) +
  geom_smooth(data = nmass.no3n.trend, 
              aes(y = emmean/100), 
              size = 2, se = FALSE, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 25), 
                     breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0.016, 0.032), 
                     breaks = seq(0.016, 0.032, 0.004)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"["mass"]*" (gN g"["dry_mass"]*""^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
nmass.no3n.plot

##########################################################################
## Leaf mass per area - soil NO3-N
##########################################################################
marea_no3n <- lmer(marea ~ soil.no3n.norm + mineral.pH + (1 | site),
              data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(marea_no3n)
qqnorm(residuals(marea_no3n))
qqline(residuals(marea_no3n))
hist(residuals(marea_no3n))
shapiro.test(residuals(marea_no3n))
outlierTest(marea_no3n)

# Model output
summary(marea_no3n)
Anova(marea_no3n)
r.squaredGLMM(marea_no3n)

# Plot
marea.no3n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.no3n.norm, y = marea)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 25), 
                     breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("M"[area]*" (g"["dry_mass"]*" m"["leaf"]*""^"-2"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
marea.no3n.plot

##########################################################################
## Narea - soil NO3-N
##########################################################################
narea_no3n <- lmer(narea ~ soil.no3n.norm + mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea_no3n)
qqnorm(residuals(narea_no3n))
qqline(residuals(narea_no3n))
hist(residuals(narea_no3n))
shapiro.test(residuals(narea_no3n))
outlierTest(narea_no3n)

# Model output
summary(narea_no3n)
Anova(narea_no3n)
r.squaredGLMM(narea_no3n)

# Pairwise comparisons
test(emtrends(narea_no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
narea.no3n.trend <- data.frame(
  emmeans(narea_no3n, ~1, "soil.no3n.norm", 
          at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
narea.no3n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.no3n.norm, y = narea)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = narea.no3n.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = narea.no3n.trend, 
              aes(y = emmean), 
              size = 2, se = FALSE, color = "black") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
narea.no3n.plot

##########################################################################
## Anet - soil NO3-N
##########################################################################
a400_no3n <- lmer(sqrt(a400) ~ soil.no3n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(a400_no3n)
qqnorm(residuals(a400_no3n))
qqline(residuals(a400_no3n))
hist(residuals(a400_no3n))
shapiro.test(residuals(a400_no3n))
outlierTest(a400_no3n)

# Model output
summary(a400_no3n)
Anova(a400_no3n)
r.squaredGLMM(a400_no3n)

# Plot
a400.no3n.plot <- ggplot(data = plot_data, 
                         aes(x = soil.no3n.norm, y = a400)) +
  geom_jitter(aes(shape = treatment, fill = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("A"["net"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
a400.no3n.plot

##########################################################################
## gsw - soil NO3-N
##########################################################################
gs_no3n <- lmer(log(gsw) ~ soil.no3n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(gs_no3n)
qqnorm(residuals(gs_no3n))
qqline(residuals(gs_no3n))
hist(residuals(gs_no3n))
shapiro.test(residuals(gs_no3n))
outlierTest(gs_no3n)

# Model output
summary(gs_no3n)
Anova(gs_no3n)
r.squaredGLMM(gs_no3n)

# Plot
gs.no3n.plot <- ggplot(data = plot_data, 
                       aes(x = soil.no3n.norm, y = a400)) +
  geom_jitter(aes(shape = treatment, fill = treatment), 
              width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("A"["net"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
a400.no3n.plot

##########################################################################
## Vcmax25 - soil NO3-N
##########################################################################
vcmax_no3n <- lmer(vcmax25 ~ soil.no3n.norm + mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(vcmax_no3n)
qqnorm(residuals(vcmax_no3n))
qqline(residuals(vcmax_no3n))
hist(residuals(vcmax_no3n))
shapiro.test(residuals(vcmax_no3n))
outlierTest(vcmax_no3n)

# Model output
summary(vcmax_no3n)
Anova(vcmax_no3n)
r.squaredGLMM(vcmax_no3n)

# Pairwise comparisons
test(emtrends(vcmax_no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
vcmax.no3n.trend <- data.frame(
  emmeans(vcmax_no3n, ~1, "soil.no3n.norm",
          at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
vcmax.no3n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.no3n.norm, y = vcmax25)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = vcmax.no3n.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = vcmax.no3n.trend, 
              aes(y = emmean), 
              size = 2, se = FALSE, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("V"["cmax25"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
vcmax.no3n.plot

##########################################################################
## Jmax25 - soil NO3-N
##########################################################################
jmax_no3n <- lmer(jmax25 ~ soil.no3n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(jmax_no3n)
qqnorm(residuals(jmax_no3n))
qqline(residuals(jmax_no3n))
hist(residuals(jmax_no3n))
shapiro.test(residuals(jmax_no3n))
outlierTest(jmax_no3n)

# Model output
summary(jmax_no3n)
Anova(jmax_no3n)
r.squaredGLMM(jmax_no3n)

## Pairwise comparisons
test(emtrends(jmax_no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
jmax.no3n.trend <- data.frame(
  emmeans(jmax_no3n, ~1, "soil.no3n.norm",
          at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
jmax.no3n.plot <- ggplot(data = plot_data, 
                         aes(x = soil.no3n.norm, y = jmax25)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = jmax.no3n.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = jmax.no3n.trend, 
              aes(y = emmean), 
              size = 2, se = FALSE, color = "black") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("J"["max25"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment",
       fill = "Treatment") +
  pubtheme
jmax.no3n.plot


##########################################################################
## Jmax25:Vcmax25 - soil NO3-N
##########################################################################
vjmax_no3n <- lmer(jmax.vcmax ~ soil.no3n.norm + mineral.pH + (1 | site),
                   data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vjmax_no3n)
qqnorm(residuals(vjmax_no3n))
qqline(residuals(vjmax_no3n))
hist(residuals(vjmax_no3n))
shapiro.test(residuals(vjmax_no3n))
outlierTest(vjmax_no3n)

# Model output
summary(vjmax_no3n)
Anova(vjmax_no3n)
r.squaredGLMM(vjmax_no3n)

# Plot
jvmax.no3n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.no3n.norm, y = jmax.vcmax)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(1, 3), breaks = seq(1, 3, 0.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("J"["max25"]*": V"["cmax25"]*" (unitless)")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
jvmax.no3n.plot

##########################################################################
## chi - soil NO3-N
##########################################################################
chi_no3n <- lmer(chi ~ soil.no3n.norm + mineral.pH + (1 | site),
                 data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(chi_no3n)
qqnorm(residuals(chi_no3n))
qqline(residuals(chi_no3n))
hist(residuals(chi_no3n))
shapiro.test(residuals(chi_no3n))
outlierTest(chi_no3n)

# Model output
summary(chi_no3n)
Anova(chi_no3n)
r.squaredGLMM(chi_no3n)

## Pairwise comparisons
test(emtrends(chi_no3n,  ~1, var = "soil.no3n.norm"))

# Plot prep
chi.no3.trend <- data.frame(emmeans(chi_no3n, ~1, "soil.no3n.norm", 
                                    at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
chi.no3n.plot <- ggplot(data = plot_data, 
                        aes(x = soil.no3n.norm, y = chi)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = chi.no3.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = chi.no3.trend, 
              aes(y = emmean), 
              size = 2, se = FALSE, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold(chi*" (unitless)")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
chi.no3n.plot

##########################################################################
## PNUE - soil NO3-N
##########################################################################
pnue_no3n <- lmer(pnue ~ soil.no3n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue_no3n)
qqnorm(residuals(pnue_no3n))
qqline(residuals(pnue_no3n))
hist(residuals(pnue_no3n))
shapiro.test(residuals(pnue_no3n))
outlierTest(pnue_no3n)

# Model output
summary(pnue_no3n)
Anova(pnue_no3n)
r.squaredGLMM(pnue_no3n)

# Pairwise comparisons
test(emtrends(pnue_no3n, ~1, var = "soil.no3n.norm"))

# Plot
pnue.no3n.plot <- ggplot(data = plot_data, 
                         aes(x = soil.no3n.norm, y = pnue)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("PNUE (μmol CO"["2"]*" mol"^"-1"*"N s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme +
  theme(axis.title.y = element_text(size = 13))
pnue.no3n.plot

##########################################################################
## Narea.chi - soil NO3-N
##########################################################################
narea.chi_no3n <- lmer(narea.chi ~ soil.no3n.norm + mineral.pH + (1 | site),
                       data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea.chi_no3n)
qqnorm(residuals(narea.chi_no3n))
qqline(residuals(narea.chi_no3n))
hist(residuals(narea.chi_no3n))
shapiro.test(residuals(narea.chi_no3n))
outlierTest(narea.chi_no3n)

# Model output
summary(narea.chi_no3n)
Anova(narea.chi_no3n)
r.squaredGLMM(narea.chi_no3n)

# Post-hoc tests
test(emtrends(narea.chi_no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
narea.chi.no3n.trend <- data.frame(
  emmeans(narea.chi_no3n, ~1, "soil.no3n.norm", 
          at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
narea.chi.no3n.plot <- ggplot(data = plot_data, 
                              aes(x = soil.no3n.norm, y = narea.chi)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = narea.chi.no3n.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = narea.chi.no3n.trend, 
              aes(y = emmean), 
              size = 2, se = FALSE, color = "black") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"["area"]*":"*chi*" (gN m"["leaf"]*""^"-2"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
narea.chi.no3n.plot

##########################################################################
## Narea:gs - soil NO3-N
##########################################################################
narea.gs_no3n <- lmer(log(narea.gs) ~ soil.no3n.norm + mineral.pH + (1 | site),
                       data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea.gs_no3n)
qqnorm(residuals(narea.gs_no3n))
qqline(residuals(narea.gs_no3n))
hist(residuals(narea.gs_no3n))
shapiro.test(residuals(narea.gs_no3n))
outlierTest(narea.gs_no3n)

# Model output
summary(narea.gs_no3n)
Anova(narea.gs_no3n)
r.squaredGLMM(narea.gs_no3n)

# Post-hoc tests
test(emtrends(narea.gs_no3n, ~1, var = "soil.no3n.norm"))

# Plot
narea.gs.no3n.plot <- ggplot(data = plot_data, 
                             aes(x = soil.no3n.norm, y = narea.gs)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"["area"]*": g"["s"]*" (gN m"["leaf"]*""^"-2"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
narea.gs.no3n.plot

##########################################################################
## Vcmax25.chi - soil NO3-N
##########################################################################
data$vcmax.chi[85] <- NA

vcmax.chi_no3n <- lmer(vcmax.chi ~ soil.no3n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vcmax.chi_no3n)
qqnorm(residuals(vcmax.chi_no3n))
qqline(residuals(vcmax.chi_no3n))
hist(residuals(vcmax.chi_no3n))
shapiro.test(residuals(vcmax.chi_no3n))
outlierTest(vcmax.chi_no3n)

# Model output
summary(vcmax.chi_no3n)
Anova(vcmax.chi_no3n)
r.squaredGLMM(vcmax.chi_no3n)

# Post-hoc tests
test(emtrends(vcmax.chi_no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
vcmax.chi.no3n.trend <- data.frame(
  emmeans(vcmax.chi_no3n, ~1, "soil.no3n.norm", 
          at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
vcmax.chi.no3n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.no3n.norm, y = vcmax.chi)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = vcmax.chi.no3n.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = vcmax.chi.no3n.trend, 
              aes(y = emmean), 
              size = 2, se = FALSE, color = "black") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(-8, 120), breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("V"["cmax25"]*":"*chi*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
vcmax.chi.no3n.plot

##########################################################################
## Vcmax25:gs - soil NO3-N
##########################################################################
vcmax.gs_no3n <- lmer(log(vcmax.gs) ~ soil.no3n.norm + mineral.pH + (1 | site),
                      data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vcmax.gs_no3n)
qqnorm(residuals(vcmax.gs_no3n))
qqline(residuals(vcmax.gs_no3n))
hist(residuals(vcmax.gs_no3n))
shapiro.test(residuals(vcmax.gs_no3n))
outlierTest(vcmax.gs_no3n)

# Model output
summary(vcmax.gs_no3n)
Anova(vcmax.gs_no3n)
r.squaredGLMM(vcmax.gs_no3n)

# Post-hoc tests
test(emtrends(vcmax.gs_no3n, ~1, var = "soil.no3n.norm"))

##########################################################################
## PNUE.chi - soil NO3-N
##########################################################################
pnue.chi_no3n <- lmer(pnue.chi ~ soil.no3n.norm + mineral.pH + (1 | site),
                      data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue.chi_no3n)
qqnorm(residuals(pnue.chi_no3n))
qqline(residuals(pnue.chi_no3n))
hist(residuals(pnue.chi_no3n))
shapiro.test(residuals(pnue.chi_no3n))
outlierTest(pnue.chi_no3n)

# Model output
summary(pnue.chi_no3n)
Anova(pnue.chi_no3n)
r.squaredGLMM(pnue.chi_no3n)

# Post-hoc tests
test(emtrends(pnue.chi_no3n, ~1, var = "soil.no3n.norm"))

# Plot
pnue.chi.no3n.plot <- ggplot(data = plot_data, 
                             aes(x = soil.no3n.norm, y = pnue.chi)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 130), breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("PNUE:"*chi*" (μmol CO"["2"]*" mol N"^"-1"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme +
  theme(axis.title.y = element_text(size = 13))
pnue.chi.no3n.plot

##########################################################################
## PNUE.chi - soil NO3-N
##########################################################################
pnue.gs_no3n <- lmer(pnue.gs ~ soil.no3n.norm + mineral.pH + (1 | site),
                     data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue.gs_no3n)
qqnorm(residuals(pnue.gs_no3n))
qqline(residuals(pnue.gs_no3n))
hist(residuals(pnue.gs_no3n))
shapiro.test(residuals(pnue.gs_no3n))
outlierTest(pnue.gs_no3n)

# Model output
summary(pnue.gs_no3n)
Anova(pnue.gs_no3n)
r.squaredGLMM(pnue.gs_no3n)

# Post-hoc tests
test(emtrends(pnue.gs_no3n, ~1, var = "soil.no3n.norm"))

# Plot
pnue.chi.no3n.plot <- ggplot(data = plot_data, 
                             aes(x = soil.no3n.norm, y = pnue.chi)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 130), breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("PNUE:"*chi*" (μmol CO"["2"]*" mol N"^"-1"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme +
  theme(axis.title.y = element_text(size = 13))
pnue.chi.no3n.plot

##########################################################################
##########################################################################
## Ammonium-specific analyses
##########################################################################
##########################################################################

##########################################################################
## Nleaf
##########################################################################
nmass_nh4n <- lmer(leaf.n ~ soil.nh4n.norm + mineral.pH + (1 | site), 
               data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(nmass_nh4n)
qqnorm(residuals(nmass_nh4n))
qqline(residuals(nmass_nh4n))
hist(residuals(nmass_nh4n))
shapiro.test(residuals(nmass_nh4n))
outlierTest(nmass_nh4n)

# Model output
summary(nmass_nh4n)
Anova(nmass_nh4n)
r.squaredGLMM(nmass_nh4n)

# Post-hoc tests
test(emtrends(nmass_nh4n, ~1, var = "soil.nh4n.norm"))

# Plot
nmass.nh4n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.nh4n.norm, y = leaf.n/100)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0.016, 0.032), 
                     breaks = seq(0.016, 0.032, 0.004)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"["mass"]*" (gN g"["dry_mass"]*""^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
nmass.nh4n.plot

##########################################################################
## Leaf mass per area
##########################################################################
marea_nh4 <- lmer(log(marea) ~ soil.nh4n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(marea_nh4)
qqnorm(residuals(marea_nh4))
qqline(residuals(marea_nh4))
hist(residuals(marea_nh4))
shapiro.test(residuals(marea_nh4))
outlierTest(marea_nh4)

# Model output
summary(marea_nh4)
Anova(marea_nh4)
r.squaredGLMM(marea_nh4)

# Post-hoc tests
test(emtrends(marea_nh4, ~1, var = "soil.nh4n.norm"))

# Plot prep
marea.nh4n.trend <- data.frame(
  emmeans(marea_nh4, ~1, "soil.nh4n.norm", 
          at = list(soil.nh4n.norm = seq(0, 7, 0.1)),
          type = "response"))

# Plot
marea.nh4n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.nh4n.norm, y = marea)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = marea.nh4n.trend, 
              aes(y = response, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = marea.nh4n.trend, 
              aes(y = response), 
              size = 2, se = FALSE, color = "black") +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("M"[area]*" (g"["dry_mass"]*" m"["leaf"]*""^"-2"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
marea.nh4n.plot

##########################################################################
## Narea
##########################################################################
narea_nh4 <- lmer(narea ~ soil.nh4n.norm + mineral.pH + (1 | site), 
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea_nh4)
qqnorm(residuals(narea_nh4))
qqline(residuals(narea_nh4))
hist(residuals(narea_nh4))
shapiro.test(residuals(narea_nh4))
outlierTest(narea_nh4)

# Model output
summary(narea_nh4)
Anova(narea_nh4)
r.squaredGLMM(narea_nh4)

# Pairwise comparisons
test(emtrends(narea_nh4, ~1, var = "soil.nh4n.norm"))

# Plot prep
narea.nh4n.trend <- data.frame(
  emmeans(narea_nh4, ~1, "soil.nh4n.norm", 
          at = list(soil.nh4n.norm = seq(0, 7, 0.1))))

# Plot
narea.nh4n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.nh4n.norm, y = narea)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = narea.nh4n.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = narea.nh4n.trend, 
              aes(y = emmean), 
              size = 2, se = FALSE, color = "black") +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"[area]*" (g N m"["leaf"]*""^"-2"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
narea.nh4n.plot

##########################################################################
## Anet
##########################################################################
a400_nh4 <- lmer(sqrt(a400) ~ soil.nh4n.norm + mineral.pH + (1 | site),
                 data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(a400_nh4)
qqnorm(residuals(a400_nh4))
qqline(residuals(a400_nh4))
hist(residuals(a400_nh4))
shapiro.test(residuals(a400_nh4))
outlierTest(a400_nh4)

# Model output
summary(a400_nh4)
Anova(a400_nh4)
r.squaredGLMM(a400_nh4)

# Plot
a400.nh4n.plot <- ggplot(data = plot_data, 
                         aes(x = soil.nh4n.norm, y = a400)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("A"["net"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
a400.nh4n.plot

##########################################################################
## gs
##########################################################################
gs_nh4 <- lmer(log(gsw) ~ soil.nh4n.norm + mineral.pH + (1 | site),
                 data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(gs_nh4)
qqnorm(residuals(gs_nh4))
qqline(residuals(gs_nh4))
hist(residuals(gs_nh4))
shapiro.test(residuals(gs_nh4))
outlierTest(gs_nh4)

# Model output
summary(gs_nh4)
Anova(gs_nh4)
r.squaredGLMM(gs_nh4)

# Plot
gs.nh4n.plot <- ggplot(data = plot_data, 
                       aes(x = soil.nh4n.norm, y = gsw)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.025)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("A"["net"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
gs.nh4n.plot

##########################################################################
## Vcmax25
##########################################################################
vcmax_nh4 <- lmer(log(vcmax25) ~ soil.nh4n.norm + mineral.pH + (1 | site), 
                  data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(vcmax_nh4)
qqnorm(residuals(vcmax_nh4))
qqline(residuals(vcmax_nh4))
hist(residuals(vcmax_nh4))
shapiro.test(residuals(vcmax_nh4))
outlierTest(vcmax_nh4)

# Model output
summary(vcmax_nh4)
Anova(vcmax_nh4)
r.squaredGLMM(vcmax_nh4)

# Pairwise comparisons
test(emtrends(vcmax_nh4, ~1, var = "soil.nh4n.norm"))

# Plot
vcmax.nh4n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.nh4n.norm, y = vcmax25)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("V"["cmax25"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
vcmax.nh4n.plot

##########################################################################
## Jmax25
##########################################################################
jmax_nh4 <- lmer(jmax25 ~ soil.nh4n.norm + mineral.pH + (1 | site),
                 data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(jmax_nh4)
qqnorm(residuals(jmax_nh4))
qqline(residuals(jmax_nh4))
hist(residuals(jmax_nh4))
shapiro.test(residuals(jmax_nh4))
outlierTest(jmax_nh4)

# Model output
summary(jmax_nh4)
Anova(jmax_nh4)
r.squaredGLMM(jmax_nh4)

## Pairwise comparisons
test(emtrends(jmax, ~1, var = "soil.nh4n.norm"))

# Plot
jmax.nh4n.plot <- ggplot(data = plot_data, 
                         aes(x = soil.nh4n.norm, y = jmax25)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("J"["max25"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
jmax.nh4n.plot

##########################################################################
## Jmax25:Vcmax25
##########################################################################
vjmax_nh4 <- lmer(jmax.vcmax ~ soil.nh4n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vjmax_nh4)
qqnorm(residuals(vjmax_nh4))
qqline(residuals(vjmax_nh4))
hist(residuals(vjmax_nh4))
shapiro.test(residuals(vjmax_nh4))
outlierTest(vjmax_nh4)

# Model output
summary(vjmax_nh4)
Anova(vjmax_nh4)
r.squaredGLMM(vjmax_nh4)

# Plot
jvmax.nh4n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.nh4n.norm, y = jmax.vcmax)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(1, 3), breaks = seq(1, 3, 0.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("J"["max25"]*": V"["cmax25"]*" (unitless)")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
jvmax.nh4n.plot

##########################################################################
## chi
##########################################################################
chi_nh4 <- lmer(chi ~ soil.nh4n.norm + mineral.pH + (1 | site),
                data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(chi_nh4)
qqnorm(residuals(chi_nh4))
qqline(residuals(chi_nh4))
hist(residuals(chi_nh4))
shapiro.test(residuals(chi_nh4))
outlierTest(chi_nh4)

# Model output
summary(chi_nh4)
Anova(chi_nh4)
r.squaredGLMM(chi_nh4)

## Pairwise comparisons
test(emtrends(chi_nh4,  ~1, var = "soil.nh4n.norm"))

# Plot prep
chi.nh4.trend <- data.frame(
  emmeans(chi_nh4, ~1, "soil.nh4n.norm", 
          at = list(soil.nh4n.norm = seq(0, 7, 0.1))))

# Plot
chi.nh4n.plot <- ggplot(data = plot_data, 
                        aes(x = soil.nh4n.norm, y = chi)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = chi.nh4.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = chi.nh4.trend, 
              aes(y = emmean), 
              size = 2, se = FALSE, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold(chi*" (unitless)")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
chi.nh4n.plot

##########################################################################
## PNUE
##########################################################################
pnue_nh4 <- lmer(pnue ~ soil.nh4n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue_nh4)
qqnorm(residuals(pnue_nh4))
qqline(residuals(pnue_nh4))
hist(residuals(pnue_nh4))
shapiro.test(residuals(pnue_nh4))
outlierTest(pnue_nh4)

# Model output
summary(pnue_nh4)
Anova(pnue_nh4)
r.squaredGLMM(pnue_nh4)

# Pairwise comparisons
test(emtrends(pnue_nh4, ~1, var = "soil.nh4n.norm"))

# Plot prep
pnue.nh4.trend <- data.frame(
  emmeans(pnue_nh4, ~1, "soil.nh4n.norm", 
          at = list(soil.nh4n.norm = seq(0, 7, 0.1))))

# Plot
pnue.nh4n.plot <- ggplot(data = plot_data, 
                         aes(x = soil.nh4n.norm, y = pnue)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = pnue.nh4.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = pnue.nh4.trend, 
              aes(y = emmean), 
              size = 2, se = FALSE, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("PNUE (μmol CO"["2"]*" mol"^"-1"*"N s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme +
  theme(axis.title.y = element_text(size = 13))
pnue.nh4n.plot

##########################################################################
## Narea.chi
##########################################################################
narea.chi_nh4 <- lmer(narea.chi ~ soil.nh4n.norm + mineral.pH + (1 | site),
                      data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea.chi_nh4)
qqnorm(residuals(narea.chi_nh4))
qqline(residuals(narea.chi_nh4))
hist(residuals(narea.chi_nh4))
shapiro.test(residuals(narea.chi_nh4))
outlierTest(narea.chi_nh4)

# Model output
summary(narea.chi_nh4)
Anova(narea.chi_nh4)
r.squaredGLMM(narea.chi_nh4)

# Post-hoc tests
test(emtrends(narea.chi_nh4, ~1, var = "soil.nh4n.norm"))

# Plot prep
narea.chi.nh4n.trend <- data.frame(
  emmeans(narea.chi_nh4, ~1, "soil.nh4n.norm", 
          at = list(soil.nh4n.norm = seq(0, 7, 0.1))))

# Plot
narea.chi.nh4n.plot <- ggplot(data = plot_data, 
                              aes(x = soil.nh4n.norm, y = narea.chi)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = narea.chi.nh4n.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = narea.chi.nh4n.trend, 
              aes(y = emmean), 
              size = 2, se = FALSE, color = "black") +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"["area"]*":"*chi*" (gN m"["leaf"]*""^"-2"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
narea.chi.nh4n.plot

##########################################################################
## Narea:gs
##########################################################################
narea.gs_nh4 <- lmer(log(narea.gs) ~ soil.nh4n.norm + mineral.pH + (1 | site),
                      data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea.gs_nh4)
qqnorm(residuals(narea.gs_nh4))
qqline(residuals(narea.gs_nh4))
hist(residuals(narea.gs_nh4))
shapiro.test(residuals(narea.gs_nh4))
outlierTest(narea.gs_nh4)

# Model output
summary(narea.gs_nh4)
Anova(narea.gs_nh4)
r.squaredGLMM(narea.gs_nh4)

##########################################################################
## Vcmax25.chi
##########################################################################
vcmax.chi_nh4 <- lmer(vcmax.chi ~ soil.nh4n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vcmax.chi_nh4)
qqnorm(residuals(vcmax.chi_nh4))
qqline(residuals(vcmax.chi_nh4))
hist(residuals(vcmax.chi_nh4))
shapiro.test(residuals(vcmax.chi_nh4))
outlierTest(vcmax.chi_nh4)

# Model output
summary(vcmax.chi_nh4)
Anova(vcmax.chi_nh4)
r.squaredGLMM(vcmax.chi_nh4)

# Post-hoc tests
test(emtrends(vcmax.chi_nh4, ~1, var = "soil.nh4n.norm"))

# Plot
vcmax.chi.nh4n.plot <- ggplot(data = plot_data, 
                              aes(x = soil.nh4n.norm, y = vcmax.chi)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(-8, 120), breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("V"["cmax25"]*":"*chi*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme
vcmax.chi.nh4n.plot

##########################################################################
## Vcmax25:gs
##########################################################################
vcmax.gs_nh4 <- lmer(log(vcmax.gs) ~ soil.nh4n.norm + mineral.pH + (1 | site),
                      data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vcmax.gs_nh4)
qqnorm(residuals(vcmax.gs_nh4))
qqline(residuals(vcmax.gs_nh4))
hist(residuals(vcmax.gs_nh4))
shapiro.test(residuals(vcmax.gs_nh4))
outlierTest(vcmax.gs_nh4)

# Model output
summary(vcmax.gs_nh4)
Anova(vcmax.gs_nh4)
r.squaredGLMM(vcmax.gs_nh4)

# Post-hoc tests
test(emtrends(vcmax.gs_nh4, ~1, var = "soil.nh4n.norm"))

##########################################################################
## PNUE.chi
##########################################################################
pnue.chi_nh4 <- lmer(pnue.chi ~ soil.nh4n.norm + mineral.pH + (1 | site),
                     data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue.chi_nh4)
qqnorm(residuals(pnue.chi_nh4))
qqline(residuals(pnue.chi_nh4))
hist(residuals(pnue.chi_nh4))
shapiro.test(residuals(pnue.chi_nh4))
outlierTest(pnue.chi_nh4)

# Model output
summary(pnue.chi_nh4)
Anova(pnue.chi_nh4)
r.squaredGLMM(pnue.chi_nh4)

# Post-hoc tests
test(emtrends(pnue.chi_nh4, ~1, var = "soil.nh4n.norm"))

# Plot prep
pnue.chi.nh4n.trend <- data.frame(
  emmeans(pnue.chi_nh4, ~1, "soil.nh4n.norm", 
          at = list(soil.nh4n.norm = seq(0, 7, 0.1))))

# Plot
pnue.chi.nh4n.plot <- ggplot(data = plot_data, 
                              aes(x = soil.nh4n.norm, y = pnue.chi)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = pnue.chi.nh4n.trend, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_smooth(data = pnue.chi.nh4n.trend, 
              aes(y = emmean), 
              size = 2, se = FALSE, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0, 130), breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#009E73"),
                    labels = c("C" = "no N; no S",
                               "NO3" = "+ N; no S",
                               "AS" = "+ N; + S",
                               "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NH"["4"]*"-N (μg NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("PNUE:"*chi*" (μmol CO"["2"]*" mol N"^"-1"*" s"^"-1"*")")),
       shape = "Treatment", fill = "Treatment") +
  pubtheme +
  theme(axis.title.y = element_text(size = 13))
pnue.chi.nh4n.plot

##########################################################################
## PNUE:gs
##########################################################################
pnue.gs_nh4 <- lmer(pnue.gs ~ soil.nh4n.norm + mineral.pH + (1 | site),
                    data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue.gs_nh4)
qqnorm(residuals(pnue.gs_nh4))
qqline(residuals(pnue.gs_nh4))
hist(residuals(pnue.gs_nh4))
shapiro.test(residuals(pnue.gs_nh4))
outlierTest(pnue.gs_nh4)

# Model output
summary(pnue.gs_nh4)
Anova(pnue.gs_nh4)
r.squaredGLMM(pnue.gs_nh4)

# Post-hoc tests
test(emtrends(pnue.gs_nh4, ~1, var = "soil.nh4n.norm"))

##########################################################################
## Make plots
##########################################################################

# Fig S3 - leaf N
png("../../nitrogen_pH/working_drafts/figs/NxS_figS3_leafn_no3_nh4.png",
    width = 10.5, height = 12, units = 'in', res = 600)
ggarrange(narea.no3n.plot, narea.nh4n.plot,
          nmass.no3n.plot, nmass.nh4n.plot,
          marea.no3n.plot, marea.nh4n.plot,
          ncol = 2, nrow = 3, align = "hv",
          common.legend = TRUE, legend = "right",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
          font.label = list(size = 18, face = "bold"))
dev.off()


# Fig S4 - leaf gas exchange
png("../../nitrogen_pH/working_drafts/figs/NxS_figS4_leafbiochem_no3_nh4.png",
    width = 10.5, height = 12, units = 'in', res = 600)
ggarrange(a400.no3n.plot, a400.nh4n.plot,
          vcmax.no3n.plot, vcmax.nh4n.plot,
          jmax.no3n.plot, jmax.nh4n.plot,
          ncol = 2, nrow = 3, align = "hv",
          common.legend = TRUE, legend = "right", 
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
          font.label = list(size = 18, face = "bold"))
dev.off()

# Fig S5 - nitrogen-water use tradeoffs
png("../../nitrogen_pH/working_drafts/figs/NxS_figS5_pnueiwue.png",
    width = 10.5, height = 16, units = 'in', res = 600)
ggarrange(chi.no3n.plot, chi.nh4n.plot,
          pnue.no3n.plot, pnue.nh4n.plot,
          narea.chi.no3n.plot, narea.chi.nh4n.plot,
          vcmax.chi.no3n.plot, vcmax.chi.nh4n.plot,
          ncol = 2, nrow = 4, align = "hv",
          common.legend = TRUE, legend = "right", 
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"),
          font.label = list(size = 18, face = "bold"))
dev.off()
