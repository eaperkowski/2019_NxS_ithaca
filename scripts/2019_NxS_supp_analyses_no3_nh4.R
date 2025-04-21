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
data <- read.csv("../data/2019_NxS_datasheet.csv", stringsAsFactors = FALSE,
                 na.strings = "NA") %>%
  mutate(treatment = factor(treatment, levels = c("C", "NO3", "AS", "S")))
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
leaf.n.no3n <- lmer(leaf.n ~ soil.no3n.norm + mineral.pH + (1 | site), 
               data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(leaf.n.no3n)
qqnorm(residuals(leaf.n.no3n))
qqline(residuals(leaf.n.no3n))
hist(residuals(leaf.n.no3n))
shapiro.test(residuals(leaf.n.no3n))
outlierTest(leaf.n.no3n)

# Model output
summary(leaf.n.no3n)
Anova(leaf.n.no3n)
r.squaredGLMM(leaf.n.no3n)

# Post-hoc tests
test(emtrends(leaf.n.no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
nmass.no3n.trend <- data.frame(emmeans(leaf.n.no3n, ~1, "soil.no3n.norm", 
                                       at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
nmass.no3n.plot <- ggplot(data = plot_data, 
                          aes(x = soil.no3n.norm, y = leaf.n/100)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_ribbon(data = nmass.no3n.trend, 
              aes(y = emmean/100,
                  ymin = lower.CL/100, ymax = upper.CL/100),
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
marea.no3n <- lmer(marea ~ soil.no3n.norm + mineral.pH + (1 | site),
              data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(marea.no3n)
qqnorm(residuals(marea.no3n))
qqline(residuals(marea.no3n))
hist(residuals(marea.no3n))
shapiro.test(residuals(marea.no3n))
outlierTest(marea.no3n)

# Model output
summary(marea.no3n)
Anova(marea.no3n)
r.squaredGLMM(marea.no3n)

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
narea.no3n <- lmer(narea ~ soil.no3n.norm + mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea.no3n)
qqnorm(residuals(narea.no3n))
qqline(residuals(narea.no3n))
hist(residuals(narea.no3n))
shapiro.test(residuals(narea.no3n))
outlierTest(narea.no3n)

# Model output
summary(narea.no3n)
Anova(narea.no3n)
r.squaredGLMM(narea.no3n)

# Pairwise comparisons
test(emtrends(narea.no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
narea.no3n.trend <- data.frame(emmeans(narea.no3n, ~1, "soil.no3n.norm", 
                                       at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
narea.no3n.plot <- ggplot(data = plot_data, aes(x = soil.no3n.norm, y = narea)) +
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
data$a400[data$a400 < 0.2] <- NA

a400.no3n <- lmer(sqrt(a400) ~ soil.no3n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(a400.no3n)
qqnorm(residuals(a400.no3n))
qqline(residuals(a400.no3n))
hist(residuals(a400.no3n))
shapiro.test(residuals(a400.no3n))
outlierTest(a400.no3n)

# Model output
summary(a400.no3n)
Anova(a400.no3n)
r.squaredGLMM(a400.no3n)

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
gs.no3n <- lmer(log(gsw) ~ soil.no3n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(gs.no3n)
qqnorm(residuals(gs.no3n))
qqline(residuals(gs.no3n))
hist(residuals(gs.no3n))
shapiro.test(residuals(gs.no3n))
outlierTest(gs.no3n)

# Model output
summary(gs.no3n)
Anova(gs.no3n)
r.squaredGLMM(gs.no3n)

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
vcmax.no3n <- lmer(vcmax25 ~ soil.no3n.norm + mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(vcmax.no3n)
qqnorm(residuals(vcmax.no3n))
qqline(residuals(vcmax.no3n))
hist(residuals(vcmax.no3n))
shapiro.test(residuals(vcmax.no3n))
outlierTest(vcmax.no3n)

# Model output
summary(vcmax.no3n)
Anova(vcmax.no3n)
r.squaredGLMM(vcmax.no3n)

# Pairwise comparisons
test(emtrends(vcmax.no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
vcmax.no3n.trend <- data.frame(emmeans(vcmax.no3n, ~1, "soil.no3n.norm", 
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
jmax.no3n <- lmer(jmax25 ~ soil.no3n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(jmax.no3n)
qqnorm(residuals(jmax.no3n))
qqline(residuals(jmax.no3n))
hist(residuals(jmax.no3n))
shapiro.test(residuals(jmax.no3n))
outlierTest(jmax.no3n)

# Model output
summary(jmax.no3n)
Anova(jmax.no3n)
r.squaredGLMM(jmax.no3n)

## Pairwise comparisons
test(emtrends(jmax.no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
jmax.no3n.trend <- data.frame(emmeans(jmax.no3n, ~1, "soil.no3n.norm", 
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
       shape = "Treatment", fill = "Treatment") +
  pubtheme
jmax.no3n.plot


##########################################################################
## Jmax25:Vcmax25 - soil NO3-N
##########################################################################
vjmax.no3n <- lmer(jmax.vcmax ~ soil.no3n.norm + mineral.pH + (1 | site),
              data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vjmax.no3n)
qqnorm(residuals(vjmax.no3n))
qqline(residuals(vjmax.no3n))
hist(residuals(vjmax.no3n))
shapiro.test(residuals(vjmax.no3n))
outlierTest(vjmax.no3n)

# Model output
summary(vjmax.no3n)
Anova(vjmax.no3n)
r.squaredGLMM(vjmax.no3n)

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
chi.no3n <- lmer(chi ~ soil.no3n.norm + mineral.pH + (1 | site),
            data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(chi.no3n)
qqnorm(residuals(chi.no3n))
qqline(residuals(chi.no3n))
hist(residuals(chi.no3n))
shapiro.test(residuals(chi.no3n))
outlierTest(chi.no3n)

# Model output
summary(chi.no3n)
Anova(chi.no3n)
r.squaredGLMM(chi.no3n)

## Pairwise comparisons
test(emtrends(chi.no3n,  ~1, var = "soil.no3n.norm"))

# Plot prep
chi.no3.trend <- data.frame(emmeans(chi.no3n, ~1, "soil.no3n.norm", 
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
pnue.no3n <- lmer(pnue ~ soil.no3n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue.no3n)
qqnorm(residuals(pnue.no3n))
qqline(residuals(pnue.no3n))
hist(residuals(pnue.no3n))
shapiro.test(residuals(pnue.no3n))
outlierTest(pnue.no3n)

# Model output
summary(pnue.no3n)
Anova(pnue.no3n)
r.squaredGLMM(pnue.no3n)

# Pairwise comparisons
test(emtrends(pnue.no3n, ~1, var = "soil.no3n.norm"))

# Plot
pnue.no3n.plot <- ggplot(data = plot_data, aes(x = soil.no3n.norm, y = pnue)) +
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
narea.chi.no3n <- lmer(narea.chi ~ soil.no3n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea.chi.no3n)
qqnorm(residuals(narea.chi.no3n))
qqline(residuals(narea.chi.no3n))
hist(residuals(narea.chi.no3n))
shapiro.test(residuals(narea.chi.no3n))
outlierTest(narea.chi.no3n)

# Model output
summary(narea.chi.no3n)
Anova(narea.chi.no3n)
r.squaredGLMM(narea.chi.no3n)

# Post-hoc tests
test(emtrends(narea.chi.no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
narea.chi.no3n.trend <- data.frame(emmeans(narea.chi.no3n, ~1, "soil.no3n.norm", 
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
## Vcmax25.chi - soil NO3-N
##########################################################################
data$vcmax.chi[85] <- NA

vcmax.chi.no3n <- lmer(vcmax.chi ~ soil.no3n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vcmax.chi.no3n)
qqnorm(residuals(vcmax.chi.no3n))
qqline(residuals(vcmax.chi.no3n))
hist(residuals(vcmax.chi.no3n))
shapiro.test(residuals(vcmax.chi.no3n))
outlierTest(vcmax.chi.no3n)

# Model output
summary(vcmax.chi.no3n)
Anova(vcmax.chi.no3n)
r.squaredGLMM(vcmax.chi.no3n)

# Post-hoc tests
test(emtrends(vcmax.chi.no3n, ~1, var = "soil.no3n.norm"))

# Plot prep
vcmax.chi.no3n.trend <- data.frame(
  emmeans(vcmax.chi.no3n, ~1, "soil.no3n.norm", 
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
## PNUE.chi - soil NO3-N
##########################################################################
pnue.chi.no3n <- lmer(pnue.chi ~ soil.no3n.norm + mineral.pH + (1 | site),
                      data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue.chi.no3n)
qqnorm(residuals(pnue.chi.no3n))
qqline(residuals(pnue.chi.no3n))
hist(residuals(pnue.chi.no3n))
shapiro.test(residuals(pnue.chi.no3n))
outlierTest(pnue.chi.no3n)

# Model output
summary(pnue.chi.no3n)
Anova(pnue.chi.no3n)
r.squaredGLMM(pnue.chi.no3n)

# Post-hoc tests
test(emtrends(pnue.chi.no3n, ~1, var = "soil.no3n.norm"))

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
leaf.n.nh4n <- lmer(leaf.n ~ soil.nh4n.norm + mineral.pH + (1 | site), 
               data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(leaf.n.nh4n)
qqnorm(residuals(leaf.n.nh4n))
qqline(residuals(leaf.n.nh4n))
hist(residuals(leaf.n.nh4n))
shapiro.test(residuals(leaf.n.nh4n))
outlierTest(leaf.n.nh4n)

# Model output
summary(leaf.n.nh4n)
Anova(leaf.n.nh4n)
r.squaredGLMM(leaf.n.nh4n)

# Post-hoc tests
test(emtrends(leaf.n.nh4n, ~1, var = "soil.nh4n.norm"))

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
marea.nh4 <- lmer(log(marea) ~ soil.nh4n.norm + mineral.pH + (1 | site),
              data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(marea.nh4)
qqnorm(residuals(marea.nh4))
qqline(residuals(marea.nh4))
hist(residuals(marea.nh4))
shapiro.test(residuals(marea.nh4))
outlierTest(marea.nh4)

# Model output
summary(marea.nh4)
Anova(marea.nh4)
r.squaredGLMM(marea.nh4)

# Post-hoc tests
test(emtrends(marea.nh4, ~1, var = "soil.nh4n.norm"))

# Plot prep
marea.nh4n.trend <- data.frame(
  emmeans(marea.nh4, ~1, "soil.nh4n.norm", 
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
narea.nh4 <- lmer(narea ~ soil.nh4n.norm + mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea.nh4)
qqnorm(residuals(narea.nh4))
qqline(residuals(narea.nh4))
hist(residuals(narea.nh4))
shapiro.test(residuals(narea.nh4))
outlierTest(narea.nh4)

# Model output
summary(narea.nh4)
Anova(narea.nh4)
r.squaredGLMM(narea.nh4)

# Pairwise comparisons
test(emtrends(narea.nh4, ~1, var = "soil.nh4n.norm"))

# Plot prep
narea.nh4n.trend <- data.frame(
  emmeans(narea.nh4, ~1, "soil.nh4n.norm", 
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
a400.nh4 <- lmer(sqrt(a400) ~ soil.nh4n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(a400.nh4)
qqnorm(residuals(a400.nh4))
qqline(residuals(a400.nh4))
hist(residuals(a400.nh4))
shapiro.test(residuals(a400.nh4))
outlierTest(a400.nh4)

# Model output
summary(a400.nh4)
Anova(a400.nh4)
r.squaredGLMM(a400.nh4)

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
gs.nh4 <- lmer(log(gsw) ~ soil.nh4n.norm + mineral.pH + (1 | site),
                 data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(gs.nh4)
qqnorm(residuals(gs.nh4))
qqline(residuals(gs.nh4))
hist(residuals(gs.nh4))
shapiro.test(residuals(gs.nh4))
outlierTest(gs.nh4)

# Model output
summary(gs.nh4)
Anova(gs.nh4)
r.squaredGLMM(gs.nh4)

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
## Vcmax25
##########################################################################
vcmax.nh4 <- lmer(log(vcmax25) ~ soil.nh4n.norm + mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(vcmax.nh4)
qqnorm(residuals(vcmax.nh4))
qqline(residuals(vcmax.nh4))
hist(residuals(vcmax.nh4))
shapiro.test(residuals(vcmax.nh4))
outlierTest(vcmax.nh4)

# Model output
summary(vcmax.nh4)
Anova(vcmax.nh4)
r.squaredGLMM(vcmax.nh4)

# Pairwise comparisons
test(emtrends(vcmax.nh4, ~1, var = "soil.nh4n.norm"))

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
jmax.nh4 <- lmer(jmax25 ~ soil.nh4n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(jmax.nh4)
qqnorm(residuals(jmax.nh4))
qqline(residuals(jmax.nh4))
hist(residuals(jmax.nh4))
shapiro.test(residuals(jmax.nh4))
outlierTest(jmax.nh4)

# Model output
summary(jmax.nh4)
Anova(jmax.nh4)
r.squaredGLMM(jmax.nh4)

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
vjmax.nh4 <- lmer(jmax.vcmax ~ soil.nh4n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vjmax.nh4)
qqnorm(residuals(vjmax.nh4))
qqline(residuals(vjmax.nh4))
hist(residuals(vjmax.nh4))
shapiro.test(residuals(vjmax.nh4))
outlierTest(vjmax.nh4)

# Model output
summary(vjmax.nh4)
Anova(vjmax.nh4)
r.squaredGLMM(vjmax.nh4)

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
chi.nh4 <- lmer(chi ~ soil.nh4n.norm + mineral.pH + (1 | site),
                data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(chi.nh4)
qqnorm(residuals(chi.nh4))
qqline(residuals(chi.nh4))
hist(residuals(chi.nh4))
shapiro.test(residuals(chi.nh4))
outlierTest(chi.nh4)

# Model output
summary(chi.nh4)
Anova(chi.nh4)
r.squaredGLMM(chi.nh4)

## Pairwise comparisons
test(emtrends(chi.nh4,  ~1, var = "soil.nh4n.norm"))

# Plot prep
chi.nh4.trend <- data.frame(
  emmeans(chi.nh4, ~1, "soil.nh4n.norm", 
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
pnue.nh4 <- lmer(pnue ~ soil.nh4n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue.nh4)
qqnorm(residuals(pnue.nh4))
qqline(residuals(pnue.nh4))
hist(residuals(pnue.nh4))
shapiro.test(residuals(pnue.nh4))
outlierTest(pnue.nh4)

# Model output
summary(pnue.nh4)
Anova(pnue.nh4)
r.squaredGLMM(pnue.nh4)

# Pairwise comparisons
test(emtrends(pnue.nh4, ~1, var = "soil.nh4n.norm"))

# Plot prep
pnue.nh4.trend <- data.frame(
  emmeans(pnue.nh4, ~1, "soil.nh4n.norm", 
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
narea.chi.nh4 <- lmer(narea.chi ~ soil.nh4n.norm + mineral.pH + (1 | site),
                      data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea.chi.nh4)
qqnorm(residuals(narea.chi.nh4))
qqline(residuals(narea.chi.nh4))
hist(residuals(narea.chi.nh4))
shapiro.test(residuals(narea.chi.nh4))
outlierTest(narea.chi.nh4)

# Model output
summary(narea.chi.nh4)
Anova(narea.chi.nh4)
r.squaredGLMM(narea.chi.nh4)

# Post-hoc tests
test(emtrends(narea.chi.nh4, ~1, var = "soil.nh4n.norm"))

# Plot prep
narea.chi.nh4n.trend <- data.frame(
  emmeans(narea.chi.nh4, ~1, "soil.nh4n.norm", 
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
## Vcmax25.chi
##########################################################################
vcmax.chi.nh4 <- lmer(vcmax.chi ~ soil.nh4n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vcmax.chi.nh4)
qqnorm(residuals(vcmax.chi.nh4))
qqline(residuals(vcmax.chi.nh4))
hist(residuals(vcmax.chi.nh4))
shapiro.test(residuals(vcmax.chi.nh4))
outlierTest(vcmax.chi.nh4)

# Model output
summary(vcmax.chi.nh4)
Anova(vcmax.chi.nh4)
r.squaredGLMM(vcmax.chi.nh4)

# Post-hoc tests
test(emtrends(vcmax.chi.nh4, ~1, var = "soil.nh4n.norm"))

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
## PNUE.chi
##########################################################################
pnue.chi.nh4 <- lmer(pnue.chi ~ soil.nh4n.norm + mineral.pH + (1 | site),
                      data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue.chi.nh4)
qqnorm(residuals(pnue.chi.nh4))
qqline(residuals(pnue.chi.nh4))
hist(residuals(pnue.chi.nh4))
shapiro.test(residuals(pnue.chi.nh4))
outlierTest(pnue.chi.nh4)

# Model output
summary(pnue.chi.nh4)
Anova(pnue.chi.nh4)
r.squaredGLMM(pnue.chi.nh4)

# Post-hoc tests
test(emtrends(pnue.chi.nh4, ~1, var = "soil.nh4n.norm"))

# Plot prep
pnue.chi.nh4n.trend <- data.frame(
  emmeans(pnue.chi.nh4, ~1, "soil.nh4n.norm", 
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
