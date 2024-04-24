##########################################################################
## Load libraries and import data
##########################################################################
# Load libraries
library(lme4)
library(emmeans)
library(car)
library(tidyverse)
library(MuMIn)

# Remove emmeans digit limiter
emm_options(opt.digits = FALSE)

# Import datasheet
data <- read.csv("../data/2019_NxS_datasheet.csv", stringsAsFactors = FALSE,
                 na.strings = "NA") %>%
  mutate(anet.mass = a400 / marea,
         vcmax.mass = vcmax25 / marea,
         jmax.mass = jmax25 / marea,
         n.trt = ifelse(treatment == "AS" | treatment == "NO3", "n.added",
                        "no.n"),
         s.trt = ifelse(treatment == "AS" | treatment == "S", "s.added",
                        "no.s"),
         treatment = factor(treatment, levels = c("C", "NO3", "AS", "S")))

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
leaf.n.no3 <- lmer(leaf.n ~ soil.no3n.norm + mineral.pH + (1 | site), 
               data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(leaf.n.no3)
qqnorm(residuals(leaf.n.no3))
qqline(residuals(leaf.n.no3))
hist(residuals(leaf.n.no3))
shapiro.test(residuals(leaf.n.no3))
outlierTest(leaf.n.no3)

# Model output
summary(leaf.n.no3)
Anova(leaf.n.no3)
r.squaredGLMM(leaf.n.no3)

# Post-hoc tests
test(emtrends(leaf.n.no3, ~1, var = "soil.no3n.norm"))

# Plot prep
nmass.no3n.trend <- data.frame(emmeans(leaf.n.no3, ~1, "soil.no3n.norm", 
                                       at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
nmass.no3n.plot <- ggplot(data = plot_data, aes(x = soil.no3n.norm, y = leaf.n/100)) +
  geom_point(data = plot_data, aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = nmass.no3n.trend, 
              aes(y = emmean/100,
                  ymin = lower.CL/100, ymax = upper.CL/100),
              alpha = 0.3) +
  geom_smooth(data = nmass.no3n.trend, aes(y = emmean/100), size = 2, se = FALSE,
              color = "black", method = 'lm', linetype = "dashed") +
  scale_x_continuous(limits = c(0, 25), 
                     breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0.015, 0.035), 
                     breaks = seq(0.015, 0.035, 0.005)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"["mass"]*" (gN g"["dry_mass"]*""^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
nmass.no3n.plot

##########################################################################
## Leaf mass per area - soil NO3-N
##########################################################################
marea.no3 <- lmer(marea ~ soil.no3n.norm + mineral.pH + (1 | site),
              data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(marea.no3)
qqnorm(residuals(marea.no3))
qqline(residuals(marea.no3))
hist(residuals(marea.no3))
shapiro.test(residuals(marea.no3))
outlierTest(marea.no3)

# Model output
summary(marea.no3)
Anova(marea.no3)
r.squaredGLMM(marea.no3)

# Plot
marea.no3n.plot <- ggplot(data = plot_data, aes(x = soil.no3n.norm, y = marea)) +
  geom_point(aes(shape = treatment), size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 25), 
                     breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("M"[area]*" (g"["dry_mass"]*" m"["leaf"]*""^"-2"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
marea.no3n.plot

##########################################################################
## Narea - soil NO3-N
##########################################################################
narea.no3 <- lmer(narea ~ soil.no3n.norm + mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea.no3)
qqnorm(residuals(narea.no3))
qqline(residuals(narea.no3))
hist(residuals(narea.no3))
shapiro.test(residuals(narea.no3))
outlierTest(narea.no3)

# Model output
summary(narea.no3)
Anova(narea.no3)
r.squaredGLMM(narea.no3)

# Pairwise comparisons
test(emtrends(narea.no3, ~1, var = "soil.no3n.norm"))

# Plot prep
narea.no3n.trend <- data.frame(emmeans(narea.no3, ~1, "soil.no3n.norm", 
                                       at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
narea.no3n.plot <- ggplot(data = plot_data, aes(x = soil.no3n.norm, y = narea)) +
  geom_point(data = plot_data, aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = narea.no3n.trend, 
              aes(y = emmean,
                  ymin = lower.CL, ymax = upper.CL),
              alpha = 0.3) +
  geom_smooth(data = narea.no3n.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black", method = 'lm') +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 2.4), breaks = seq(0, 2.4, 0.6)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("N"["area"]*" (gN m"["leaf"]*""^"-2"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
narea.no3n.plot

##########################################################################
## Anet - soil NO3-N
##########################################################################
data$a400[data$a400 < 0.2] <- NA

a400.no3 <- lmer(sqrt(a400) ~ soil.no3n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(a400.no3)
qqnorm(residuals(a400.no3))
qqline(residuals(a400.no3))
hist(residuals(a400.no3))
shapiro.test(residuals(a400.no3))
outlierTest(a400.no3)

# Model output
summary(a400.no3)
Anova(a400.no3)
r.squaredGLMM(a400.no3)

# Plot
a400.no3n.plot <- ggplot(data = plot_data, aes(x = soil.no3n.norm, y = a400)) +
  geom_jitter(aes(shape = treatment), width = 1, size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("A"["net"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
a400.no3n.plot

##########################################################################
## Vcmax25 - soil NO3-N
##########################################################################
vcmax.no3 <- lmer(vcmax25 ~ soil.no3n.norm + mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(vcmax.no3)
qqnorm(residuals(vcmax.no3))
qqline(residuals(vcmax.no3))
hist(residuals(vcmax.no3))
shapiro.test(residuals(vcmax.no3))
outlierTest(vcmax.no3)

# Model output
summary(vcmax.no3)
Anova(vcmax.no3)
r.squaredGLMM(vcmax.no3)

# Pairwise comparisons
test(emtrends(vcmax.no3, ~1, var = "soil.no3n.norm"))

# Plot prep
vcmax.no3n.trend <- data.frame(emmeans(vcmax.no3, ~1, "soil.no3n.norm", 
                                       at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
vcmax.no3n.plot <- ggplot(data = plot_data, aes(x = soil.no3n.norm, y = vcmax25)) +
  geom_point(data = plot_data, aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = vcmax.no3n.trend, 
              aes(y = emmean,
                  ymin = lower.CL, ymax = upper.CL),
              alpha = 0.3) +
  geom_smooth(data = vcmax.no3n.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black", method = 'lm', linetype = "dashed") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("V"["cmax25"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
vcmax.no3n.plot



##########################################################################
## Jmax25 - soil NO3-N
##########################################################################
jmax.no3 <- lmer(jmax25 ~ soil.no3n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(jmax.no3)
qqnorm(residuals(jmax.no3))
qqline(residuals(jmax.no3))
hist(residuals(jmax.no3))
shapiro.test(residuals(jmax.no3))
outlierTest(jmax.no3)

# Model output
summary(jmax.no3)
Anova(jmax.no3)
r.squaredGLMM(jmax.no3)

## Pairwise comparisons
test(emtrends(jmax.no3, ~1, var = "soil.no3n.norm"))

# Plot prep
jmax.no3n.trend <- data.frame(emmeans(jmax.no3, ~1, "soil.no3n.norm", 
                                      at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
jmax.no3n.plot <- ggplot(data = plot_data, aes(x = soil.no3n.norm, y = jmax25)) +
  geom_point(data = plot_data, aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = jmax.no3n.trend, 
              aes(y = emmean,
                  ymin = lower.CL, ymax = upper.CL),
              alpha = 0.3) +
  geom_smooth(data = jmax.no3n.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black", method = 'lm') +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("J"["max25"]*" (μmol m"^"-2"*" s"^"-1"*")")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
jmax.no3n.plot


##########################################################################
## Jmax25:Vcmax25 - soil NO3-N
##########################################################################
vjmax.no3 <- lmer(jmax.vcmax ~ soil.no3n.norm + mineral.pH + (1 | site),
              data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vjmax.no3)
qqnorm(residuals(vjmax.no3))
qqline(residuals(vjmax.no3))
hist(residuals(vjmax.no3))
shapiro.test(residuals(vjmax.no3))
outlierTest(vjmax.no3)

# Model output
summary(vjmax.no3)
Anova(vjmax.no3)
r.squaredGLMM(vjmax.no3)

# Plot
jvmax.no3n.plot <- ggplot(data = plot_data, aes(x = soil.no3n.norm, y = jmax.vcmax)) +
  geom_point(data = plot_data, aes(shape = treatment), size = 4, alpha = 0.75) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(1, 3), breaks = seq(1, 3, 0.5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold("J"["max25"]*": V"["cmax25"]*" (unitless)")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
jvmax.no3n.plot

##########################################################################
## chi - soil NO3-N
##########################################################################
chi.no3 <- lmer(chi ~ soil.no3n.norm + mineral.pH + (1 | site),
            data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(chi.no3)
qqnorm(residuals(chi.no3))
qqline(residuals(chi.no3))
hist(residuals(chi.no3))
shapiro.test(residuals(chi.no3))
outlierTest(chi.no3)

# Model output
summary(chi.no3)
Anova(chi.no3)
r.squaredGLMM(chi.no3)

## Pairwise comparisons
test(emtrends(chi.no3,  ~1, var = "soil.no3n.norm"))

# Plot prep
chi.no3.trend <- data.frame(emmeans(chi.no3, ~1, "soil.no3n.norm", 
                                    at = list(soil.no3n.norm = seq(0, 25, 0.1))))

# Plot
jmax.no3n.plot <- ggplot(data = plot_data, aes(x = soil.no3n.norm, y = chi)) +
  geom_point(data = plot_data, aes(shape = treatment), size = 4, alpha = 0.75) +
  geom_ribbon(data = chi.no3.trend, 
              aes(y = emmean,
                  ymin = lower.CL, ymax = upper.CL),
              alpha = 0.3) +
  geom_smooth(data = chi.no3.trend, aes(y = emmean), size = 2, se = FALSE,
              color = "black", method = 'lm') +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("C" = "no N; no S",
                                "NO3" = "+ N; no S",
                                "AS" = "+ N; + S",
                                "S" = "no N; + S")) +
  labs(x = expression(bold("Soil NO"["3"]*"-N (μg NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")")),
       y = expression(bold(chi*" (unitless)")),
       shape = "Treatment") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  pubtheme
jmax.no3n.plot

##########################################################################
## PNUE - soil NO3-N
##########################################################################
data$pnue[data$pnue < 0] <- NA

pnue.no3 <- lmer(sqrt(pnue) ~ soil.no3n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue.no3)
qqnorm(residuals(pnue.no3))
qqline(residuals(pnue.no3))
hist(residuals(pnue.no3))
shapiro.test(residuals(pnue.no3))
outlierTest(pnue)

# Model output
summary(pnue.no3)
Anova(pnue.no3)
r.squaredGLMM(pnue.no3)

# Pairwise comparisons
test(emtrends(pnue.no3, ~1, var = "soil.no3n.norm"))

##########################################################################
## Narea.chi - soil NO3-N
##########################################################################
narea.chi.no3 <- lmer(narea.chi ~ soil.no3n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea.chi.no3)
qqnorm(residuals(narea.chi.no3))
qqline(residuals(narea.chi.no3))
hist(residuals(narea.chi.no3))
shapiro.test(residuals(narea.chi.no3))
outlierTest(narea.chi.no3)

# Model output
summary(narea.chi.no3)
Anova(narea.chi.no3)
r.squaredGLMM(narea.chi.no3)

# Post-hoc tests
test(emtrends(narea.chi.no3, ~1, var = "soil.no3n.norm"))

##########################################################################
## Vcmax25.chi - soil NO3-N
##########################################################################
data$vcmax.chi[85] <- NA

vcmax.chi.no3 <- lmer(vcmax.chi ~ soil.no3n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vcmax.chi.no3)
qqnorm(residuals(vcmax.chi.no3))
qqline(residuals(vcmax.chi.no3))
hist(residuals(vcmax.chi.no3))
shapiro.test(residuals(vcmax.chi.no3))
outlierTest(vcmax.chi.no3)

# Model output
summary(vcmax.chi.no3)
Anova(vcmax.chi.no3)
r.squaredGLMM(vcmax.chi.no3)

# Post-hoc tests
test(emtrends(vcmax.chi, ~1, var = "soil.no3n.norm"))


##########################################################################
##########################################################################
## Ammonium-specific analyses
##########################################################################
##########################################################################


##########################################################################
## Nleaf
##########################################################################
leaf.n.nh4 <- lmer(leaf.n ~ soil.nh4n.norm + mineral.pH + (1 | site), 
               data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(leaf.n.nh4)
qqnorm(residuals(leaf.n.nh4))
qqline(residuals(leaf.n.nh4))
hist(residuals(leaf.n.nh4))
shapiro.test(residuals(leaf.n.nh4))
outlierTest(leaf.n.nh4)

# Model output
summary(leaf.n.nh4)
Anova(leaf.n.nh4)
r.squaredGLMM(leaf.n.nh4)

# Post-hoc tests
test(emtrends(leaf.n.nh4, ~1, var = "soil.nh4n.norm"))

##########################################################################
## Leaf mass per area
##########################################################################
marea.nh4 <- lmer(marea ~ soil.nh4n.norm + mineral.pH + (1 | site),
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

##########################################################################
## Vcmax25
##########################################################################
vcmax.nh4 <- lmer(vcmax25 ~ soil.nh4n.norm + mineral.pH + (1 | site), 
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


##########################################################################
## PNUE
##########################################################################
data$pnue[data$pnue < 0] <- NA

pnue.nh4 <- lmer(sqrt(pnue) ~ soil.nh4n.norm + mineral.pH + (1 | site),
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

##########################################################################
## Vcmax25.chi
##########################################################################
data$vcmax.chi[85] <- NA

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