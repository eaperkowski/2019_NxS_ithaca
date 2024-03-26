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
                        "no.s"))

##########################################################################
## Nleaf - soil N
##########################################################################
leaf.n <- lmer(leaf.n ~ soil.n.norm + mineral.pH + (1 | site), 
               data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(leaf.n)
qqnorm(residuals(leaf.n))
qqline(residuals(leaf.n))
hist(residuals(leaf.n))
shapiro.test(residuals(leaf.n))
outlierTest(leaf.n)

# Model output
summary(leaf.n)
Anova(leaf.n)
r.squaredGLMM(leaf.n)

# Post-hoc tests
test(emtrends(leaf.n, ~1, var = "soil.n.norm"))

##########################################################################
## Leaf mass per area - soil N
##########################################################################
marea <- lmer(marea ~ soil.n.norm + mineral.pH + (1 | site),
              data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(marea)
qqnorm(residuals(marea))
qqline(residuals(marea))
hist(residuals(marea))
shapiro.test(residuals(marea))
outlierTest(marea)

# Model output
summary(marea)
Anova(marea)
r.squaredGLMM(marea)

# Post-hoc tests
test(emtrends(marea, ~1, var = "soil.n.norm"))

##########################################################################
## Narea - soil N
##########################################################################
narea <- lmer(narea ~ soil.n.norm + mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea)
qqnorm(residuals(narea))
qqline(residuals(narea))
hist(residuals(narea))
shapiro.test(residuals(narea))
outlierTest(narea)

# Model output
summary(narea)
Anova(narea)
r.squaredGLMM(narea)

# Pairwise comparisons
test(emtrends(narea, ~1, var = "soil.n.norm"))

##########################################################################
## Anet - soil N
##########################################################################
data$a400[data$a400 < 0.2] <- NA

a400 <- lmer(sqrt(a400) ~ soil.n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(a400)
qqnorm(residuals(a400))
qqline(residuals(a400))
hist(residuals(a400))
shapiro.test(residuals(a400))
outlierTest(a400)

# Model output
summary(a400)
Anova(a400)
r.squaredGLMM(a400)

##########################################################################
## Anet - Narea
##########################################################################
anet.narea <- lmer(sqrt(a400) ~ narea +  (1 | site), 
                   data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(anet.narea)
qqnorm(residuals(anet.narea))
qqline(residuals(anet.narea))
hist(residuals(anet.narea))
shapiro.test(residuals(anet.narea))
outlierTest(anet.narea)

# Model output
summary(anet.narea)
Anova(anet.narea)
r.squaredGLMM(anet.narea)

##########################################################################
## Vcmax25 - soil N
##########################################################################
vcmax <- lmer(vcmax25 ~ soil.n.norm + mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(vcmax)
qqnorm(residuals(vcmax))
qqline(residuals(vcmax))
hist(residuals(vcmax))
shapiro.test(residuals(vcmax))
outlierTest(vcmax)

# Model output
summary(vcmax)
Anova(vcmax)
r.squaredGLMM(vcmax)

# Pairwise comparisons
test(emtrends(vcmax, ~1, var = "soil.n.norm"))

##########################################################################
## Vcmax25 - leaf N
##########################################################################
vcmax.nleaf <- lmer(vcmax25 ~ narea + (1 | site), 
                    data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(vcmax.nleaf)
qqnorm(residuals(vcmax.nleaf))
qqline(residuals(vcmax.nleaf))
hist(residuals(vcmax.nleaf))
shapiro.test(residuals(vcmax.nleaf))
outlierTest(vcmax.nleaf)

# Model output
summary(vcmax.nleaf)
Anova(vcmax.nleaf)
r.squaredGLMM(vcmax.nleaf)

# Test Narea-Vcmax25 slope
test(emtrends(vcmax.nleaf, ~1, "narea"))

##########################################################################
## Jmax25 - soil N
##########################################################################
jmax <- lmer(jmax25 ~ soil.n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(jmax)
qqnorm(residuals(jmax))
qqline(residuals(jmax))
hist(residuals(jmax))
shapiro.test(residuals(jmax))
outlierTest(jmax)

# Model output
summary(jmax)
Anova(jmax)
r.squaredGLMM(jmax)

## Pairwise comparisons
test(emtrends(jmax, ~1, var = "soil.n.norm"))

##########################################################################
## Jmax25 - leaf N
##########################################################################
jmax.nleaf <- lmer(jmax25 ~ narea + (1 | site), 
                   data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(jmax.nleaf)
qqnorm(residuals(jmax.nleaf))
qqline(residuals(jmax.nleaf))
hist(residuals(jmax.nleaf))
shapiro.test(residuals(jmax.nleaf))
outlierTest(jmax.nleaf)

# Model output
summary(jmax.nleaf)
Anova(jmax.nleaf)
r.squaredGLMM(jmax.nleaf)

# Test Narea-Jmax25 slope
test(emtrends(jmax.nleaf, ~1, "narea"))

##########################################################################
## Jmax25:Vcmax25 - soil N
##########################################################################
vjmax <- lmer(jmax.vcmax ~ soil.n.norm + mineral.pH + (1 | site),
              data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vjmax)
qqnorm(residuals(vjmax))
qqline(residuals(vjmax))
hist(residuals(vjmax))
shapiro.test(residuals(vjmax))
outlierTest(vjmax)

# Model output
summary(vjmax)
Anova(vjmax)
r.squaredGLMM(vjmax)

##########################################################################
## Jmax25:Vcmax25 - leaf N
##########################################################################
vjmax.nleaf <- lmer(jmax.vcmax ~ narea + (1 | site), 
                    data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(vjmax.nleaf)
qqnorm(residuals(vjmax.nleaf))
qqline(residuals(vjmax.nleaf))
hist(residuals(vjmax.nleaf))
shapiro.test(residuals(vjmax.nleaf))
outlierTest(vjmax.nleaf)

# Model output
summary(vjmax.nleaf)
Anova(vjmax.nleaf)
r.squaredGLMM(vjmax.nleaf)

##########################################################################
## chi - soil N
##########################################################################
chi <- lmer(chi ~ soil.n.norm + mineral.pH + (1 | site),
            data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(chi)
qqnorm(residuals(chi))
qqline(residuals(chi))
hist(residuals(chi))
shapiro.test(residuals(chi))
outlierTest(chi)

# Model output
summary(chi)
Anova(chi)
r.squaredGLMM(chi)

## Pairwise comparisons
test(emtrends(chi,  ~1, var = "soil.n.norm"))

##########################################################################
## chi - leaf N
##########################################################################
data$chi[66] <- NA

chi.nleaf <- lmer(chi ~ narea + (1 | site), 
                  data = subset(data, nrcs.code == "ACSA3"))

# Check model assumptions
plot(chi.nleaf)
qqnorm(residuals(chi.nleaf))
qqline(residuals(chi.nleaf))
hist(residuals(chi.nleaf))
shapiro.test(residuals(chi.nleaf))
outlierTest(chi.nleaf)

# Model output
summary(chi.nleaf)
Anova(chi.nleaf)
r.squaredGLMM(chi.nleaf)

# Test Narea-Anet slope
test(emtrends(chi.nleaf, ~1, "narea"))

##########################################################################
## PNUE - soil N
##########################################################################
data$pnue[data$pnue < 0] <- NA

pnue <- lmer(sqrt(pnue) ~ soil.n.norm + mineral.pH + (1 | site),
             data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(pnue)
qqnorm(residuals(pnue))
qqline(residuals(pnue))
hist(residuals(pnue))
shapiro.test(residuals(pnue))
outlierTest(pnue)

# Model output
summary(pnue)
Anova(pnue)
r.squaredGLMM(pnue)

# Pairwise comparisons
test(emtrends(pnue, ~1, var = "soil.n.norm"))

##########################################################################
## Narea.chi - soil N
##########################################################################
narea.chi <- lmer(narea.chi ~ soil.n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(narea.chi)
qqnorm(residuals(narea.chi))
qqline(residuals(narea.chi))
hist(residuals(narea.chi))
shapiro.test(residuals(narea.chi))
outlierTest(narea.chi)

# Model output
summary(narea.chi)
Anova(narea.chi)
r.squaredGLMM(narea.chi)

# Post-hoc tests
test(emtrends(narea.chi, ~1, var = "soil.n.norm"))

##########################################################################
## Vcmax25.chi - soil N
##########################################################################
data$vcmax.chi[85] <- NA

vcmax.chi <- lmer(vcmax.chi ~ soil.n.norm + mineral.pH + (1 | site),
                  data = subset(data, nrcs.code == "ACSA3"))

# Check normality assumptions
plot(vcmax.chi)
qqnorm(residuals(vcmax.chi))
qqline(residuals(vcmax.chi))
hist(residuals(vcmax.chi))
shapiro.test(residuals(vcmax.chi))
outlierTest(vcmax.chi)

# Model output
summary(vcmax.chi)
Anova(vcmax.chi)
r.squaredGLMM(vcmax.chi)

# Post-hoc tests
test(emtrends(vcmax.chi, ~1, var = "soil.n.norm"))
    