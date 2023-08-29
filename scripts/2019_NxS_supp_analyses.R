##########################################################################
## Load libraries and import data
##########################################################################
# Libraries
library(lme4)
library(emmeans)
library(car)
library(tidyverse)
library(MuMIn)
library(multcomp)
library(multcompView)

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
# Leaf temp. effect on Anet - linear regression
##########################################################################
data$a400[data$a400 < 0.2] <- NA

a400 <- lmer(sqrt(a400) ~ leaf.temp + (1 | nrcs.code) + (1 | site), 
             data = subset(data, nrcs.code == "ACRU" |
                             nrcs.code == "ACSA" |
                             nrcs.code == "QURU" |
                             nrcs.code == "FAGR" |
                             nrcs.code == "FRAM"))

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

##########################################################################
# Leaf temp. effect on Anet - log polynomial
##########################################################################
a400.nls <- nls(formula = log(a400) ~ a + b*leaf.temp + c*(leaf.temp^2),
                start = list(a = 9.422, b = -0.572, c = 0.0099),
                data = subset(data, nrcs.code == "ACRU" |
                                nrcs.code == "ACSA" |
                                nrcs.code == "QURU" |
                                nrcs.code == "FAGR" |
                                nrcs.code == "FRAM"))
plot(residuals(a400.nls))
coef(a400.nls)

##########################################################################
# Leaf temp. effect on gsw - linear regression
##########################################################################
gs400 <- lmer(log(gsw) ~ leaf.temp + (1 | nrcs.code) +
                (1 | site), 
              data = subset(data, nrcs.code == "ACRU" |
                              nrcs.code == "ACSA" |
                              nrcs.code == "QURU" |
                              nrcs.code == "FAGR" |
                              nrcs.code == "FRAM"))

# Check model assumptions
plot(gs400)
qqnorm(residuals(gs400))
qqline(residuals(gs400))
hist(residuals(gs400))
shapiro.test(residuals(gs400))
outlierTest(gs400)

# Model output
summary(gs400)
Anova(gs400)

##########################################################################
# Leaf temp. effect on gsw - log polynomial
##########################################################################
gs400.nls <- nls(formula = log(gsw) ~ a + b*leaf.temp + c*(leaf.temp^2),
                 start = list(a = -0.17, b = -0.18, c = 0.0027),
                 data = data)
coef(gs400.nls)

##########################################################################
##########################################################################
## Isolating effects of soil pH on measured leaf traits in A. saccharum
##########################################################################
##########################################################################

##########################################################################
## Nmass - soil N for only C and S plots
##########################################################################
leaf.n <- lmer(leaf.n ~ mineral.pH + (1 | site), 
               data = subset(data, nrcs.code == "ACSA3" &
                               (treatment == "C" | treatment == "S")))

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

##########################################################################
## Leaf mass per area - soil N
##########################################################################
marea <- lmer(marea ~ mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3" &
                              (treatment == "C" | treatment == "S")))

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
test(emtrends(marea, ~1, var = "mineral.pH"))

##########################################################################
## Narea - soil N
##########################################################################
narea <- lmer(narea ~ mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3" &
                              (treatment == "C" | treatment == "S")))

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
test(emtrends(narea, ~1, var = "mineral.pH"))

##########################################################################
## Anet,area - soil N
##########################################################################
data$a400[data$a400 < 0.2] <- NA

a400 <- lmer(log(a400) ~ mineral.pH +  (1 | site), 
             data = subset(data, nrcs.code == "ACSA3" &
                             (treatment == "C" | treatment == "S")))

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
## Vcmax25 - soil N
##########################################################################
vcmax <- lmer(vcmax25 ~ mineral.pH + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3" &
                              (treatment == "C" | treatment == "S")))

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
emtrends(vcmax, ~1,  var = "mineral.pH")

##########################################################################
## Jmax25 - soil N
##########################################################################
jmax <- lmer(jmax25 ~ mineral.pH + (1 | site), 
             data = subset(data, nrcs.code == "ACSA3" &
                             (treatment == "C" | treatment == "S")))

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

##########################################################################
## Jmax25:Vcmax25 - soil N
##########################################################################
vjmax <- lmer(log(jmax.vcmax) ~ mineral.pH + (1 | site),
              data = subset(data, nrcs.code == "ACSA3" &
                              (treatment == "C" | treatment == "S")))

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
## chi - soil N
##########################################################################
chi <- lmer(chi ~ mineral.pH + (1 | site), 
            data = subset(data, nrcs.code == "ACSA3" &
                            (treatment == "C" | treatment == "S")))

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

##########################################################################
## PNUE - soil N
##########################################################################
data$pnue[data$pnue < 0] <- NA

pnue <- lmer(log(pnue) ~ mineral.pH + (1 | site), 
             data = subset(data, nrcs.code == "ACSA3" &
                             (treatment == "C" | treatment == "S")))

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

##########################################################################
## Narea.chi - soil N
##########################################################################
narea.chi <- lmer(narea.chi ~ mineral.pH + (1 | site), 
                  data = subset(data, nrcs.code == "ACSA3" &
                                  (treatment == "C" | treatment == "S")))

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
test(emtrends(narea.chi, ~1, var = "mineral.pH"))

##########################################################################
## Vcmax25.chi - soil N
##########################################################################
vcmax.chi <- lmer(vcmax.chi ~ mineral.pH + (1 | site), 
                  data = subset(data, nrcs.code == "ACSA3" &
                                  (treatment == "C" | treatment == "S")))

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
test(emtrends(vcmax.chi, ~1, var = "mineral.pH"))

##########################################################################
##########################################################################
## Effects of N availability and soil pH on traits across ALL species
##########################################################################
##########################################################################

##########################################################################
## Nleaf - soil N
##########################################################################
leaf.n <- lmer(leaf.n ~ soil.n.norm + mineral.pH + nrcs.code + 
                 (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                             nrcs.code == "ACSA3" |
                                             nrcs.code == "QURU" |
                                             nrcs.code == "FAGR" | 
                                             nrcs.code == "FRAM2"))

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
emmeans(leaf.n, pairwise~nrcs.code)

##########################################################################
## Leaf mass per area - soil N
##########################################################################
marea <- lmer(marea ~ soil.n.norm + mineral.pH + nrcs.code + 
                (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                            nrcs.code == "ACSA3" |
                                            nrcs.code == "QURU" |
                                            nrcs.code == "FAGR" | 
                                            nrcs.code == "FRAM2"))

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

# Pairwise comparisons
test(emtrends(marea, ~1, var = "soil.n.norm"))
emmeans(marea, pairwise~nrcs.code)

##########################################################################
## Narea - soil N
##########################################################################
data$narea[12] <- NA
narea <- lmer(narea ~ soil.n.norm + mineral.pH + nrcs.code +
                (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                            nrcs.code == "ACSA3" |
                                            nrcs.code == "QURU" |
                                            nrcs.code == "FAGR" | 
                                            nrcs.code == "FRAM2"))

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
emmeans(narea, pairwise~nrcs.code)


##########################################################################
## Net photosynthesis - soil N
##########################################################################
data$a400[data$a400 < 0.2] <- NA
a400 <- lmer(sqrt(a400) ~ soil.n.norm + mineral.pH + nrcs.code + 
               (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                           nrcs.code == "ACSA3" |
                                           nrcs.code == "QURU" |
                                           nrcs.code == "FAGR" | 
                                           nrcs.code == "FRAM2"))

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

# Pairwise comparisons
test(emtrends(a400, ~1, var = "soil.n.norm"))
emmeans(a400, pairwise~nrcs.code)

##########################################################################
## Net photosynthesis-leaf N
##########################################################################
data$a400[data$a400 < 0.2] <- NA

a400.nleaf <- lmer(sqrt(a400) ~ narea + (1 | nrcs.code) + (1 | site), 
                   data = subset(data, nrcs.code == "ACRU" |
                                   nrcs.code == "ACSA3" |
                                   nrcs.code == "QURU" |
                                   nrcs.code == "FAGR" | 
                                   nrcs.code == "FRAM2"))

# Check model assumptions
plot(a400.nleaf)
qqnorm(residuals(a400.nleaf))
qqline(residuals(a400.nleaf))
hist(residuals(a400.nleaf))
shapiro.test(residuals(a400.nleaf))
outlierTest(a400.nleaf)

# Model output
summary(a400.nleaf)
Anova(a400.nleaf)
r.squaredGLMM(a400.nleaf)

# Test Narea-Anet slope
test(emtrends(a400.nleaf, ~1, "narea"))

##########################################################################
## Vcmax25 - soil N
##########################################################################
data$vcmax25[11] <- NA
vcmax <- lmer(vcmax25 ~ soil.n.norm + mineral.pH + nrcs.code  +
                (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                            nrcs.code == "ACSA3" |
                                            nrcs.code == "QURU" |
                                            nrcs.code == "FAGR" | 
                                            nrcs.code == "FRAM2"))

# Check normality assumptions
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
emmeans(vcmax, pairwise~nrcs.code)

##########################################################################
## Vcmax25 - leaf N
##########################################################################
vcmax.nleaf <- lmer(vcmax25 ~ narea + (1 | nrcs.code) +  (1 | site), 
                    data = subset(data, nrcs.code == "ACRU" |
                                    nrcs.code == "ACSA3" |
                                    nrcs.code == "QURU" |
                                    nrcs.code == "FAGR" | 
                                    nrcs.code == "FRAM2"))

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

# Test Narea-Anet slope
test(emtrends(vcmax.nleaf, ~1, "narea"))

##########################################################################
## Jmax25 - soil N
##########################################################################
data$jmax25[11] <- NA

jmax <- lmer(jmax25 ~ soil.n.norm + mineral.pH + nrcs.code + 
               (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                           nrcs.code == "ACSA3" |
                                           nrcs.code == "QURU" |
                                           nrcs.code == "FAGR" | 
                                           nrcs.code == "FRAM2"))

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
emmeans(jmax, pairwise~nrcs.code)

##########################################################################
## Jmax25 - leaf N
##########################################################################
jmax.nleaf <- lmer(jmax25 ~ narea + (1 | nrcs.code) + (1 | site), 
                   data = subset(data, nrcs.code == "ACRU" |
                                   nrcs.code == "ACSA3" |
                                   nrcs.code == "QURU" |
                                   nrcs.code == "FAGR" | 
                                   nrcs.code == "FRAM2"))

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

# Test Narea-Anet slope
test(emtrends(jmax.nleaf, ~1, "narea"))

##########################################################################
## Jmax25:Vcmax25 - soil N
##########################################################################
vjmax <- lmer(log(jmax.vcmax) ~ soil.n.norm + mineral.pH + nrcs.code + 
                (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                            nrcs.code == "ACSA3" |
                                            nrcs.code == "QURU" |
                                            nrcs.code == "FAGR" | 
                                            nrcs.code == "FRAM2"))

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

## Pairwise comparisons
test(emtrends(vjmax, ~1, var = "soil.n.norm"))
emmeans(vjmax, pairwise~nrcs.code)

##########################################################################
## Jmax25:Vcmax25 - leaf N
##########################################################################
vjmax.nleaf <- lmer(log(jmax.vcmax) ~ narea + (1 | nrcs.code) + (1 | site), 
                    data = subset(data, nrcs.code == "ACRU" |
                                    nrcs.code == "ACSA3" |
                                    nrcs.code == "QURU" |
                                    nrcs.code == "FAGR" | 
                                    nrcs.code == "FRAM2"))

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

# Test Narea-Anet slope
test(emtrends(vjmax.nleaf, ~1, "narea"))

##########################################################################
## chi - soil N
##########################################################################
chi <- lmer(chi ~ soil.n.norm + mineral.pH + nrcs.code + 
              (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                          nrcs.code == "ACSA3" |
                                          nrcs.code == "QURU" |
                                          nrcs.code == "FAGR" |
                                          nrcs.code == "FRAM2"))

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
emmeans(chi, pairwise~nrcs.code)

##########################################################################
## chi - leaf N
##########################################################################
data$chi[66] <- NA

chi.nleaf <- lmer(chi ~ narea + (1 | nrcs.code) + (1 | site), 
                  data = subset(data, nrcs.code == "ACRU" |
                                  nrcs.code == "ACSA3" |
                                  nrcs.code == "QURU" |
                                  nrcs.code == "FAGR" | 
                                  nrcs.code == "FRAM2"))

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

pnue <- lmer(sqrt(pnue) ~ soil.n.norm + mineral.pH + nrcs.code +
               (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                           nrcs.code == "ACSA3" |
                                           nrcs.code == "QURU" |
                                           nrcs.code == "FAGR" |
                                           nrcs.code == "FRAM2"))

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
emmeans(pnue, pairwise~nrcs.code)

##########################################################################
## Narea.chi - soil N
##########################################################################
data$narea.chi[12] <- NA

narea.chi <- lmer(narea.chi ~ soil.n.norm + mineral.pH + nrcs.code +
                    (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                                nrcs.code == "ACSA3" |
                                                nrcs.code == "QURU" |
                                                nrcs.code == "FAGR" |
                                                nrcs.code == "FRAM2"))

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
emmeans(narea.chi, pairwise~nrcs.code)

##########################################################################
## Vcmax25.chi - soil N
##########################################################################
data$vcmax.chi[c(11, 37)] <- NA

vcmax.chi <- lmer(vcmax.chi ~ soil.n.norm + mineral.pH + nrcs.code +
                    (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                                nrcs.code == "ACSA3" |
                                                nrcs.code == "QURU" |
                                                nrcs.code == "FAGR" |
                                                nrcs.code == "FRAM2"))

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
cld(emmeans(vcmax.chi, pairwise~nrcs.code))

##########################################################################
## Vcmax25.chi - leaf N
##########################################################################
vcmax.chi.nleaf <- lmer(vcmax.chi ~ narea + (1 | nrcs.code) + (1 | site), 
                        data = subset(data, nrcs.code == "ACRU" |
                                        nrcs.code == "ACSA3" |
                                        nrcs.code == "QURU" |
                                        nrcs.code == "FAGR" | 
                                        nrcs.code == "FRAM2"))

# Check model assumptions
plot(vcmax.chi.nleaf)
qqnorm(residuals(vcmax.chi.nleaf))
qqline(residuals(vcmax.chi.nleaf))
hist(residuals(vcmax.chi.nleaf))
shapiro.test(residuals(vcmax.chi.nleaf))
outlierTest(vcmax.chi.nleaf)

# Model output
summary(vcmax.chi.nleaf)
Anova(vcmax.chi.nleaf)
r.squaredGLMM(vcmax.chi.nleaf)

# Test Narea-Anet slope
test(emtrends(vcmax.chi.nleaf, ~1, "narea"))
