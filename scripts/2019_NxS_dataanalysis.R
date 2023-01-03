##########################################################################
## Load libraries and import data
##########################################################################
# Libraries
library(lme4)
library(emmeans)
library(car)
library(tidyverse)
library(dplyr)
library(MuMIn)
library(multcomp)
library(multcompView)

emm_options(opt.digits = FALSE)

# Import datasheet
data <- read.csv("../data/2019_NxS_datasheet.csv",
                 stringsAsFactors = FALSE,
                 na.strings = "NA")

n.spp <- data %>%
  # filter(nrcs.code == "ACRU" | nrcs.code == "ACSA3" |
  #          nrcs.code == "QURU" | nrcs.code == "FAGR" | 
  #          nrcs.code == "FRAM2") %>%
  # group_by(site, treatment, nrcs.code) %>%
  group_by(nrcs.code) %>%
  summarize(n = n())

## Mean and SE of soil nitrogen availabilityfor discussion
mean(unique(data$soil.n.norm))
sd(unique(data$soil.n.norm))

## Mean and SE of soil pH for discussion
mean(unique(data$mineral.pH))
sd(unique(data$mineral.pH))

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
emmeans(leaf.n, ~1, "soil.n.norm", at = list(soil.n.norm = 0))
cld(emmeans(leaf.n, pairwise~nrcs.code))


# Emmean output for fig making
leaf.n.pairwise <- data.frame(variable = "leaf.n",
                              cld(emmeans(leaf.n, 
                                          ~nrcs.code),
                                  Letters = LETTERS))

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
cld(emmeans(marea, pairwise~nrcs.code))

## Figure slopes, intercepts, percent change values
emtrends(marea, ~1,var = "soil.n.norm")
emmeans(marea,~1, "soil.n.norm", at = list(soil.n.norm = 0))
emmeans(marea, pairwise~nrcs.code)

# Emmean output for fig making
marea.pairwise <- data.frame(variable = "marea",
                             cld(emmeans(marea, ~nrcs.code),
                                 Letters = LETTERS))

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
emmeans(narea, ~1, at = list(soil.n.norm = 0))
emmeans(narea, pairwise~nrcs.code)

# Emmean output for fig making
narea.pairwise <- data.frame(variable = "narea",
                             cld(emmeans(narea, ~nrcs.code),
                                 Letters = LETTERS))

##########################################################################
## N.rubisco - soil N
##########################################################################
p.rubisco <- lmer(p.rubisco ~ soil.n.norm + mineral.pH + nrcs.code +
                    (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                                nrcs.code == "ACSA3" |
                                                nrcs.code == "QURU" |
                                                nrcs.code == "FAGR" |
                                                nrcs.code == "FRAM2"))

# Check normality assumptions
plot(p.rubisco)
qqnorm(residuals(p.rubisco))
qqline(residuals(p.rubisco))
hist(residuals(p.rubisco))
shapiro.test(residuals(p.rubisco))
outlierTest(p.rubisco)

# Model output
summary(p.rubisco)
Anova(p.rubisco)
r.squaredGLMM(p.rubisco)

# Species pairwise comparison
test(emtrends(p.rubisco, ~1, "soil.n.norm"))
emmeans(p.rubisco, pairwise~nrcs.code)

# Emmean output for fig making
nrub.pairwise <- data.frame(variable = "p.rubisco", 
                            cld(emmeans(p.rubisco, ~nrcs.code),
                                Letters = LETTERS))

##########################################################################
## N.bioenergetics - soil N
##########################################################################
data$p.bioe[11] <- NA

p.bioe <- lmer(p.bioe ~ soil.n.norm + mineral.pH + nrcs.code +
                 (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                             nrcs.code == "ACSA3" |
                                             nrcs.code == "QURU" |
                                             nrcs.code == "FAGR" |
                                             nrcs.code == "FRAM2"))

# Check normality assumptions
plot(p.bioe)
qqnorm(residuals(p.bioe))
qqline(residuals(p.bioe))
hist(residuals(p.bioe))
shapiro.test(residuals(p.bioe))
outlierTest(p.bioe)

# Model output
summary(p.bioe)
Anova(p.bioe)
r.squaredGLMM(p.bioe)

# Emmean output for fig making
nbio.pairwise <- data.frame(variable = "p.bioe", 
                            cld(emmeans(p.bioe, ~nrcs.code),
                                Letters = LETTERS))

##########################################################################
## N.photosynthesis - soil N
##########################################################################
data$p.photo[11] <- NA

p.photo <- lmer(p.photo ~ soil.n.norm + mineral.pH +  nrcs.code +
                  (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                              nrcs.code == "ACSA3" |
                                              nrcs.code == "QURU" |
                                              nrcs.code == "FAGR" |
                                              nrcs.code == "FRAM2"))

# Check normality assumptions
plot(p.photo)
qqnorm(residuals(p.photo))
qqline(residuals(p.photo))
hist(residuals(p.photo))
shapiro.test(residuals(p.photo))
outlierTest(p.photo)

# Model output
summary(p.photo)
Anova(p.photo)
r.squaredGLMM(p.photo)

# species pairwise comparison
cld(emmeans(p.photo, pairwise~nrcs.code))

# Emmean output for fig making
npho.pairwise <- data.frame(variable = "p.photo", 
                            cld(emmeans(p.photo, ~nrcs.code),
                                Letters = LETTERS))

##########################################################################
## N in structure
##########################################################################
data$p.structure[49] <- NA

p.structure <- lmer(p.structure ~ soil.n.norm + mineral.pH +  nrcs.code +
                  (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                              nrcs.code == "ACSA3" |
                                              nrcs.code == "QURU" |
                                              nrcs.code == "FAGR" |
                                              nrcs.code == "FRAM2"))

# Check normality assumptions
plot(p.structure)
qqnorm(residuals(p.structure))
qqline(residuals(p.structure))
hist(residuals(p.structure))
shapiro.test(residuals(p.structure))
outlierTest(p.structure)

# Model output
summary(p.structure)
Anova(p.structure)
r.squaredGLMM(p.structure)

# species pairwise comparison
cld(emmeans(p.structure, pairwise~nrcs.code))

# Emmean output for fig making
nstr.pairwise <- data.frame(variable = "p.structure", 
                            cld(emmeans(p.structure, ~nrcs.code),
                                Letters = LETTERS))


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
test(emtrends(a400, ~1, var = "mineral.pH"))
cld(emmeans(a400, pairwise~nrcs.code))

# Emmean output for fig making
a400.pairwise <- data.frame(variable = "a400", 
                            cld(emmeans(a400, ~nrcs.code),
                                Letters = LETTERS))

##########################################################################
## Net photosynthesis-leaf N
##########################################################################
data$a400[data$a400 < 0.2] <- NA

a400.nleaf <- lmer(sqrt(a400) ~ narea + (1 | nrcs.code) +  (1 | site), 
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
emmeans(vcmax, ~1, at = list(soil.n.norm = 0))
emmeans(vcmax, pairwise~nrcs.code)

# Pairwise comparisons
emtrends(vcmax, ~1,  var = "soil.n.norm")
emtrends(vcmax, ~1,  var = "mineral.pH")

# Emmean output for fig making
vcmax.pairwise <- data.frame(variable = "vcmax25", 
                             cld(emmeans(vcmax, ~nrcs.code),
                                 Letters = LETTERS))

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
emmeans(jmax, ~mineral.pH, at = list(soil.n.norm = 0))
emmeans(jmax, pairwise~nrcs.code)

# Emmean output for fig making
jmax.pairwise <- data.frame(variable = "jmax25", 
                            cld(emmeans(jmax, ~nrcs.code),
                                Letters = LETTERS))

##########################################################################
## Jmax25 - leaf N
##########################################################################
jmax.nleaf <- lmer(jmax25 ~ narea + (1 | nrcs.code) +  (1 | site), 
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
emmeans(vjmax, ~1, at = list(soil.n.norm = 0))
emmeans(vjmax, pairwise~nrcs.code)

## Fig int and slope
test(emtrends(vjmax, ~1, var = "soil.n.norm"))
emmeans(vjmax, ~mineral.pH, at = list(soil.n.norm = 0))

# Emmean output for fig making
vjmax.pairwise <- data.frame(variable = "jmax:vcmax", 
                             cld(emmeans(vjmax, ~nrcs.code, type = "response"),
                                 Letters = LETTERS))
names(vjmax.pairwise)[3] <- "emmean"

##########################################################################
## Jmax25:Vcmax25 - leaf N
##########################################################################
vjmax.nleaf <- lmer(log(jmax.vcmax) ~ narea + (1 | nrcs.code) +  (1 | site), 
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
emmeans(chi, ~1, at = list(soil.n.norm = 0))
emmeans(chi, pairwise~nrcs.code)

# Emmean output for fig making
chi.pairwise <- data.frame(variable = "chi", cld(emmeans(chi, ~nrcs.code),
                           Letters = LETTERS))

##########################################################################
## chi - leaf N
##########################################################################
data$chi[66] <- NA

chi.nleaf <- lmer(chi ~ narea + (1 | nrcs.code) +  (1 | site), 
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
test(emtrends(gs.nleaf, ~1, "narea"))

##########################################################################
## Stomatal limitation - soil N
##########################################################################
data$stom.lim[c(68, 75)] <- NA

l <- lmer(stom.lim ~ soil.n.norm + mineral.pH + nrcs.code +
            (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                        nrcs.code == "ACSA3" |
                                        nrcs.code == "QURU" |
                                        nrcs.code == "FAGR" |
                                        nrcs.code == "FRAM2"))

# Check normality assumptions
plot(l)
qqnorm(residuals(l))
qqline(residuals(l))
hist(residuals(l))
shapiro.test(residuals(l))
outlierTest(l)

# Model output
summary(l)
Anova(l)
r.squaredGLMM(l)

# Pairwise comparisons
test(emtrends(l, ~1, var = "soil.n.norm"))
emmeans(l, ~1, at = list(soil.n.norm = 0))
emmeans(l, pairwise~nrcs.code)

# Emmean output for fig making
stomlim.pairwise <- data.frame(variable = "stom.lim",
                               cld(emmeans(l,  ~nrcs.code),
                               Letters = LETTERS))

##########################################################################
## Stomatal limitation - leaf N
##########################################################################
l.nleaf <- lmer(stom.lim ~ narea + (1 | nrcs.code) +  (1 | site), 
                  data = subset(data, nrcs.code == "ACRU" |
                                  nrcs.code == "ACSA3" |
                                  nrcs.code == "QURU" |
                                  nrcs.code == "FAGR" | 
                                  nrcs.code == "FRAM2"))

# Check model assumptions
plot(l.nleaf)
qqnorm(residuals(l.nleaf))
qqline(residuals(l.nleaf))
hist(residuals(l.nleaf))
shapiro.test(residuals(l.nleaf))
outlierTest(l.nleaf)

# Model output
summary(l.nleaf)
Anova(l.nleaf)
r.squaredGLMM(l.nleaf)

# Test Narea-Anet slope
test(emtrends(gs.nleaf, ~1, "narea"))

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
emmeans(pnue, ~1, at = list(soil.n.norm = 0))

test(emtrends(pnue, ~1, var = "mineral.pH"))
emmeans(pnue, pairwise~nrcs.code)

# Fig making intercept and slope
test(emtrends(pnue, ~1, var = "soil.n.norm", regrid = "response"))
emmeans(pnue, ~1, at = list(soil.n.norm = 0), regrid = "response")

# Emmean output for fig making
pnue.pairwise <- data.frame(variable = "pnue", 
                            cld(emmeans(pnue, ~nrcs.code, 
                                        type = "response"),
                                Letters = LETTERS))
names(pnue.pairwise)[3] <- "emmean"


pnue.0 <- data.frame(emmeans(pnue, ~1, at = list(soil.n.norm = 0)))[2]
pnue.27 <- data.frame(emmeans(pnue, ~1, at = list(soil.n.norm = 27.109)))[2]
emmeans(pnue, ~1, at = list(soil.n.norm = 27.109))

## Percent change due to soil N
(pnue.0 - pnue.27) / pnue.0 * 100

##########################################################################
## iWUE - soil N
##########################################################################
data$iwue[data$iwue < 0] <- NA
data$iwue[c(16, 102)] <- NA

iwue <- lmer(iwue ~ soil.n.norm + mineral.pH + nrcs.code +
               (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                           nrcs.code == "ACSA3" |
                                           nrcs.code == "QURU" |
                                           nrcs.code == "FAGR" |
                                           nrcs.code == "FRAM2"))

# Check normality assumptions
plot(iwue)
qqnorm(residuals(iwue))
qqline(residuals(iwue))
hist(residuals(iwue))
shapiro.test(residuals(iwue))
outlierTest(iwue)

# Model output
summary(iwue)
Anova(iwue)
r.squaredGLMM(iwue)

# Pairwise comparisons
emmeans(iwue, pairwise~nrcs.code)

# Emmean output for fig making
iwue.pairwise <- data.frame(variable = "iwue", 
                            cld(emmeans(iwue, ~nrcs.code),
                                Letters = LETTERS))

##########################################################################
## iWUE - leaf N
##########################################################################
iwue.nleaf <- lmer(iwue ~ narea + (1 | nrcs.code) +  (1 | site), 
                data = subset(data, nrcs.code == "ACRU" |
                                nrcs.code == "ACSA3" |
                                nrcs.code == "QURU" |
                                nrcs.code == "FAGR" | 
                                nrcs.code == "FRAM2"))

# Check model assumptions
plot(iwue.nleaf)
qqnorm(residuals(iwue.nleaf))
qqline(residuals(iwue.nleaf))
hist(residuals(iwue.nleaf))
shapiro.test(residuals(iwue.nleaf))
outlierTest(iwue.nleaf)

# Model output
summary(iwue.nleaf)
Anova(iwue.nleaf)
r.squaredGLMM(iwue.nleaf)

# Test Narea-Anet slope
test(emtrends(gs.nleaf, ~1, "narea"))

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
emmeans(narea.chi, ~1, var = "soil.n.norm")
emmeans(narea.chi, pairwise~nrcs.code)

# Emmean output for fig making
narea.chi.pairwise <- data.frame(variable = "narea.chi", 
                                cld(emmeans(narea.chi, ~nrcs.code),
                                    Letters = LETTERS))

## Percent change
narea.chi.0 <- data.frame(emmeans(narea.chi, ~1, at = list(soil.n.norm = 0)))[2]
narea.chi.27 <- data.frame(emmeans(narea.chi, ~1, at = list(soil.n.norm = 27.109)))[2]

emmeans(pnue, pairwise~1, var = "soil.n.norm", 
        at = list(soil.n.norm = c(0, 27.109)))

## Percent change due to soil N
(narea.chi.0 - narea.chi.27) / narea.chi.0 * 100

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

# Emmean output for fig making
vcmax.chi.pairwise <- data.frame(variable = "vcmax.chi", 
                                cld(emmeans(vcmax.chi, ~nrcs.code),
                                    Letters = LETTERS))

##########################################################################
## Vcmax25.chi - leaf N
##########################################################################
vcmax.chi.nleaf <- lmer(vcmax.chi ~ narea + (1 | nrcs.code) +  (1 | site), 
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

##########################################################################
## Create data.frame for all species emmean outputs for figures
##########################################################################
spp.diff <- narea.pairwise %>%
  full_join(leaf.n.pairwise) %>%
  full_join(marea.pairwise) %>%
  full_join(a400.pairwise) %>%
  full_join(vcmax.pairwise) %>%
  full_join(jmax.pairwise) %>%
  full_join(vjmax.pairwise) %>%
  full_join(gs.pairwise) %>%
  full_join(chi.pairwise) %>%
  full_join(stomlim.pairwise) %>%
  full_join(pnue.pairwise) %>%
  full_join(nrub.pairwise) %>%
  full_join(nbio.pairwise) %>%
  full_join(npho.pairwise) %>%
  full_join(nstr.pairwise) %>%
  full_join(iwue.pairwise) %>%
  full_join(narea.chi.pairwise) %>%
  full_join(vcmax.chi.pairwise) %>%
  mutate(.group = trimws(.group, which = "both"))
spp.diff

write.csv(spp.diff, "../data_sheets/NxS_figs_emmeanOutputs.csv")

##########################################################################
## Create Table 1: summary table for Narea, Nmass, and Marea Anova results
##########################################################################
narea.coefs <- data.frame(summary(narea)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.),
         coef.narea = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.narea) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH")) %>%
  print(., row.names = FALSE)
narea.table <- data.frame(Anova(narea)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(narea.coefs) %>%
  dplyr::select(treatment, df = Df, coef.narea, 
                chisq.narea = Chisq, pval.narea = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)",  "soil.n.norm", 
                                       "mineral.pH", "nrcs.code")),
         across(chisq.narea:pval.narea, round, digits = 3),
         chisq.narea = ifelse(chisq.narea < 0.001 & chisq.narea >= 0, 
                             "<0.001", chisq.narea),
         pval.narea = ifelse(pval.narea < 0.001 & pval.narea >= 0, 
                            "<0.001", pval.narea)) %>%
  arrange(treatment)

nmass.coefs <- data.frame(summary(leaf.n)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.),
         coef.nmass = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.nmass) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH")) %>%
  print(., row.names = FALSE)
nmass.table <- data.frame(Anova(leaf.n)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(nmass.coefs) %>%
  dplyr::select(treatment, df = Df, coef.nmass, 
                chisq.nmass = Chisq, pval.nmass = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", "nrcs.code")),
         across(chisq.nmass:pval.nmass, round, digits = 3),
         chisq.nmass = ifelse(chisq.nmass < 0.001 & chisq.nmass >= 0, 
                              "<0.001", chisq.nmass),
         pval.nmass = ifelse(pval.nmass < 0.001 & pval.nmass >= 0, 
                             "<0.001", pval.nmass)) %>%
  arrange(treatment)

marea.coefs <- data.frame(summary(marea)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.),
         coef.marea = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.marea) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH")) %>%
  print(., row.names = FALSE)
marea.table <- data.frame(Anova(marea)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(marea.coefs) %>%
  dplyr::select(treatment, df = Df, coef.marea, 
                chisq.marea = Chisq, pval.marea = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", "nrcs.code")),
         across(chisq.marea:pval.marea, round, digits = 3),
         chisq.marea = ifelse(chisq.marea < 0.001 & chisq.marea >= 0, 
                              "<0.001", chisq.marea),
         pval.marea = ifelse(pval.marea < 0.001 & pval.marea >= 0, 
                             "<0.001", pval.marea)) %>%
  arrange(treatment)

table1 <- narea.table %>% full_join(nmass.table) %>% full_join(marea.table) %>%
  replace(is.na(.), "-")

write.csv(table1, file = "../working_drafts/tables/NxS_table1_leafN.csv",
          row.names = FALSE)

##########################################################################
## Create Table 2: summary table for gas exchange
##########################################################################
anet.coefs <-  data.frame(summary(a400)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(summary(a400.nleaf)$coefficient)) %>%
  mutate(coef.anet = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.anet) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH", "(narea_intercept)", "narea")) %>%
  print(., row.names = FALSE)
anet.table <- data.frame(Anova(a400)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(Anova(a400.nleaf))) %>%
  mutate(treatment = c("soil.n.norm", "mineral.pH", "nrcs.code", "narea")) %>%
  full_join(anet.coefs) %>%
  dplyr::select(treatment, df = Df, coef.anet, 
                chisq.anet = Chisq, pval.anet = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", 
                                       "nrcs.code", "(narea_intercept)", "narea")),
         across(chisq.anet:pval.anet, round, digits = 3),
         chisq.anet = ifelse(chisq.anet < 0.001 & chisq.anet >= 0, 
                              "<0.001", chisq.anet),
         pval.anet = ifelse(pval.anet < 0.001 & pval.anet >= 0, 
                             "<0.001", pval.anet)) %>%
  arrange(treatment)

vcmax.coefs <-  data.frame(summary(vcmax)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(summary(vcmax.nleaf)$coefficient)) %>%
  mutate(coef.vcmax = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.vcmax) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH", "(narea_intercept)", "narea")) %>%
  print(., row.names = FALSE)
vcmax.table <- data.frame(Anova(vcmax)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(Anova(vcmax.nleaf))) %>%
  mutate(treatment = c("soil.n.norm", "mineral.pH", "nrcs.code", "narea")) %>%
  full_join(vcmax.coefs) %>%
  dplyr::select(treatment, df = Df, coef.vcmax, 
                chisq.vcmax = Chisq, pval.vcmax = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", 
                                       "nrcs.code", "(narea_intercept)", "narea")),
         across(chisq.vcmax:pval.vcmax, round, digits = 3),
         chisq.vcmax = ifelse(chisq.vcmax < 0.001 & chisq.vcmax >= 0, 
                             "<0.001", chisq.vcmax),
         pval.vcmax = ifelse(pval.vcmax < 0.001 & pval.vcmax >= 0, 
                            "<0.001", pval.vcmax)) %>%
  arrange(treatment)

jmax.coefs <-  data.frame(summary(jmax)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(summary(jmax.nleaf)$coefficient)) %>%
  mutate(coef.jmax = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.jmax) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH", "(narea_intercept)", "narea")) %>%
  print(., row.names = FALSE)
jmax.table <- data.frame(Anova(jmax)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(Anova(jmax.nleaf))) %>%
  mutate(treatment = c("soil.n.norm", "mineral.pH", "nrcs.code", "narea")) %>%
  full_join(jmax.coefs) %>%
  dplyr::select(treatment, df = Df, coef.jmax, 
                chisq.jmax = Chisq, pval.jmax = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", 
                                       "nrcs.code", "(narea_intercept)", "narea")),
         across(chisq.jmax:pval.jmax, round, digits = 3),
         chisq.jmax = ifelse(chisq.jmax < 0.001 & chisq.jmax >= 0, 
                              "<0.001", chisq.jmax),
         pval.jmax = ifelse(pval.jmax < 0.001 & pval.jmax >= 0, 
                             "<0.001", pval.jmax)) %>%
  arrange(treatment)

jvmax.coefs <-  data.frame(summary(vjmax)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(summary(vjmax.nleaf)$coefficient)) %>%
  mutate(coef.jvmax = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.jvmax) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH", "(narea_intercept)", "narea")) %>%
  print(., row.names = FALSE)
jvmax.table <- data.frame(Anova(vjmax)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(Anova(vjmax.nleaf))) %>%
  mutate(treatment = c("soil.n.norm", "mineral.pH", "nrcs.code", "narea")) %>%
  full_join(jvmax.coefs) %>%
  dplyr::select(treatment, df = Df, coef.jvmax, 
                chisq.jvmax = Chisq, pval.jvmax = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", 
                                       "nrcs.code", "(narea_intercept)", "narea")),
         across(chisq.jvmax:pval.jvmax, round, digits = 3),
         chisq.jvmax = ifelse(chisq.jvmax < 0.001 & chisq.jvmax >= 0, 
                             "<0.001", chisq.jvmax),
         pval.jvmax = ifelse(pval.jvmax < 0.001 & pval.jvmax >= 0, 
                            "<0.001", pval.jvmax)) %>%
  arrange(treatment)

stomlim.coefs <-  data.frame(summary(l)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(summary(l.nleaf)$coefficient)) %>%
  mutate(coef.stomlim = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.stomlim) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH", "(narea_intercept)", "narea")) %>%
  print(., row.names = FALSE)
stomlim.table <- data.frame(Anova(l)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(Anova(l.nleaf))) %>%
  mutate(treatment = c("soil.n.norm", "mineral.pH", "nrcs.code", "narea")) %>%
  full_join(stomlim.coefs) %>%
  dplyr::select(treatment, df = Df, coef.stomlim, 
                chisq.stomlim = Chisq, pval.stomlim = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", 
                                       "nrcs.code", "(narea_intercept)", "narea")),
         across(chisq.stomlim:pval.stomlim, round, digits = 3),
         chisq.stomlim = ifelse(chisq.stomlim < 0.001 & chisq.stomlim >= 0, 
                              "<0.001", chisq.stomlim),
         pval.stomlim = ifelse(pval.stomlim < 0.001 & pval.stomlim >= 0, 
                             "<0.001", pval.stomlim)) %>%
  arrange(treatment)

table2 <- anet.table %>% full_join(vcmax.table) %>% full_join(jmax.table) %>%
  full_join(jvmax.table) %>% full_join(stomlim.table) %>%
  replace(is.na(.), "-")
write.csv(table2,
          "../working_drafts/tables/NxS_table2_leafPhoto.csv", row.names = FALSE)

##########################################################################
## Create Table 3: summary table for prop N in photosynthesis
##########################################################################
p.photo.coefs <- data.frame(summary(p.photo)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.),
         coef.p.photo = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.p.photo) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH")) %>%
  print(., row.names = FALSE)
p.photo.table <- data.frame(Anova(p.photo)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(p.photo.coefs) %>%
  dplyr::select(treatment, df = Df, coef.p.photo, 
                chisq.p.photo = Chisq, pval.p.photo = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", "nrcs.code")),
         across(chisq.p.photo:pval.p.photo, round, digits = 3),
         chisq.p.photo = ifelse(chisq.p.photo < 0.001 & chisq.p.photo >= 0, 
                              "<0.001", chisq.p.photo),
         pval.p.photo = ifelse(pval.p.photo < 0.001 & pval.p.photo >= 0, 
                             "<0.001", pval.p.photo)) %>%
  arrange(treatment)

p.rubisco.coefs <- data.frame(summary(p.rubisco)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.),
         coef.p.rubisco = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.p.rubisco) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH")) %>%
  print(., row.names = FALSE)
p.rubisco.table <- data.frame(Anova(p.rubisco)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(p.rubisco.coefs) %>%
  dplyr::select(treatment, df = Df, coef.p.rubisco, 
                chisq.p.rubisco = Chisq, pval.p.rubisco = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", "nrcs.code")),
         across(chisq.p.rubisco:pval.p.rubisco, round, digits = 3),
         chisq.p.rubisco = ifelse(chisq.p.rubisco < 0.001 & chisq.p.rubisco >= 0, 
                                "<0.001", chisq.p.rubisco),
         pval.p.rubisco = ifelse(pval.p.rubisco < 0.001 & pval.p.rubisco >= 0, 
                               "<0.001", pval.p.rubisco)) %>%
  arrange(treatment)

p.bioe.coefs <- data.frame(summary(p.bioe)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.),
         coef.p.bioe = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.p.bioe) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH")) %>%
  print(., row.names = FALSE)
p.bioe.table <- data.frame(Anova(p.bioe)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(p.bioe.coefs) %>%
  dplyr::select(treatment, df = Df, coef.p.bioe, 
                chisq.p.bioe = Chisq, pval.p.bioe = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", "nrcs.code")),
         across(chisq.p.bioe:pval.p.bioe, round, digits = 3),
         chisq.p.bioe = ifelse(chisq.p.bioe < 0.001 & chisq.p.bioe >= 0, 
                                  "<0.001", chisq.p.bioe),
         pval.p.bioe = ifelse(pval.p.bioe < 0.001 & pval.p.bioe >= 0, 
                                 "<0.001", pval.p.bioe)) %>%
  arrange(treatment)

p.structure.coefs <- data.frame(summary(p.structure)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.),
         coef.p.structure = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.p.structure) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH")) %>%
  print(., row.names = FALSE)
p.structure.table <- data.frame(Anova(p.structure)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(p.structure.coefs) %>%
  dplyr::select(treatment, df = Df, coef.p.structure, 
                chisq.p.structure = Chisq, pval.p.structure = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", "nrcs.code")),
         across(chisq.p.structure:pval.p.structure, round, digits = 3),
         chisq.p.structure = ifelse(chisq.p.structure < 0.001 & chisq.p.structure >= 0, 
                               "<0.001", chisq.p.structure),
         pval.p.structure = ifelse(pval.p.structure < 0.001 & pval.p.structure >= 0, 
                              "<0.001", pval.p.structure)) %>%
  arrange(treatment)

table3 <- p.photo.table %>% full_join(p.rubisco.table) %>% 
  full_join(p.bioe.table) %>% full_join(p.structure.table) %>%
  replace(is.na(.), "-")

write.csv(table3, file = "../working_drafts/tables/NxS_table3_propN.csv",
          row.names = FALSE)

##########################################################################
## Create Table 4: summary table for N-H2O tradeoffs
##########################################################################
chi.coefs <-  data.frame(summary(chi)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(summary(chi.nleaf)$coefficient)) %>%
  mutate(coef.chi = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.chi) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH", "(narea_intercept)", "narea")) %>%
  print(., row.names = FALSE)
chi.table <- data.frame(Anova(chi)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(Anova(chi.nleaf))) %>%
  mutate(treatment = c("soil.n.norm", "mineral.pH", "nrcs.code", "narea")) %>%
  full_join(chi.coefs) %>%
  dplyr::select(treatment, df = Df, coef.chi, 
                chisq.chi = Chisq, pval.chi = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", 
                                       "nrcs.code", "(narea_intercept)", "narea")),
         across(chisq.chi:pval.chi, round, digits = 3),
         chisq.chi = ifelse(chisq.chi < 0.001 & chisq.chi >= 0, 
                             "<0.001", chisq.chi),
         pval.chi = ifelse(pval.chi < 0.001 & pval.chi >= 0, 
                            "<0.001", pval.chi)) %>%
  arrange(treatment)

pnue.coefs <-  data.frame(summary(pnue)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.)) %>%
  mutate(coef.pnue = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.pnue) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH")) %>%
  print(., row.names = FALSE)
pnue.table <- data.frame(Anova(pnue)) %>% 
  mutate(treatment = row.names(.)) %>%
  mutate(treatment = c("soil.n.norm", "mineral.pH", "nrcs.code")) %>%
  full_join(pnue.coefs) %>%
  dplyr::select(treatment, df = Df, coef.pnue, 
                chisq.pnue = Chisq, pval.pnue = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", 
                                       "nrcs.code")),
         across(chisq.pnue:pval.pnue, round, digits = 3),
         chisq.pnue = ifelse(chisq.pnue < 0.001 & chisq.pnue >= 0, 
                            "<0.001", chisq.pnue),
         pval.pnue = ifelse(pval.pnue < 0.001 & pval.pnue >= 0, 
                           "<0.001", pval.pnue)) %>%
  arrange(treatment)

narea.chi.coefs <-  data.frame(summary(narea.chi)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.)) %>%
  mutate(coef.narea.chi = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.narea.chi) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH")) %>%
  print(., row.names = FALSE)
narea.chi.table <- data.frame(Anova(narea.chi)) %>% 
  mutate(treatment = row.names(.)) %>%
  mutate(treatment = c("soil.n.norm", "mineral.pH", "nrcs.code")) %>%
  full_join(narea.chi.coefs) %>%
  dplyr::select(treatment, df = Df, coef.narea.chi, 
                chisq.narea.chi = Chisq, pval.narea.chi = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", 
                                       "nrcs.code")),
         across(chisq.narea.chi:pval.narea.chi, round, digits = 3),
         chisq.narea.chi = ifelse(chisq.narea.chi < 0.001 & chisq.narea.chi >= 0, 
                             "<0.001", chisq.narea.chi),
         pval.narea.chi = ifelse(pval.narea.chi < 0.001 & pval.narea.chi >= 0, 
                            "<0.001", pval.narea.chi)) %>%
  arrange(treatment)

vcmax.chi.coefs <-  data.frame(summary(vcmax.chi)$coefficient[c(1:3),]) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(summary(vcmax.chi.nleaf)$coefficient)) %>%
  mutate(coef.vcmax.chi = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.vcmax.chi) %>%
  mutate(treatment = c("(Intercept)", "soil.n.norm", "mineral.pH", "(narea_intercept)", "narea")) %>%
  print(., row.names = FALSE)
vcmax.chi.table <- data.frame(Anova(vcmax.chi)) %>% 
  mutate(treatment = row.names(.)) %>%
  full_join(data.frame(Anova(vcmax.chi.nleaf))) %>%
  mutate(treatment = c("soil.n.norm", "mineral.pH", "nrcs.code", "narea")) %>%
  full_join(vcmax.chi.coefs) %>%
  dplyr::select(treatment, df = Df, coef.vcmax.chi, 
                chisq.vcmax.chi = Chisq, pval.vcmax.chi = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "soil.n.norm", "mineral.pH", 
                                       "nrcs.code", "(narea_intercept)", "narea")),
         across(chisq.vcmax.chi:pval.vcmax.chi, round, digits = 3),
         chisq.vcmax.chi = ifelse(chisq.vcmax.chi < 0.001 & chisq.vcmax.chi >= 0, 
                            "<0.001", chisq.vcmax.chi),
         pval.vcmax.chi = ifelse(pval.vcmax.chi < 0.001 & pval.vcmax.chi >= 0, 
                           "<0.001", pval.vcmax.chi)) %>%
  arrange(treatment)

table4 <- chi.table %>% full_join(pnue.table) %>% full_join(narea.chi.table) %>%
  full_join(vcmax.chi.table) %>%
  replace(is.na(.), "-")
write.csv(table4,
          "../working_drafts/tables/NxS_table4_pnue_iwue.csv", row.names = FALSE)










    