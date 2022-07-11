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
data <- read.csv("../data_sheets/NxS_datasheet.csv",
                 stringsAsFactors = FALSE,
                 na.strings = "NA")

n.spp <- data %>%
  filter(nrcs.code == "ACRU" | nrcs.code == "ACSA3" |
           nrcs.code == "QURU" | nrcs.code == "FAGR" | 
           nrcs.code == "FRAM2") %>%
  group_by(site, treatment, nrcs.code) %>%
  summarize(n = n())

## Mean and SE of soil nitrogen availabilityfor discussion
mean(unique(data$soil.n.total.day))
sd(unique(data$soil.n.total.day)) / sqrt(12)

## Mean and SE of soil pH for discussion
mean(unique(data$mineral.pH))
sd(unique(data$mineral.pH)) / sqrt(12)

##########################################################################
## Mean soil N against categorical N and S treatment
##########################################################################
soil.n <- lmer(soil.n.total.day ~ n.trt * s.trt + (1 | site),
               data = data)

# Check model assumptions (no S-W tests due to low n)
plot(soil.n)
hist(residuals(soil.n))
shapiro.test(residuals(soil.n))

# Model output
summary(soil.n)
Anova(soil.n)

# Pairwise comparisons
emmeans(soil.n, pairwise~n.trt)
emmeans(soil.n, pairwise~s.trt)
cld(emmeans(soil.n, pairwise~n.trt * s.trt), Letters = letters)

# Emmean data output
soil.n.pairwise <- data.frame(variable = "soil.n",
                              cld(emmeans(soil.n, 
                                          ~n.trt*s.trt, type = "response"),
                                  Letters = LETTERS)) 
soil.n.pairwise$.group <- trimws(soil.n.pairwise$.group, which = "both")
soil.n.pairwise$.group <- c("B", "A", "A", "C")

##########################################################################
## Mineral pH against categorical N and S treatment
##########################################################################
soil.pH <- lmer(mineral.pH ~ n.trt * s.trt + (1 | site),
                data = data)

# Check model assumptions (no S-W tests due to low n)
plot(soil.pH)
hist(residuals(soil.pH))

# Model output
summary(soil.pH)
Anova(soil.pH)

# Pairwise comparisons
emmeans(soil.pH, pairwise~n.trt)
emmeans(soil.pH, pairwise~s.trt)
cld(emmeans(soil.pH, pairwise~n.trt * s.trt))

soil.pH.pairwise <- data.frame(variable = "soil.pH",
                              cld(emmeans(soil.pH, 
                                          ~n.trt*s.trt),
                                  Letters = LETTERS))
soil.pH.pairwise$.group <- trimws(soil.pH.pairwise$.group, which = "both")
soil.pH.pairwise$.group <- c("C", "C", "B", "A")


soil.pairwise <- soil.n.pairwise %>%
        full_join(soil.pH.pairwise)

write.csv(soil.pairwise, "../data_sheets/NxS_figs_soilemmeanOutputs.csv")

##########################################################################
## Nleaf
##########################################################################
leaf.n <- lmer(leaf.n ~ soil.n.total.day + mineral.pH + nrcs.code + 
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
test(emtrends(leaf.n, ~1, var = "soil.n.total.day"))
emmeans(leaf.n, ~1, at = list(mean.soil.n = 0))
cld(emmeans(leaf.n, pairwise~nrcs.code))


# Emmean output for fig making
leaf.n.pairwise <- data.frame(variable = "leaf.n",
                              cld(emmeans(leaf.n, 
                                          ~nrcs.code),
                                  Letters = LETTERS))

##########################################################################
## SLA
##########################################################################
sla <- lmer(log(sla) ~ soil.n.total.day + mineral.pH + nrcs.code + 
              (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                          nrcs.code == "ACSA3" |
                                          nrcs.code == "QURU" |
                                          nrcs.code == "FAGR" | 
                                          nrcs.code == "FRAM2"))

# Check model assumptions
plot(sla)
qqnorm(residuals(sla))
qqline(residuals(sla))
hist(residuals(sla))
shapiro.test(residuals(sla))
outlierTest(sla)

# Model output
summary(sla)
Anova(sla)
r.squaredGLMM(sla)

# Pairwise comparisons
test(emtrends(sla, ~1, var = "soil.n.total.day"))
cld(emmeans(sla, pairwise~nrcs.code))

## Figure slopes, intercepts, percent change values
emtrends(sla, ~1,var = "soil.n.total.day", transform = "response")
emmeans(sla,~1, at = list(soil.n.total.day = 0), transform = "response")
emmeans(sla, pairwise~nrcs.code, transform = "response")

# Emmean output for fig making
sla.pairwise <- data.frame(variable = "sla",
                           cld(emmeans(sla, ~nrcs.code, type = "response"),
                               Letters = LETTERS))
names(sla.pairwise)[3] <- "emmean"

##########################################################################
## Narea
##########################################################################
data$narea[12] <- NA

narea <- lmer(narea ~ soil.n.total.day + mineral.pH + nrcs.code +
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
test(emtrends(narea, ~1, var = "soil.n.total.day"))
emmeans(narea, ~1, at = list(soil.n.total.day = 0))
cld(emmeans(narea, pairwise~nrcs.code))

# Emmean output for fig making
narea.pairwise <- data.frame(variable = "narea",
                             cld(emmeans(narea, ~nrcs.code),
                                 Letters = LETTERS))

##########################################################################
## Net photosynthesis (area basis)
##########################################################################
data$a400[data$a400 < 0.2] <- NA

a400 <- lmer(log(a400) ~ soil.n.total.day + mineral.pH + nrcs.code + 
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
test(emtrends(a400, ~1, var = "soil.n.total.day"))
test(emtrends(a400, ~1, var = "mineral.pH"))
cld(emmeans(a400, pairwise~nrcs.code))

# Emmean output for fig making
a400.pairwise <- data.frame(variable = "a400", 
                            cld(emmeans(a400, ~nrcs.code, type = "response"),
                                Letters = LETTERS))
names(a400.pairwise)[3] <- "emmean"

##########################################################################
## Vcmax.area standardized to 25degC 
##########################################################################
data$vcmax25[11] <- NA

vcmax <- lmer(sqrt(vcmax25) ~ soil.n.total.day + mineral.pH + nrcs.code + 
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
test(emtrends(vcmax, ~1, var = "soil.n.total.day"))
emmeans(vcmax, ~1, at = list(soil.n.total.day = 0))
cld(emmeans(vcmax, pairwise~nrcs.code))

# Pairwise comparisons
emtrends(vcmax, ~1,  var = "soil.n.total.day", transform = "response")
emmeans(vcmax, ~1, at = list(soil.n.total.day = 0), transform = "response")

# Emmean output for fig making
vcmax.pairwise <- data.frame(variable = "vcmax25", 
                             cld(emmeans(vcmax, ~nrcs.code, type = "response"),
                                 Letters = LETTERS))
names(vcmax.pairwise)[3] <- "emmean"

##########################################################################
## Jmax standardized to 25degC
##########################################################################
data$jmax25[11] <- NA

jmax <- lmer(jmax25 ~ soil.n.total.day + mineral.pH + nrcs.code + 
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
test(emtrends(jmax, ~1, var = "soil.n.total.day"))
emmeans(jmax, ~mineral.pH, at = list(mean.soil.n = 0))
cld(emmeans(jmax, pairwise~nrcs.code))

# Emmean output for fig making
jmax.pairwise <- data.frame(variable = "jmax25", 
                            cld(emmeans(jmax, ~nrcs.code),
                                Letters = LETTERS))

##########################################################################
## Jmax25:Vcmax25
##########################################################################
vjmax <- lmer(log(jmax.vcmax) ~ soil.n.total.day + mineral.pH + nrcs.code + 
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
test(emtrends(vjmax, ~1, var = "soil.n.total.day"))
emmeans(vjmax, ~1, at = list(soil.n.total.day = 0))
emmeans(vjmax, pairwise~nrcs.code)

## Fig int and slope
test(emtrends(vjmax, ~1, var = "soil.n.total.day", transform = "response"))
emmeans(vjmax, ~mineral.pH, at = list(soil.n.total.day = 0), transform = "response")

# Emmean output for fig making
vjmax.pairwise <- data.frame(variable = "jmax:vcmax", 
                             cld(emmeans(vjmax, ~nrcs.code, type = "response"),
                                 Letters = LETTERS))
names(vjmax.pairwise)[3] <- "emmean"


##########################################################################
## Vcmax as a function of Narea 
##########################################################################
vcmax.fx.narea <- lmer(vcmax25 ~ narea + (1 | site), 
                       data = subset(data, nrcs.code == "ACRU" |
                                       nrcs.code == "ACSA3" |
                                       nrcs.code == "QURU" |
                                       nrcs.code == "FAGR" | 
                                       nrcs.code == "FRAM2"))

# Check normality assumptions
plot(vcmax.fx.narea)
qqnorm(residuals(vcmax.fx.narea))
qqline(residuals(vcmax.fx.narea))
hist(residuals(vcmax.fx.narea))
shapiro.test(residuals(vcmax.fx.narea))
outlierTest(vcmax.fx.narea)

# Model output
summary(vcmax.fx.narea)
Anova(vcmax.fx.narea)
r.squaredGLMM(vcmax.fx.narea)

## Pairwise comparisons
test(emtrends(vjmax, ~1, var = "soil.n.total.day"))
emmeans(vjmax, ~1, at = list(soil.n.total.day = 0))
emmeans(vjmax, pairwise~nrcs.code)


##########################################################################
## stomatal conductance (gs)
##########################################################################
gs <- lmer(log(gsw) ~ soil.n.total.day + mineral.pH + nrcs.code + 
             (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                         nrcs.code == "ACSA3" |
                                         nrcs.code == "QURU" |
                                         nrcs.code == "FAGR" | 
                                         nrcs.code == "FRAM2"))

# Check normality assumptions
plot(gs)
qqnorm(residuals(gs))
qqline(residuals(gs))
hist(residuals(gs))
shapiro.test(residuals(gs))
outlierTest(gs)

# Model output
summary(gs)
Anova(gs)
r.squaredGLMM(gs)

## Pairwise comparisons
test(emtrends(gs, ~1, var = "soil.n.total.day"))
emmeans(gs, pairwise~nrcs.code, type = "response")

## Figure slopes and intercepts
emtrends(gs, ~1, var = "soil.n.total.day", transform = "response")
emmeans(gs, ~1, at = list(soil.n.total.day = 0), transform = "response")

# Emmean output for fig making
gs.pairwise <- data.frame(variable = "gs", 
                          cld(emmeans(gs, ~nrcs.code, type = "response"),
                              Letters = LETTERS))
names(gs.pairwise)[3] <- "emmean"

##########################################################################
## chi
##########################################################################
chi <- lmer(chi ~ soil.n.total.day + mineral.pH + nrcs.code +
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
test(emtrends(chi,  ~1, var = "soil.n.total.day"))
emmeans(chi, ~1, at = list(soil.n.total.day = 0))
emmeans(chi, pairwise~nrcs.code)

# Emmean output for fig making
chi.pairwise <- data.frame(variable = "chi", 
                           cld(emmeans(chi, ~nrcs.code),
                               Letters = LETTERS))

##########################################################################
## Stomatal Limitation
##########################################################################
data$stom.lim[c(10, 23, 75)] <- NA

l <- lmer(stom.lim ~ soil.n.total.day + mineral.pH + nrcs.code +
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
test(emtrends(l, ~1, var = "soil.n.total.day"))
emmeans(l, ~1, at = list(soil.n.total.day = 0))
emmeans(l, pairwise~nrcs.code)

# Emmean output for fig making
stomlim.pairwise <- data.frame(variable = "stom.lim",
                               cld(emmeans(l,  ~nrcs.code),
                               Letters = LETTERS))

##########################################################################
## Photosynthetic nitrogen-use efficiency
##########################################################################
data$pnue[data$pnue < 0] <- NA
data$pnue[c(16, 102)] <- NA

pnue <- lmer(sqrt(pnue) ~ soil.n.total.day + mineral.pH + nrcs.code +
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
test(emtrends(pnue, ~1, var = "soil.n.total.day"))
emmeans(pnue, ~1, at = list(soil.n.total.day = 0))

test(emtrends(pnue, ~1, var = "mineral.pH"))
emmeans(pnue, pairwise~nrcs.code)

# Fig making intercept and slope
emtrends(pnue, ~1, var = "soil.n.total.day", transform = "response")
emmeans(pnue, ~1, at = list(soil.n.total.day = 0), transform = "response")

# Emmean output for fig making
pnue.pairwise <- data.frame(variable = "pnue", 
                            cld(emmeans(pnue, ~nrcs.code, type = "response"),
                                Letters = LETTERS))
names(pnue.pairwise)[3] <- "emmean"

##########################################################################
## Water-use efficiency (gas exchange)
##########################################################################
# Remove outliers (MEM residual Bonferroni correction p<0.05)
data$iwue[data$iwue < 0] <- NA
data$iwue[c(16, 102)] <- NA

iwue <- lmer(iwue ~ soil.n.total.day + mineral.pH + nrcs.code +
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
## Narea-gs
##########################################################################
narea.gs <- lmer(log(narea.gs) ~ soil.n.total.day + mineral.pH + nrcs.code +
                   (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                               nrcs.code == "ACSA3" |
                                               nrcs.code == "QURU" |
                                               nrcs.code == "FAGR" |
                                               nrcs.code == "FRAM2"))

# Check normality assumptions
plot(narea.gs)
qqnorm(residuals(narea.gs))
qqline(residuals(narea.gs))
hist(residuals(narea.gs))
shapiro.test(residuals(narea.gs))
outlierTest(narea.gs)

# Model output
summary(narea.gs)
Anova(narea.gs)
r.squaredGLMM(narea.gs)

# Post-hoc tests
test(emtrends(narea.gs, ~1, var = "soil.n.total.day"))
emmeans(narea.gs, ~1, var = "soil.n.total.day")
emmeans(narea.gs, pairwise~nrcs.code)

# Emmean output for fig making
narea.gs.pairwise <- data.frame(variable = "narea.gs", 
                                cld(emmeans(narea.gs, ~nrcs.code,
                                            type = "response"),
                                    Letters = LETTERS))
names(narea.gs.pairwise)[3] <- "emmean"

##########################################################################
## Vcmax-gs
##########################################################################
data$vcmax.gs[11] <- NA

vcmax.gs <- lmer(log(vcmax.gs) ~ soil.n.total.day + mineral.pH + nrcs.code +
                   (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                               nrcs.code == "ACSA3" |
                                               nrcs.code == "QURU" |
                                               nrcs.code == "FAGR" |
                                               nrcs.code == "FRAM2"))

# Check normality assumptions
plot(vcmax.gs)
qqnorm(residuals(vcmax.gs))
qqline(residuals(vcmax.gs))
hist(residuals(vcmax.gs))
shapiro.test(residuals(vcmax.gs))
outlierTest(vcmax.gs)

# Model output
summary(vcmax.gs)
Anova(vcmax.gs)
r.squaredGLMM(vcmax.gs)

# Pairwise comparisons
test(emtrends(vcmax.gs, ~mineral.pH, var = "soil.n.total.day"))
emmeans(vcmax.gs, pairwise~nrcs.code)

## Figure slopes and intercepts
emtrends(vcmax.gs, ~1, var = "soil.n.total.day", transform = "response")
emmeans(vcmax.gs, ~1, at = list(soil.n.total.day = 0), transform = "response")

# Emmean output for fig making
vcmax.gs.pairwise <- data.frame(variable = "vcmax.gs", 
                                cld(emmeans(vcmax.gs, ~nrcs.code, 
                                            type = "response"),
                                    Letters = LETTERS))
names(vcmax.gs.pairwise)[3] <- "emmean"

##########################################################################
## Rate of change in basal area between 2011 and 2019
##########################################################################
data$ba.2011.2019[data$ba.2011.2019 <= 0] <- NA

ba <- lmer(sqrt(ba.2011.2019) ~ soil.n.total.day + mineral.pH + nrcs.code + 
             (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                         nrcs.code == "ACSA3" |
                                         nrcs.code == "QURU" |
                                         nrcs.code == "FAGR" |
                                         nrcs.code == "FRAM2"))

# Check normality assumptions
plot(ba)
qqnorm(residuals(ba))
qqline(residuals(ba))
hist(residuals(ba))
shapiro.test(residuals(ba))
outlierTest(ba)

# Model results
summary(ba)
Anova(ba)
r.squaredGLMM(ba)

## Pairwise comparisons
test(emtrends(ba, ~1, var = "soil.n.total.day"))
cld(emmeans(ba, pairwise~nrcs.code))

## Figure making
emtrends(ba, ~1, var = "soil.n.total.day", transform = "response")
emmeans(ba, ~1, at = list(soil.n.total.day = 0), transform = "response")
emmeans(ba, pairwise~nrcs.code)

# Emmean output for fig making
ba.pairwise <- data.frame(variable = "basal.area",
                          cld(emmeans(ba, ~nrcs.code, type = "response"),
                              Letters = LETTERS))
names(ba.pairwise)[3] <- "emmean"

##########################################################################
## Biomass since 2011
##########################################################################
## Replace outlier and zero values with "NA"
data$growth.2011.2019[data$growth.2011.2019 <= 0] <- NA

growth <- lmer(sqrt(growth.2011.2019) ~ soil.n.total.day + mineral.pH + nrcs.code +
                 (1 | site), data = subset(data, nrcs.code == "ACRU" |
                                             nrcs.code == "ACSA3" |
                                             nrcs.code == "QURU" |
                                             nrcs.code == "FAGR" |
                                             nrcs.code == "FRAM2"))

# Visualize data and check normality assumptions
plot(growth)
qqnorm(residuals(growth))
qqline(residuals(growth))
hist(residuals(growth))
shapiro.test(residuals(growth))
outlierTest(growth)

# Model results
summary(growth)
Anova(growth)
r.squaredGLMM(growth)

# Post-hoc analyses
test(emtrends(growth, ~1, var = 'mineral.pH'))
cld(emmeans(growth, pairwise~nrcs.code))

# Emmean output for fig making
growth.pairwise <- data.frame(variable = "rgr",
                              cld(emmeans(growth,  ~nrcs.code, type = "response"),
                                  Letters = LETTERS))
names(growth.pairwise)[3] <- "emmean"

##########################################################################
## Create data.frame for all species emmean outputs for figures
##########################################################################
spp.diff <- narea.pairwise %>%
  full_join(leaf.n.pairwise) %>%
  full_join(sla.pairwise) %>%
  full_join(a400.pairwise) %>%
  full_join(vcmax.pairwise) %>%
  full_join(jmax.pairwise) %>%
  full_join(vjmax.pairwise) %>%
  full_join(gs.pairwise) %>%
  full_join(chi.pairwise) %>%
  full_join(stomlim.pairwise) %>%
  full_join(pnue.pairwise) %>%
  full_join(iwue.pairwise) %>%
  full_join(narea.gs.pairwise) %>%
  full_join(vcmax.gs.pairwise) %>%
  full_join(ba.pairwise) %>%
  full_join(growth.pairwise) %>%
  mutate(.group = trimws(.group, which = "both"))
spp.diff

write.csv(spp.diff, "../data_sheets/NxS_figs_emmeanOutputs.csv")
                       