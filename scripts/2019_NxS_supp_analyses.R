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
library(ggpubr)

emm_options(opt.digits = FALSE)

# Import datasheet
data <- read.csv("../data/2019_NxS_datasheet.csv", stringsAsFactors = FALSE,
                 na.strings = "NA") %>%
  mutate(treatment = factor(treatment, levels = c("C", "NO3", "AS", "S")))

## Central figure theme
pubtheme <- theme_bw(base_size = 16) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(size = 1.5, fill = NA),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 16),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor.y = element_blank(),
        legend.text.align = 0)

##########################################################################
## Import and clean soil dataset
##########################################################################
soils <- data %>% group_by(site, n.trt, s.trt) %>%
  dplyr::summarize(plot.n = mean(soil.n.norm),
                   plot.pH = mean(mineral.pH),
                   plot.no3n = mean(soil.no3n.norm),
                   plot.nh4n = mean(soil.nh4n.norm)) %>%
  mutate(n.trt = factor(n.trt, levels = c("no.n", "n.added")),
         s.trt = factor(s.trt, levels = c("no.s", "s.added")),
         trt = str_c(n.trt, "_", s.trt),
         trt = factor(trt, levels = c("no.n_s.added",
                                      "no.n_no.s",
                                      "n.added_no.s",
                                      "n.added_s.added")))

soils.data.total <- read.csv("../data/2019_NxS_resinbag_datasheet_rep.csv") %>%
  mutate(n.trt = ifelse(treatment == "AS" | treatment == "NO3", "n.added",
                        "no.n"),
         s.trt = ifelse(treatment == "AS" | treatment == "S", "s.added",
                        "no.s"),
         n.trt = factor(n.trt, levels = c("no.n", "n.added")),
         s.trt = factor(s.trt, levels = c("no.s", "s.added")))

##########################################################################
# Leaf temp. effect on Anet - linear regression
##########################################################################
a400 <- lmer(sqrt(a400) ~ leaf.temp + (1 | site), 
             data = subset(data, nrcs.code == "ACSA3"))

head(data)

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
                data = subset(data, nrcs.code == "ACSA3"))
plot(residuals(a400.nls))
coef(a400.nls)

##########################################################################
# Leaf temp. effect on gsw - linear regression
##########################################################################
gs400 <- lmer(log(gsw) ~ leaf.temp + (1 | site), 
              data = subset(data, nrcs.code == "ACSA3"))

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
# Plot for Anet, gsw (Fig. S1)
##########################################################################
# Trendline prep
a400.linear.pred <- data.frame(emmeans(a400, ~1, "leaf.temp",
                                       at = list(leaf.temp = seq(21.7, 31.8, 0.1)),
                                       type = "response"))
a400.nls.pred <- data.frame(emmeans(a400.nls, ~1, "leaf.temp", 
                                    at = list(leaf.temp = seq(21.7, 31.8, 0.1)),
                                    type = "response"))

# Anet plot
a400_plot <- ggplot(data = subset(data, nrcs.code == "ACSA3"),
                 aes(x = leaf.temp, y = a400)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_smooth(data = a400.linear.pred, 
              aes(x = leaf.temp, y = response),
              lty = 2, size = 2, color = "black") +
  geom_smooth(data = a400.nls.pred, 
            aes(x = leaf.temp, y = response),
            lty = 2, size = 2, color = "blue") +
  scale_x_continuous(limits = c(21, 33), breaks = seq(21, 33, 3)) +
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
  labs(x = expression(bold("Leaf temperature ("*degree*"C)")),
       y = expression(bold("A"["net"]~ "(μmol m"^"-2"~"s"^"-1"~")")),
       fill = "Treatment", shape = "Treatment") +
  pubtheme
a400_plot

# Gs plot
gs400.linear.pred <- data.frame(emmeans(gs400, ~1, "leaf.temp",
                                        at = list(leaf.temp = seq(21.7, 31.8, 0.1)),
                                        type = "response"))
gs400.nls.pred <- data.frame(emmeans(gs400.nls, ~1, "leaf.temp", 
                                     at = list(leaf.temp = seq(21.7, 31.8, 0.1)),
                                     type = "response"))
gsw_plot <- ggplot(data = subset(data, nrcs.code == "ACSA3"),
                 aes(x = leaf.temp, y = gsw)) +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, alpha = 0.75) +
  geom_smooth(data = gs400.linear.pred, 
              aes(x = leaf.temp, y = response),
              lty = 2, size = 2, color = "black") +
  geom_smooth(data = gs400.nls.pred, 
              aes(x = leaf.temp, y = response),
              lty = 2, size = 2, color = "blue") +
  scale_x_continuous(limits = c(21, 33), breaks = seq(21, 33, 3)) +
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
  labs(x = expression(bold("Leaf temperature ("*degree*"C)")),
       y = expression(bold("g"["s"]~ "(mol m"^"-2"~"s"^"-1"~")")),
       fill = "Treatment", shape = "Treatment") +
  pubtheme
gsw_plot

png("../../nitrogen_pH/working_drafts/figs/NxS_figS1_leaftemp.png",
    width = 10, height = 3.5, units = "in", res = 600)
ggarrange(a400_plot, gsw_plot, ncol = 2, nrow = 1,
          common.legend = TRUE, legend = "right",
          labels = c("(a)", "(b)", "(c)"),
          font.label = list(size = 18, face = "bold"))
dev.off()


##########################################################################
## Soil N variance across plots (Table S4; Fig S2)
##########################################################################
table_s2 <- soils.data.total %>%
  group_by(site, n.trt, s.trt) %>%
  summarize(plot.n.mean = mean(totalN.normalized, na.rm = TRUE),
            plot.n.count = n(),
            plot.n.se = plot.n.mean / sqrt(plot.n.count),
            plot.n.lci = plot.n.mean - (1.96 * plot.n.se),
            plot.n.uci = plot.n.mean + (1.96 * plot.n.se)) %>%
  mutate(trt = str_c(n.trt, "_", s.trt))

soils.data.total$totalN.normalized[c(3,19)] <- NA

soiln_trt <- soils.data.total %>%
  group_by(n.trt, s.trt) %>%
  summarize(plot.n.mean = mean(totalN.normalized, na.rm = TRUE),
            plot.n.count = n(),
            plot.n.se = plot.n.mean / sqrt(plot.n.count),
            plot.n.lci = plot.n.mean - (1.96 * plot.n.se),
            plot.n.uci = plot.n.mean + (1.96 * plot.n.se),
            plot.no3.mean = mean(no3.normalized, na.rm = TRUE),
            plot.no3.se = plot.no3.mean / sqrt(plot.n.count),
            plot.no3.lci = plot.no3.mean - (1.96 * plot.no3.se),
            plot.no3.uci = plot.no3.mean + (1.96 * plot.no3.se),
            plot.nh4.mean = mean(nh4.normalized, na.rm = TRUE),
            plot.nh4.se = plot.nh4.mean / sqrt(plot.n.count),
            plot.nh4.lci = plot.nh4.mean - (1.96 * plot.nh4.se),
            plot.nh4.uci = plot.nh4.mean + (1.96 * plot.nh4.se)) %>%
  mutate(trt = str_c(n.trt, "_", s.trt),
         trt = factor(trt, levels = c("no.n_s.added", "no.n_no.s", 
                      "n.added_no.s", "n.added_s.added"))) %>%
  dplyr::select(-plot.n.count)

##########################################################################
## Treatment - soil N
##########################################################################
nitrogen <- lmer(plot.n ~ n.trt * s.trt + (1 | site), 
                 data = soils)

# Check model assumptions
plot(nitrogen)
qqnorm(residuals(nitrogen))
qqline(residuals(nitrogen))
hist(residuals(nitrogen))
shapiro.test(residuals(nitrogen))
outlierTest(nitrogen)

# Model output
summary(nitrogen)
Anova(nitrogen)
r.squaredGLMM(nitrogen)

# Post-hoc tests
emmeans(nitrogen, pairwise~n.trt)

# Plot prep
n_compact <- cld(emmeans(nitrogen, pairwise~n.trt*s.trt),
                 Letters = letters) %>% 
  data.frame() %>%
  mutate(.group = trimws(.group, "both"),
         trt = str_c(n.trt, s.trt, sep = "_"),
         trt = factor(trt, levels = c("no.n_s.added", "no.n_no.s", 
                                      "n.added_no.s", "n.added_s.added")))

##########################################################################
## Treatment - soil pH
##########################################################################
pH <- lmer(plot.pH ~ n.trt * s.trt + (1 | site), data = soils)

# Check model assumptions
plot(pH)
qqnorm(residuals(pH))
qqline(residuals(pH))
hist(residuals(pH))
shapiro.test(residuals(pH))
outlierTest(pH)

# Model output
summary(pH)
Anova(pH)
r.squaredGLMM(pH)

# Post-hoc tests
emmeans(pH, pairwise~s.trt*n.trt)
emmeans(pH, pairwise~s.trt)

# Plot prep
pH_compact <- cld(emmeans(pH, pairwise~n.trt*s.trt),
                  Letters = letters, alpha = 0.1) %>% 
  data.frame() %>%
  mutate(.group = trimws(.group, "both"),
         trt = str_c(n.trt, s.trt, sep = "_"),
         trt = factor(trt, levels = c("no.n_s.added", "no.n_no.s", 
                                      "n.added_no.s", "n.added_s.added")))

##########################################################################
## Treatment - soil no3n
##########################################################################
no3n <- lmer(plot.no3n ~ n.trt * s.trt + (1 | site), 
             data = soils)

# Check model assumptions
plot(no3n)
qqnorm(residuals(no3n))
qqline(residuals(no3n))
hist(residuals(no3n))
shapiro.test(residuals(no3n))
outlierTest(no3n)

# Model output
summary(no3n)
Anova(no3n)
r.squaredGLMM(no3n)

# Post-hoc tests
cld(emmeans(no3n, pairwise~n.trt*s.trt))
emmeans(no3n, pairwise~n.trt)

# Plot prep
no3n_compact <- cld(emmeans(no3n, pairwise~n.trt*s.trt),
                    Letters = letters) %>% 
  data.frame() %>%
  mutate(.group = trimws(.group, "both"),
         trt = str_c(n.trt, s.trt, sep = "_"),
         trt = factor(trt, levels = c("no.n_s.added", "no.n_no.s", 
                                      "n.added_no.s", "n.added_s.added")))

##########################################################################
## Treatment - soil nh4
##########################################################################
nh4n <- lmer(plot.nh4n ~ n.trt * s.trt + (1 | site), 
             data = soils)

# Check model assumptions
plot(nh4n)
qqnorm(residuals(nh4n))
qqline(residuals(nh4n))
hist(residuals(nh4n))
shapiro.test(residuals(nh4n))
outlierTest(nh4n)

# Model output
summary(nh4n)
Anova(nh4n)
r.squaredGLMM(nh4n)

# Post-hoc tests
cld(emmeans(nh4n, pairwise~n.trt*s.trt), Letters = letters)
emmeans(nh4n, pairwise~n.trt)

# Plot prep
nh4_compact <- cld(emmeans(nh4n, pairwise~n.trt*s.trt),
                   Letters = letters) %>% 
  data.frame() %>%
  mutate(.group = trimws(.group, "both"),
         trt = str_c(n.trt, s.trt, sep = "_"),
         trt = factor(trt, levels = c("no.n_s.added", "no.n_no.s", 
                                      "n.added_no.s", "n.added_s.added")))

##########################################################################
## Figure S2
##########################################################################
n_plot <- ggplot(data = soiln_trt, 
                 aes(x = trt, y = plot.n.mean)) +
  geom_errorbar(aes(ymin = plot.n.lci, ymax = plot.n.uci),
                width = 0.5, linewidth= 1) +
  geom_point(aes(shape = trt, fill = trt), 
             size = 5) +
  geom_text(data = n_compact, 
            aes(x = trt, y = 30, label = .group),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_x_discrete(labels = c("S", "C", expression("NO"["3"]), "AS")) +
  scale_shape_manual(values = c(24, 21, 22, 23),
                     labels = c("no.n_no.s" = "no N; no S",
                                "n.added_no.s" = "+ N; no S",
                                "n.added_s.added" = "+ N; + S",
                                "no.n_s.added" = "no N; + S")) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#D55E00", "#E69F00"),
                    labels = c("no.n_no.s" = "no N; no S",
                               "n.added_no.s" = "+ N; no S",
                               "n.added_s.added" = "+ N; + S",
                               "no.n_s.added" = "no N; + S")) +
  labs(x = "Treatment",
       y = expression(bold("Soil N ("*mu*"g N g"["resin"]*""^"-1"*" d"^"-1"*")"))) +
  guides(fill = "none", shape = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        panel.grid.minor.y = element_blank())

pH_plot <- ggplot(data = pH_compact, 
                  aes(x = trt, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.5, linewidth = 1) +
  geom_point(aes(fill = trt, shape = trt), size = 5) +
  geom_text(aes(x = trt, y = 5.5, label = .group),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(3.5, 5.5), breaks = seq(3.5, 5.5, 0.5)) +
  scale_x_discrete(labels = c("S", "C", expression("NO"["3"]), "AS")) +
  scale_shape_manual(values = c(24, 21, 22, 23),
                     labels = c("no.n_no.s" = "no N; no S",
                                "n.added_no.s" = "+ N; no S",
                                "n.added_s.added" = "+ N; + S",
                                "no.n_s.added" = "no N; + S")) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#D55E00", "#E69F00"),
                    labels = c("no.n_no.s" = "no N; no S",
                               "n.added_no.s" = "+ N; no S",
                               "n.added_s.added" = "+ N; + S",
                               "no.n_s.added" = "no N; + S")) +
  labs(x = "Treatment", y = "Soil pH") +
  guides(fill = "none", shape = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        panel.grid.minor.y = element_blank())

no3_plot <- ggplot(data = soiln_trt, 
                   aes(x = trt, y = plot.no3.mean)) +
  geom_errorbar(aes(ymin = plot.no3.lci, ymax = plot.no3.uci),
                width = 0.5, linewidth = 1) +
  geom_point(aes(fill = trt, shape = trt), size = 5) +
  geom_text(data = no3n_compact, 
            aes(x = trt, y = 30, label = .group),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_x_discrete(labels = c("S", "C", expression("NO"["3"]), "AS")) +
  scale_shape_manual(values = c(24, 21, 22, 23),
                     labels = c("no.n_no.s" = "no N; no S",
                                "n.added_no.s" = "+ N; no S",
                                "n.added_s.added" = "+ N; + S",
                                "no.n_s.added" = "no N; + S")) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#D55E00", "#E69F00"),
                    labels = c("no.n_no.s" = "no N; no S",
                               "n.added_no.s" = "+ N; no S",
                               "n.added_s.added" = "+ N; + S",
                               "no.n_s.added" = "no N; + S")) +
  labs(x = "Treatment",
       y = expression(bold("Soil NO"["3"]*"-N ("*mu*"g NO"["3"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")"))) +
  guides(fill = "none", shape = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        panel.grid.minor.y = element_blank())

nh4_plot <- ggplot(data = soiln_trt, aes(x = trt, y = plot.nh4.mean)) +
  geom_errorbar(aes(ymin = plot.nh4.lci, ymax = plot.nh4.uci),
                width = 0.5, linewidth = 1) +
  geom_point(aes(fill = trt, shape = trt), size = 5) +
  geom_text(data = nh4_compact, 
            aes(x = trt, y = 10, label = .group),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  scale_x_discrete(labels = c("S", "C", expression("NO"["3"]), "AS")) +
  scale_shape_manual(values = c(24, 21, 22, 23),
                     labels = c("no.n_no.s" = "no N; no S",
                                "n.added_no.s" = "+ N; no S",
                                "n.added_s.added" = "+ N; + S",
                                "no.n_s.added" = "no N; + S")) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#D55E00", "#E69F00"),
                    labels = c("no.n_no.s" = "no N; no S",
                               "n.added_no.s" = "+ N; no S",
                               "n.added_s.added" = "+ N; + S",
                               "no.n_s.added" = "no N; + S")) +
  labs(x = "Treatment",
       y = expression(bold("Soil NH"["4"]*"-N ("*mu*"g NH"["4"]*"-N g"["resin"]*""^"-1"*" d"^"-1"*")"))) +
  guides(fill = "none", shape = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        panel.grid.minor.y = element_blank())


png(filename = "../../nitrogen_pH/working_drafts/figs/NxS_figS2_treatment_soil.png",
    height = 10, width = 12, units = "in", res = 600)
ggpubr::ggarrange(n_plot, pH_plot, no3_plot, nh4_plot,
                  ncol = 2, nrow = 2, align = "hv",
                  labels = c("(a)", "(b)", "(c)", "(d)"),
                  font.label = list(size = 18, face = "bold"))
dev.off()

##########################################################################
##########################################################################
## Isolating effects of soil pH on measured leaf traits in A. saccharum
##########################################################################
##########################################################################

##########################################################################
## Nmass - soil pH for only C and S plots
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
## Leaf mass per area - soil pH for only C and S plots
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
## Narea - soil pH for only C and S plots
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
## Anet,area - soil pH for only C and S plots
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
## Vcmax25 - soil pH for only C and S plots
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
## Jmax25 - soil pH for only C and S plots
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
## Jmax25:Vcmax25 - soil pH for only C and S plots
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
## chi - soil pH for only C and S plots
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
## PNUE - soil pH for only C and S plots
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
## Narea.chi - soil pH for only C and S plots
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
## Vcmax25.chi - soil pH for only C and S plots
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
## Effects of nitrate availability and soil pH on traits
##########################################################################
##########################################################################

##########################################################################
## Nleaf - soil nitrate
##########################################################################
leaf.n <- lmer(leaf.n ~ soil.no3n.norm + mineral.pH + (1 | site), 
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
test(emtrends(leaf.n, ~1, var = "soil.no3n.norm"))

##########################################################################
## Leaf mass per area - soil nitrate
##########################################################################
marea <- lmer(marea ~ soil.no3n.norm + mineral.pH + (1 | site),
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
test(emtrends(marea, ~1, var = "soil.no3n.norm"))

##########################################################################
## Narea - soil nitrate
##########################################################################
narea <- lmer(narea ~ soil.no3n.norm + mineral.pH + (1 | site), 
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
test(emtrends(narea, ~1, var = "soil.no3n.norm"))

##########################################################################
## Anet - soil nitrate
##########################################################################
data$a400[data$a400 < 0.2] <- NA

a400 <- lmer(sqrt(a400) ~ soil.no3n.norm + mineral.pH + (1 | site),
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
## Vcmax25 - soil nitrate
##########################################################################
vcmax <- lmer(vcmax25 ~ soil.no3n.norm + mineral.pH + (1 | site), 
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
test(emtrends(vcmax, ~1, var = "soil.no3n.norm"))

##########################################################################
## Jmax25 - soil nitrate
##########################################################################
jmax <- lmer(jmax25 ~ soil.no3n.norm + mineral.pH + (1 | site),
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
test(emtrends(jmax, ~1, var = "soil.no3n.norm"))

##########################################################################
## Jmax25:Vcmax25 - soil nitrate
##########################################################################
vjmax <- lmer(jmax.vcmax ~ soil.no3n.norm + mineral.pH + (1 | site),
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
## chi - soil nitrate
##########################################################################
chi <- lmer(chi ~ soil.no3n.norm + mineral.pH + (1 | site),
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
test(emtrends(chi,  ~1, var = "soil.no3n.norm"))

##########################################################################
## PNUE - soil nitrate
##########################################################################
data$pnue[data$pnue < 0] <- NA

pnue <- lmer(sqrt(pnue) ~ soil.no3n.norm + mineral.pH + (1 | site),
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

##########################################################################
## Narea.chi - soil nitrate
##########################################################################
narea.chi <- lmer(narea.chi ~ soil.no3n.norm + mineral.pH + (1 | site),
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
test(emtrends(narea.chi, ~1, var = "soil.no3n.norm"))

##########################################################################
## Vcmax25.chi - soil nitrate
##########################################################################
data$vcmax.chi[85] <- NA

vcmax.chi <- lmer(vcmax.chi ~ soil.no3n.norm + mineral.pH + (1 | site),
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
test(emtrends(vcmax.chi, ~1, var = "soil.no3n.norm"))

##########################################################################
##########################################################################
## Effects of ammonium availability and soil pH on traits
##########################################################################
##########################################################################

##########################################################################
## Nleaf - soil ammonium
##########################################################################
leaf.n <- lmer(leaf.n ~ soil.nh4n.norm + mineral.pH + (1 | site), 
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

##########################################################################
## Leaf mass per area - soil ammonium
##########################################################################
data$marea[49] <- NA

marea <- lmer(marea ~ soil.nh4n.norm + mineral.pH + (1 | site),
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
test(emtrends(marea, ~1, var = "soil.nh4n.norm"))
test(emtrends(marea, ~1, var = "mineral.pH"))

##########################################################################
## Narea - soil ammonium
##########################################################################
narea <- lmer(narea ~ soil.nh4n.norm + mineral.pH + (1 | site), 
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
test(emtrends(narea, ~1, var = "soil.nh4n.norm"))

##########################################################################
## Anet - soil ammonium
##########################################################################
data$a400[data$a400 < 0.2] <- NA

a400 <- lmer(sqrt(a400) ~ soil.nh4n.norm + mineral.pH + (1 | site),
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
## Vcmax25 - soil ammonium
##########################################################################
vcmax <- lmer(vcmax25 ~ soil.nh4n.norm + mineral.pH + (1 | site), 
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

##########################################################################
## Jmax25 - soil ammonium
##########################################################################
jmax <- lmer(jmax25 ~ soil.nh4n.norm + mineral.pH + (1 | site),
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

##########################################################################
## Jmax25:Vcmax25 - soil ammonium
##########################################################################
vjmax <- lmer(jmax.vcmax ~ soil.nh4n.norm + mineral.pH + (1 | site),
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
## chi - soil ammonium
##########################################################################
chi <- lmer(chi ~ soil.nh4n.norm + mineral.pH + (1 | site),
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
test(emtrends(chi,  ~1, var = "soil.nh4n.norm"))

##########################################################################
## PNUE - soil ammonium
##########################################################################
data$pnue[data$pnue < 0] <- NA

pnue <- lmer(sqrt(pnue) ~ soil.nh4n.norm + mineral.pH + (1 | site),
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
test(emtrends(pnue, ~1, var = "soil.nh4n.norm"))
test(emtrends(pnue, ~1, var = "mineral.pH"))

##########################################################################
## Narea.chi - soil ammonium
##########################################################################
narea.chi <- lmer(narea.chi ~ soil.nh4n.norm + mineral.pH + (1 | site),
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
test(emtrends(narea.chi, ~1, var = "soil.nh4n.norm"))
test(emtrends(narea.chi, ~1, var = "mineral.pH"))

##########################################################################
## Vcmax25.chi - soil ammonium
##########################################################################
vcmax.chi <- lmer(log(vcmax.chi) ~ soil.nh4n.norm + mineral.pH + (1 | site),
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
