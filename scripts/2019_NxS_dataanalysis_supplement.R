#################################################
# Load libraries
#################################################
library(lme4)
library(car)
library(emmeans)
library(tidyverse)
library(ggpubr)
library(sjPlot)
library(MuMIn)

#################################################
# Load data, replace negative a400 values with NA
#################################################
data <- read.csv("../data/2019_NxS_datasheet.csv",
                 stringsAsFactors = FALSE,
                 na.strings = "NA")
data$a400[data$a400 < 0.2] <- NA

sub.data <- data %>%
  subset(nrcs.code == "ACRU" | 
           nrcs.code == "ACSA3" |
           nrcs.code == "QURU" | 
           nrcs.code == "FAGR" | 
           nrcs.code == "FRAM2")

#################################################
# Figure theme and color palette
#################################################
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
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor.y = element_blank(),
        legend.text.align = 0)

## Add colorblind friendly palette
cbbPalette <- c("#DDAA33", "#BB5566", "#004488", "#BBBBBB", "#FFFFFF")

#################################################
# Linear mixed-effects model for A400
#################################################
a400 <- lmer(log(a400) ~ leaf.temp + (1 | site), 
             data = subset(data, nrcs.code == "ACRU" |
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

#################################################
# Log polynomial fit for A400
#################################################
a400.nls <- nls(formula = log(a400) ~ a + b*leaf.temp + c*(leaf.temp^2),
                start = list(a = 9.422, b = -0.572, c = 0.0099),
                data = sub.data)

coef(a400.nls)

AICc(a400, a400.nls)
plot(residuals(a400.nls))

#################################################
# Plot for A400 (Fig. S1a)
#################################################
a400.linear <- data.frame(get_model_data(a400, type = "pred"))

a.area <- ggplot(data = sub.data,
                 aes(x = leaf.temp, y = a400)) +
  geom_point(aes(fill = nrcs.code, shape = treatment), 
             size = 3, alpha = 0.65) +
  geom_line(data = a400.linear, aes(x = leaf.temp.x, y = leaf.temp.predicted),
            size = 1, color = "blue", linetype = "dashed") +
  geom_function(fun = function(x) exp(9.422 -0.573*x + 0.010*(x^2)),
                linetype = "dashed", 
                size = 1, xlim = c(21, 32)) +
  scale_x_continuous(limits = c(21, 33), 
                     breaks = seq(21, 33, 3)) +
  scale_y_continuous(limits = c(0, 15), 
                     breaks = seq(0, 15, 5)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         linetype = "none") +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  labs(x = NULL,
       y = expression(bold("A"["400"]~ "(μmol m"^"-2"~"s"^"-1"~")")),
       fill = "Species",
       shape = "Treatment") +
  pubtheme
a.area

#################################################
# Linear mixed-effects model for gs400
#################################################
gs400 <- lmer(log(gsw) ~ leaf.temp + (1 | site), 
             data = subset(data, nrcs.code == "ACRU" |
                             nrcs.code == "ACSA3" |
                             nrcs.code == "QURU" |
                             nrcs.code == "FAGR" |
                             nrcs.code == "FRAM2"))

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
emtrends(gs400, var = "leaf.temp", type = "response")

#################################################
# Log polynomial fit for gs400
#################################################
gs400.nls <- nls(formula = log(gsw) ~ a + b*leaf.temp + c*(leaf.temp^2),
                start = list(a = -0.17, b = -0.18, c = 0.0027),
                data = sub.data)
coef(gs400.nls)
summary(gs400.nls)

AICc(gs400, gs400.nls)

#################################################
# Plot for gs400 (Fig. S1b)
#################################################
gs400.linear <- data.frame(get_model_data(gs400, type = "pred"))

gs.area <- ggplot(data = sub.data,
                  aes(x = leaf.temp, y = gsw)) +
  geom_point(aes(fill = nrcs.code, shape = treatment), 
             size = 3, alpha = 0.65) +
  geom_line(data = gs400.linear, aes(x = leaf.temp.x, y = leaf.temp.predicted),
            size = 1, linetype = "dashed", color = "blue") +
  geom_function(fun = function(x) exp(-0.1695 - 0.1857*x + 0.00266*(x^2)),
                linetype = "dashed", 
                size = 1, xlim = c(21, 32)) +
  scale_x_continuous(limits = c(21, 33), 
                     breaks = seq(21, 33, 3)) +
  scale_y_continuous(limits = c(0, 0.2), 
                     breaks = seq(0, 0.2, 0.05)) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("AS" = "+N, +S",
                                "NO3" = "+N, -S",
                                "S" = "-N, +S",
                                "C" = "-N, -S")) +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         linetype = "none") +
  scale_fill_manual(values = cbbPalette,
                    labels = c("ACRU" = expression(italic("A. rubrum")),
                               "ACSA3" = expression(italic("A. saccharum")),
                               "FRAM2" = expression(italic("F. americana")),
                               "FAGR" = expression(italic("F. grandifolia")),
                               "QURU" = expression(italic("Q. rubra")))) +
  labs(x = NULL,
       y = expression(bold("g"["s400"]~ "(mol m"^"-2"~"s"^"-1"~")")),
       fill = "Species",
       shape = "Treatment") +
  pubtheme
gs.area


#################################################
# Create Fig. S1
#################################################
png("/Users/eaperkowski/git/nitrogen_pH/working_drafts/figs/NxS_figS1_leaftemp.png",
    width = 10, height = 4.5, units = 'in', res = 600)
ggarrange(a.area, gs.area, ncol = 2, nrow = 1,
          common.legend = TRUE,
          align = "hv", legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold")) %>%
  annotate_figure(bottom = text_grob(
    expression(bold("Leaf temperature ("*degree*C*")")),
    size = 18))
dev.off()
