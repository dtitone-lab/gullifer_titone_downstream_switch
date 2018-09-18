# Script corresponding to 
# Gullifer, J. W., & Titone, D. (submitted). The impact of a momentary language 
# switch on bilingual reading. Journal of Experimental Psychology: Learning, Memory, and Cognition.

library(checkpoint)
checkpoint("2017-04-16")

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(effects)
library(lmerTest)
library(devtools)
# devtools::install_github("crsh/papaja") # install papaja package
library(papaja) 

# Load data ----------------------
# alldata used for analyses at switch region
alldata       <- read.csv("raw_data/alldata.csv", encoding = "UTF-8")

# cognatedata used for analyses at the target region corresponding to the
# cognate effect
cognatedata   <- read.csv("raw_data/cognatedata.csv", encoding = "UTF-8")

# homographdata for analyses at the target region corresponding to the
# homograph effect
homographdata <- read.csv("raw_data/homographdata.csv", encoding = "UTF-8")

# contrast coding for the switching effects
contrasts(alldata$switch_noswitch)        <- -1 * contr.sum(2) / 2
contrasts(homographdata$switch_noswitch)  <- -1 * contr.sum(2) / 2 
contrasts(cognatedata$switch_noswitch)    <- -1 * contr.sum(2) / 2

# contrasts for the word type effects
contrasts(homographdata$WordType)  <- -1 * contr.sum(2) / 2
contrasts(cognatedata$WordType)    <- -1 * contr.sum(2) / 2

# Datasets for French L1 bilinguals
alldata.f       <- alldata[alldata$L1 == "French",]
cognatedata.f   <- cognatedata[cognatedata$L1 == "French",]
homographdata.f <- homographdata[homographdata$L1 == "French",]

# Datasets for English L1 bilinguals
alldata.e       <- alldata[alldata$L1 == "English",]
cognatedata.e   <- cognatedata[cognatedata$L1 == "English",]
homographdata.e <- homographdata[homographdata$L1 =="English",]

# Aggregate data for subject-level plots ----------------------
# GD, switch region
switch.sum.gd.subj <- alldata %>%
  group_by(subject, L1, switch_noswitch) %>%
  summarize(GD = mean(FPGD_switch, na.rm = T)
  )

switch.sum.gd.avg <- switch.sum.gd.subj %>%
  group_by(L1, switch_noswitch) %>%
  summarize(N = n(),
            meanRT = mean(GD),
            sd = sd(GD),
            serr = sd(GD) / sqrt(N))%>%
  mutate(measure = "GD", switch_target = "switch")

# TRT, switch region
switch.sum.trt.subj <- alldata %>%
  group_by(subject, L1, switch_noswitch) %>%
  summarize(TRT = mean(TRT_switch, na.rm = T)
  )   

switch.sum.trt.avg <- switch.sum.trt.subj %>%
  group_by(L1, switch_noswitch) %>%
  summarize(N = n(),
            meanRT = mean(TRT),                   
            sd = sd(TRT),           
            serr = sd(TRT) / sqrt(N))%>%
  mutate(measure = "TRT", switch_target = "switch")

# Combine GD and TRT, switch region
switch.sum.avg <- bind_rows(switch.sum.gd.avg, switch.sum.trt.avg)

# GD, target cognate stimuli
cogs.target.sum.gd.subj <- cognatedata %>%
  group_by(subject, L1, switch_noswitch, WordType) %>%
  summarize(GD = mean(FPGD_target, na.rm = T)
  )

cogs.target.sum.gd.avg <- cogs.target.sum.gd.subj %>%
  group_by(L1,switch_noswitch,WordType) %>%
  summarize(N = n(),
            meanRT = mean(GD),
            sd = sd(GD),
            serr = sd(GD) / sqrt(N))%>%
  mutate(measure = "GD", switch_target = "target",dataset = "Cognates")

# TRT, target cognate stimuli
cogs.target.sum.trt.subj <- cognatedata %>%
  group_by(subject, L1, switch_noswitch, WordType) %>%
  summarize(TRT = mean(TRT_target, na.rm = T)
  )    

cogs.target.sum.trt.avg <- cogs.target.sum.trt.subj %>%
  group_by(L1, switch_noswitch, WordType) %>%
  summarize(N = n(),
            meanRT = mean(TRT),                   
            sd = sd(TRT),           
            serr = sd(TRT) / sqrt(N))%>%
  mutate(measure = "TRT", switch_target = "target", dataset = "Cognates")

# Combine the GD and TRT, target cognate stimuli
cogs.target.sum.avg <- bind_rows(cogs.target.sum.gd.avg, cogs.target.sum.trt.avg)

# GD, target homograph stimuli
homs.target.sum.gd.subj <- homographdata %>%
  group_by(subject, L1, switch_noswitch, WordType) %>%
  summarize(GD = mean(FPGD_target, na.rm = T)
  )

homs.target.sum.gd.avg <- homs.target.sum.gd.subj %>%
  group_by(L1, switch_noswitch, WordType) %>%
  summarize(N = n(),
            meanRT = mean(GD),
            sd = sd(GD),
            serr = sd(GD) / sqrt(N))%>%
  mutate(measure = "GD", switch_target = "target", dataset = "Homographs")

# TRT, target homograph stimuli
homs.target.sum.trt.subj <- homographdata %>%
  group_by(subject, L1, switch_noswitch, WordType) %>%
  summarize(TRT = mean(TRT_target, na.rm = T)
  )    

homs.target.sum.trt.avg <- homs.target.sum.trt.subj %>%
  group_by(L1, switch_noswitch, WordType) %>%
  summarize(N = n(),
            meanRT = mean(TRT),                   
            sd = sd(TRT),           
            serr = sd(TRT) / sqrt(N))%>%
  mutate(measure = "TRT", switch_target = "target", dataset = "Homographs")

# Combine the GD and TRT, target homograph stimuli
homs.target.sum.avg <- bind_rows(homs.target.sum.gd.avg, homs.target.sum.trt.avg)

# Combine the aggregate together for plotting
all.target.sum.avg <- bind_rows(cogs.target.sum.avg, homs.target.sum.avg)
all.target.sum.avg$targ_cont <- ifelse(all.target.sum.avg$WordType == "Cognate" | 
                                         all.target.sum.avg$WordType == "Homograph", "Target", "Control")
switch.sum.avg$WordType  = NA
switch.sum.avg$dataset   = "Switch"
switch.sum.avg$targ_cont = NA

plot_data<-bind_rows(switch.sum.avg, all.target.sum.avg)

# Plot Figure 1 ---------------------------------------------------------------------
# Plots data for the two experiments, all conditions
plot_data$l1_measure <- paste(plot_data$L1, plot_data$measure, sep = "_")
plot_data$l1_measure <- factor(plot_data$l1_measure, levels = c("French_GD", "French_TRT", "English_GD",
                                                              "English_TRT"))
plot_data$dataset    <- factor(plot_data$dataset, levels = c("Switch", "Homographs", "Cognates")) 
plot_data$targ_cont  <- factor(plot_data$targ_cont, levels = c("Target", "Control"))

dataset_names <- c(
  Switch     = "Switch region",
  Homographs = "Downstream target region:\nHomographs",
  Cognates   = "Downstream target region:\nCognates"
)

l1_measure_names <- c(
  French_GD   = "Experiment 1: L2 reading\nGaze duration",
  French_TRT  = "Experiment 1: L2 reading\nTotal reading time",
  English_GD  = "Experiment 2: L1 reading\nGaze duration",
  English_TRT = "Experiment 2: L1 reading\nTotal reading time"
)

ggplot(plot_data,
       aes(x = switch_noswitch, y = meanRT, ymin = meanRT-serr, ymax = meanRT+serr, group = targ_cont)) +
  geom_bar(stat = "identity", position = position_dodge(.9), aes(fill = targ_cont, colour = targ_cont)) +
  geom_errorbar(position = position_dodge(.9), width = .2) +
  facet_grid(l1_measure ~ dataset, labeller = labeller(dataset = dataset_names, l1_measure = l1_measure_names)) +
  xlab("Switching")+
  theme_apa() +
  guides(fill = guide_legend(title = "Word Type"))+
  coord_cartesian(ylim = c(225, 575)) +
  theme(plot.title = element_text(size = rel(1), lineheight = .8, face = "bold")) +
  ylab("Mean reading time (in ms)") +
  scale_fill_manual(breaks = c("Target", "Control"),
                    values = c("#808080", "#D3D3D3") , na.value = "#202020") +
  scale_colour_manual(breaks = c("Target","Control"),
                      values = c("#808080", "#D3D3D3"), na.value = "#202020", guide = F)+
  theme(legend.position = "bottom") +
  theme(strip.text.y = element_text(angle = 360))

ggsave("figures/figure 1.png", width = 12)


# Tables of means ------
switch_data <- plot_data[plot_data$switch_target=="switch",]
switch_data <- switch_data %>% select(Group=L1, switch_noswitch, Measure=measure, meanRT, sd )

swd1 <- switch_data  %>% select(-sd) %>% spread(switch_noswitch, meanRT)
swd2 <- switch_data  %>% select(-meanRT) %>% spread(switch_noswitch, sd)
switch_data <- left_join(swd1, swd2, by=c("Group", "Measure"), suffix=c(".mean",".sd"))

switch_data$Nonswitch <- paste0(round(switch_data$Nonswitch.mean), " (", round(switch_data$Nonswitch.sd), ")")
switch_data$Switch <- paste0(round(switch_data$Switch.mean), " (", round(switch_data$Switch.sd), ")")
switch_data <- switch_data %>% select(-contains("."))

write.csv(switch_data, "tables/table 3.csv", row.names = F)


downstream_data <- plot_data[plot_data$switch_target=="target",]
downstream_data <- downstream_data %>% select(Group=L1, switch_noswitch, Measure=measure, targ_cont, WordType=dataset, meanRT, sd )
downstream_data <- downstream_data %>% unite(condition,switch_noswitch, WordType) 

dd1 <- downstream_data %>% select(-sd) %>% spread(condition, meanRT)
dd2 <- downstream_data %>% select(-meanRT) %>% spread(condition, sd)
downstream_data <- left_join(dd1,dd2, by=c("Group","Measure","targ_cont"), suffix=c(".mean",".sd"))

downstream_data$Nonswitch_Cognates <- paste0(round(downstream_data$Nonswitch_Cognates.mean), " (", round(downstream_data$Nonswitch_Cognates.sd), ")")
downstream_data$Nonswitch_Homographs <- paste0(round(downstream_data$Nonswitch_Homographs.mean), " (", round(downstream_data$Nonswitch_Homographs.sd), ")")
downstream_data$Switch_Cognates <- paste0(round(downstream_data$Switch_Cognates.mean), " (", round(downstream_data$Switch_Cognates.sd), ")")
downstream_data$Switch_Homographs <- paste0(round(downstream_data$Switch_Homographs.mean), " (", round(downstream_data$Switch_Homographs.sd), ")")

downstream_data <- downstream_data %>% select(-contains("."))
write.csv(downstream_data,"tables/table 4.csv", row.names = F)


# Experiment 1 analyses ----------------------
# Switch region ----------------------
# Core model, GD ----------------------
switch.gd.base <- lmer(lFPGD_switch ~ cswitch_word_len + ctrial * switch_noswitch +
                         (1 + cswitch_word_len + ctrial + switch_noswitch | subject) + 
                         (1 + ctrial | switch_word), 
                       data = alldata.f)
summary(switch.gd.base)

# Core model, TRT ----------------------
switch.trt.base <- lmer(lTRT_switch ~ cswitch_word_len  + ctrial * switch_noswitch +
                          (1 + cswitch_word_len + ctrial * switch_noswitch | subject) + 
                          (1 + ctrial | switch_word), 
                        data = alldata.f)
summary(switch.trt.base)

# Interactions with L2 exposure, GD ----------------------
switch.gd.exposure <- lmer(lFPGD_switch ~ cswitch_word_len + ctrial * switch_noswitch * ccurrent_exposure_L2 +
                             (1 | subject) + (1 | switch_word), data = alldata.f)
summary(switch.gd.exposure)

# Interactions with L2 exposure, TRT ----------------------
switch.trt.exposure <- lmer(lTRT_switch ~ cswitch_word_len  + ctrial * switch_noswitch*ccurrent_exposure_L2 +
                              (1 |subject) + (1 |switch_word), data = alldata.f)
summary(switch.trt.exposure)

# Add interactions with AoA, GD ----------------------
switch.gd.aoa <- lmer(lFPGD_switch ~ cswitch_word_len + ctrial * switch_noswitch*caoa + ctrial * switch_noswitch * ccurrent_exposure_L2 +
                             (1 | subject) + (1 | switch_word), data = alldata.f)
summary(switch.gd.aoa)

# Add interactions with AoA, TRT ----------------------
switch.trt.aoa <- lmer(lTRT_switch ~ cswitch_word_len  + ctrial * switch_noswitch*caoa + ctrial * switch_noswitch * ccurrent_exposure_L2 +
                              (1 |subject) + (1 |switch_word), data = alldata.f)
summary(switch.trt.aoa)


# Target region, homographs ----------------------
# Core model, GD ----------------------
target.gd.homs.base <- lmer(lFPGD_target ~ cswitch_word_len + ctrial * switch_noswitch * WordType +
                              (1 + ctrial + switch_noswitch + WordType |subject) +
                              (1 + ctrial + switch_noswitch | target_word), data = homographdata.f)
summary(target.gd.homs.base)

# Core model, TRT ----------------------
target.trt.homs.base <- lmer(lTRT_target ~ cswitch_word_len + ctrial * switch_noswitch * WordType +
                               (1 + ctrial + switch_noswitch * WordType |subject) + 
                               (1 + ctrial + switch_noswitch | target_word), data = homographdata.f)
summary(target.trt.homs.base)

# Interactions with L2 exposure, GD ----------------------
target.gd.homs.exposure <- lmer(lFPGD_target ~ cswitch_word_len + 
                                  ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                                  (1 |subject) + (1 | target_word), data = homographdata.f)
summary(target.gd.homs.exposure)

# Interactions with L2 exposure, TRT ----------------------
target.trt.homs.exposure <- lmer(lTRT_target ~ cswitch_word_len + 
                                   ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                                   (1 |subject) + (1 | target_word), data = homographdata.f)
summary(target.trt.homs.exposure)

# Add interactions with AoA, GD ----------------------
target.gd.homs.aoa <- lmer(lFPGD_target ~ cswitch_word_len + 
                                  ctrial * switch_noswitch * WordType * caoa +  ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                                  (1 |subject) + (1 | target_word), data = homographdata.f)
summary(target.gd.homs.aoa)

# Add interactions with AoA, TRT ----------------------
target.trt.homs.aoa <- lmer(lTRT_target ~ cswitch_word_len + 
                                   ctrial * switch_noswitch * WordType * caoa + ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                                   (1 |subject) + (1 | target_word), data = homographdata.f)
summary(target.trt.homs.aoa)

# Plot Figure 2 ----------------------
ef.gd<-as.data.frame(Effect(c("ccurrent_exposure_L2","WordType", "ctrial"), target.gd.homs.exposure, 
                            xlevels = list(ccurrent_exposure_L2 = c(-2, -1, 0, 1, 2),
                                           ctrial = c(-1.6, 0.35, 2.3)),
                            confidence.level = .6827))
ef.gd$measure <- "GD"

ef.trt<-as.data.frame(Effect(c("ccurrent_exposure_L2","WordType","ctrial"), target.trt.homs.exposure, 
                             xlevels = list(ccurrent_exposure_L2 = c(-2, -1, 0, 1, 2),
                                            ctrial = c(-1.6, 0.35, 2.3)),
                             confidence.level = .6827))
ef.trt$measure <- "TRT"

ef <- bind_rows(ef.gd, ef.trt)

trial_names <- c(
  "-1.6" = "First third\nof experiment",
  "0.35" = "Second third\nof experiment",
  "2.3"  = "Last third\nof experiment"
)

measure_names <- c(
  GD   = "Gaze\nduration",
  TRT  = "Total reading\ntime"
)

ggplot(ef, aes(x = ccurrent_exposure_L2, y = fit, ymin = lower, ymax = upper, group = WordType)) + 
  geom_line(aes(colour = WordType)) + geom_ribbon(alpha = .7, aes(fill = WordType))+
  facet_grid(measure ~ ctrial, labeller = labeller(measure = measure_names, ctrial = trial_names)) +
  theme_apa() +
  scale_colour_grey() + scale_fill_grey() +
  coord_cartesian(ylim = c(5.0,6.5)) +
  ylab("Estimated reading time (in log-ms)")+
  xlab("Z-score current L2 exposure")+
  labs(fill = "Word type", colour = "Word type") +
  theme(plot.title = element_text(size = rel(1), lineheight = .8, face = "bold"))+
  theme(legend.position = "top")
ggsave("figures/figure 2.png", width = 6.5)

# Target region, cognates ----------------------
# Core model, GD ----------------------
target.gd.cogs.base <- lmer(lFPGD_target ~ cswitch_word_len + ctrial * switch_noswitch * WordType +
                              (1 + cswitch_word_len + ctrial + switch_noswitch * WordType |subject) + 
                              (1  + ctrial + switch_noswitch | target_word), data = cognatedata.f)
summary(target.gd.cogs.base)

# Core model, TRT ----------------------
target.trt.cogs.base <- lmer(lTRT_target ~  cswitch_word_len + ctrial * switch_noswitch * WordType +
                               (1  + ctrial + switch_noswitch + WordType |subject) + 
                               (1  + ctrial + switch_noswitch | target_word), data = cognatedata.f)
summary(target.trt.cogs.base)

# Interactions with L2 exposure, GD ----------------------
target.gd.cogs.exp <- lmer(lFPGD_target ~ cswitch_word_len + 
                             ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                             (1 |subject) + (1 | target_word), data = cognatedata.f)
summary(target.gd.cogs.exp)

# Interactions with L2 exposure, TRT ----------------------
target.trt.cogs.exposure <- lmer(lTRT_target ~  cswitch_word_len + 
                                   ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                                   (1 |subject) + (1 | target_word), data = cognatedata.f)
summary(target.trt.cogs.exposure)

# Add interactions with AoA, GD ----------------------
target.gd.cogs.aoa <- lmer(lFPGD_target ~ cswitch_word_len + 
                             ctrial * switch_noswitch * WordType * caoa + ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                             (1 |subject) + (1 | target_word), data = cognatedata.f)
summary(target.gd.cogs.aoa)

ef <- data.frame(Effect(c("caoa","ctrial","WordType","switch_noswitch"), target.gd.cogs.aoa))
ggplot(ef, aes(x=ctrial, y=fit, colour=WordType)) + facet_wrap(caoa~switch_noswitch) + geom_line()

# Add interactions with AoA, TRT ----------------------
target.trt.cogs.aoa <- lmer(lTRT_target ~  cswitch_word_len + 
                                   ctrial * switch_noswitch * WordType * caoa +  ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                                   (1 |subject) + (1 | target_word), data = cognatedata.f)
summary(target.trt.cogs.aoa)
ef <- data.frame(Effect(c("caoa","ctrial","WordType","switch_noswitch"), target.trt.cogs.aoa))
ggplot(ef, aes(x=ctrial, y=fit, colour=WordType)) + facet_wrap(caoa~switch_noswitch) + geom_line()

# Plot Figure 3 ----------------------
ef.gd <- as.data.frame(Effect(c("ccurrent_exposure_L2", "WordType"), target.gd.cogs.exp, 
                              xlevels = list(ccurrent_exposure_L2 = c(-2, -1, 0, 1, 2)),
                              confidence.level = .6827))
ef.gd$measure <- "GD"

ef.trt<-as.data.frame(Effect(c("ccurrent_exposure_L2","WordType"), target.trt.cogs.exposure, 
                             xlevels = list(ccurrent_exposure_L2 = c(-2, -1, 0, 1, 2)),
                             confidence.level = .6827))
ef.trt$measure <- "TRT"

ef <- bind_rows(ef.gd, ef.trt)

trial_names <- c(
  "-1.6" = "First third\nof experiment",
  "0.35" = "Second third\nof experiment",
  "2.3"  = "Last third\nof experiment"
)
measure_names <- c(
  GD   = "Gaze\nduration",
  TRT  = "Total reading\ntime"
)

ggplot(ef, aes(x = ccurrent_exposure_L2, y = fit, ymin = lower, ymax = upper, group = WordType)) + 
  geom_line(aes(colour = WordType)) + geom_ribbon(alpha = .7, aes(fill = WordType)) +
  facet_grid(measure ~ ., labeller = labeller(measure = measure_names)) +
  theme_apa() +
  scale_fill_grey() + scale_colour_grey() +
  coord_cartesian(ylim = c(5.0,6.5)) +
  ylab("Estimated reading time (in log-ms)") +
  xlab("Z-score current L2 exposure") +
  labs(fill = "Word type", colour = "Word type") +
  theme(plot.title = element_text(size = rel(1), lineheight = .8, face = "bold")) +
  theme(legend.position = "top")
ggsave("figures/figure 3.png", width = 4.2, height = 4.2)

# Experiment 2 analyses----------------------
# Switch region ----------------------
# Core model, GD ----------------------
switch.gd.base <- lmer(lFPGD_switch ~ cswitch_word_len + ctrial * switch_noswitch +
                         (1 + cswitch_word_len + ctrial * switch_noswitch | subject) + 
                         (1 + ctrial | switch_word), data = alldata.e)
summary(switch.gd.base)

# Core model, TRT ----------------------
switch.trt.base <- lmer(lTRT_switch ~ cswitch_word_len + ctrial * switch_noswitch +
                          (1 + cswitch_word_len + ctrial + switch_noswitch | subject) + 
                          (1 | switch_word), data = alldata.e)
summary(switch.trt.base)

# Plot Figure 4 ----------------------
ef.gd<-as.data.frame(Effect(c("ctrial","switch_noswitch"), switch.gd.base, 
                            xlevels = list(ctrial = c(-1.6,0.35,2.3)),
                            confidence.level = .6827))
ef.gd$measure <- "GD"

ef.trt<-as.data.frame(Effect(c("ctrial","switch_noswitch"), switch.trt.base, 
                             xlevels = list(ctrial = c(-1.6,0.35,2.3)),
                             confidence.level = .6827))
ef.trt$measure <- "TRT"

ef <- bind_rows(ef.gd, ef.trt)

trial_names <- c(
  "-1.6" = "First third\nof experiment",
  "0.35" = "Second third\nof experiment",
  "2.3"  = "Last third\nof experiment"
)

measure_names <- c(
  GD   = "Gaze\nduration",
  TRT  = "Total reading\ntime"
)

ggplot(ef, aes(x = switch_noswitch, y = fit, ymin = lower, ymax = upper)) + 
  geom_bar(stat = "identity", position = position_dodge(.9), fill = "#202020") +
  geom_errorbar(position = position_dodge(.9), width = .2) +
  theme_apa() +
  facet_grid(measure ~ ctrial, labeller = labeller(measure = measure_names,ctrial = trial_names))+
  coord_cartesian(ylim = c(5.0,6.5)) +
  xlab("Switching") +
  ylab("Estimated reading time (in log-ms)") +
  theme(plot.title = element_text(size = rel(1), lineheight = .8, face = "bold"))+
  theme(legend.position = "top")
ggsave("figures/figure 4.png", width = 6.5)


# Interactions with L2 exposure, GD ----------------------
switch.gd.exposure <- lmer(lFPGD_switch ~ cswitch_word_len + ctrial * switch_noswitch * ccurrent_exposure_L2 +
                             (1|subject) + (1|switch_word), data = alldata.e)
summary(switch.gd.exposure)

# Interactions with L2 exposure, TRT ----------------------
switch.trt.exposure <- lmer(lTRT_switch ~ cswitch_word_len + ctrial * switch_noswitch*ccurrent_exposure_L2 +
                              (1|subject) + (1|switch_word), data = alldata.e)
summary(switch.trt.exposure)

# Add interactions with AoA, GD ----------------------
switch.gd.aoa <- lmer(lFPGD_switch ~ cswitch_word_len + ctrial * switch_noswitch * caoa + ctrial * switch_noswitch * ccurrent_exposure_L2 +
                             (1|subject) + (1|switch_word), data = alldata.e)
summary(switch.gd.aoa)

# Add interactions with AoA, TRT ----------------------
switch.trt.aoa <- lmer(lTRT_switch ~ cswitch_word_len + ctrial * switch_noswitch * caoa  + ctrial * switch_noswitch * ccurrent_exposure_L2 +
                              (1|subject) + (1|switch_word), data = alldata.e)
summary(switch.trt.aoa)

# Target region, homographs ----------------------
# Core model, GD ----------------------
target.gd.homs.base <- lmer(lFPGD_target ~ cswitch_word_len + ctrial * switch_noswitch * WordType + 
                              (1 + switch_noswitch + WordType |subject) + (1  + switch_noswitch | target_word), 
                            data = homographdata.e)
summary(target.gd.homs.base)

# Core model, TRT ----------------------
target.trt.homs.base <- lmer(lTRT_target ~ cswitch_word_len + ctrial * switch_noswitch * WordType +
                               (1 + cswitch_word_len + ctrial + switch_noswitch * WordType |subject) + 
                               (1  + ctrial + switch_noswitch | target_word), data = homographdata.e)
summary(target.trt.homs.base)


# Plot Figure 5 ----------------------
ef.gd <- as.data.frame(Effect(c("WordType", "switch_noswitch", "ctrial"), target.gd.homs.base, 
                              xlevels = list(ctrial = c(-1.6,0.35,2.3)),
                              confidence.level = .6827))
ef.gd$measure <- "GD"
ef.trt <- as.data.frame(Effect(c("WordType", "switch_noswitch", "ctrial"), target.trt.homs.base, 
                               xlevels = list(ctrial = c(-1.6, 0.35, 2.3)),
                               confidence.level = .6827))
ef.trt$measure <- "TRT"

ef <- bind_rows(ef.gd,ef.trt)

trial_names <- c(
  "-1.6" = "First third\nof experiment",
  "0.35" = "Second third\nof experiment",
  "2.3"  = "Last third\nof experiment"
)
measure_names <- c(
  GD   = "Gaze\nduration",
  TRT  = "Total reading\ntime"
)

ggplot(ef, aes(x = switch_noswitch, y = fit, ymin = lower, ymax = upper, fill = WordType)) + 
  geom_bar(stat = "identity",position = position_dodge(.9)) +
  geom_errorbar(position = position_dodge(.9), width = .2)+
  facet_grid(measure ~ ctrial, labeller = labeller(measure = measure_names,ctrial = trial_names)) +
  theme_apa() +
  scale_fill_manual(breaks = c("Homograph", "Nonhomograph"),
                    values = c("#808080", "#D3D3D3"), na.value = "#202020") +
  coord_cartesian(ylim = c(5.0,6.5)) +
  ylab("Estimated reading time (in log-ms)")+
  xlab("Switching")+
  labs(fill = "Word type", colour = "Word type") +
  theme(plot.title = element_text(size = rel(1), lineheight = .8, face = "bold")) +
  theme(legend.position = "top")
ggsave("figures/figure 5.png", width = 6.5)

# Interactions with L2 exposure, GD ----------------------
target.gd.homs.exposure <- lmer(lFPGD_target ~ cswitch_word_len + 
                                  ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 + (1 |subject) + 
                                  (1 | target_word), data = homographdata.e)
summary(target.gd.homs.exposure)

# Interactions with L2 exposure, TRT ----------------------
target.trt.homs.exposure <- lmer(lTRT_target ~ cswitch_word_len + 
                                   ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                                   (1 |subject) + (1 | target_word), data = homographdata.e)
summary(target.trt.homs.exposure)

# Add interactions with AoA, GD ----------------------
target.gd.homs.aoa <- lmer(lFPGD_target ~ cswitch_word_len + 
                                  ctrial * switch_noswitch * WordType * caoa +  ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                             (1 |subject) +  (1 | target_word), data = homographdata.e)
summary(target.gd.homs.aoa)

# Add interactions with AoA, TRT ----------------------
target.trt.homs.aoa <- lmer(lTRT_target ~ cswitch_word_len + 
                                   ctrial * switch_noswitch * WordType * caoa + ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                                   (1 |subject) + (1 | target_word), data = homographdata.e)
summary(target.trt.homs.aoa)

# Target region, cognates ----------------------
# Core model, GD ----------------------
target.gd.cogs.base <- lmer(lFPGD_target ~ cswitch_word_len + ctrial * switch_noswitch * WordType +
                              (1 + switch_noswitch + WordType | subject) + (1 + switch_noswitch | target_word),
                            data = cognatedata.e)
summary(target.gd.cogs.base)

# Core model, TRT ----------------------
target.trt.cogs.base <- lmer(lTRT_target ~  cswitch_word_len + ctrial * switch_noswitch * WordType +
                               (1 + ctrial + switch_noswitch * WordType | subject) + 
                               (1 + switch_noswitch | target_word), data = cognatedata.e)
summary(target.trt.cogs.base)

# Interactions with L2 exposure, GD ----------------------
target.gd.cogs.exp <- lmer(lFPGD_target ~ cswitch_word_len + 
                             ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                             (1 |subject) + (1 | target_word), data = cognatedata.e)
summary(target.gd.cogs.exp)

# Interactions with L2 exposure, TRT ----------------------
target.trt.cogs.exposure <- lmer(lTRT_target ~  cswitch_word_len + 
                                   ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                                   (1 |subject) + (1 | target_word), data = cognatedata.e)
summary(target.trt.cogs.exposure)

# Add interactions with AoA, GD ----------------------
target.gd.cogs.aoa <- lmer(lFPGD_target ~ cswitch_word_len + 
                             ctrial * switch_noswitch * WordType * caoa +  ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                             (1 |subject) + (1 | target_word), data = cognatedata.e)
summary(target.gd.cogs.aoa)

# Add interactions with AoA, TRT ----------------------
target.trt.cogs.aoa <- lmer(lTRT_target ~  cswitch_word_len + 
                                   ctrial * switch_noswitch * WordType * caoa +  ctrial * switch_noswitch * WordType * ccurrent_exposure_L2 +
                                   (1 |subject) + (1 | target_word), data = cognatedata.e)
summary(target.trt.cogs.aoa)

