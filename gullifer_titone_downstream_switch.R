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
alldata       <- read.csv("data/alldata.csv")

# cognatedata used for cognate analyses
cognatedata   <- read.csv("data/cognatedata.csv")

# homographdata used for homograph analyses
homographdata <- read.csv("data/homographdata.csv")

# contrasts for the switching effects
contrasts(alldata$switch_noswitch)        <- -1 * contr.sum(2)/2
contrasts(homographdata$switch_noswitch)  <- -1 * contr.sum(2)/2
contrasts(cognatedata$switch_noswitch)    <- -1 * contr.sum(2)/2

# contrasts for the word type effects
contrasts(homographdata$WordType)  <- -1 * contr.sum(2)/2
contrasts(cognatedata$WordType)    <- -1 * contr.sum(2)/2

# Datasets for French L1 bilinguals
alldata.f       <- alldata[alldata$L1 == "French",]
cognatedata.f   <- cognatedata[cognatedata$L1 == "French",]
homographdata.f <- homographdata[homographdata$L1 == "French",]

# Datasets for English L1 bilinguals
alldata.e       <-alldata[alldata$L1 == "English",]
cognatedata.e   <-cognatedata[cognatedata$L1 == "English",]
homographdata.e <-homographdata[homographdata$L1 =="English",]

# Aggregate data for subject-level plots ----------------------
# GD, switch region
switch.sum.gd.subj <- alldata %>%
  group_by(subject, L1, switch_noswitch) %>%
  summarize(GD = mean(FPGD_switch, na.rm=T)
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
  summarize(TRT = mean(TRT_switch, na.rm=T)
  )   

switch.sum.trt.avg <- switch.sum.trt.subj %>%
  group_by(L1, switch_noswitch) %>%
  summarize(N = n(),
            meanRT = mean(TRT),                   
            sd =sd(TRT),           
            serr = sd(TRT) / sqrt(N))%>%
  mutate(measure = "TRT", switch_target = "switch")

# Combine GD and TRT, switch region
switch.sum.avg <- bind_rows(switch.sum.gd.avg, switch.sum.trt.avg)

# GD, target cognate stimuli
cogs.target.sum.gd.subj <- cognatedata %>%
  group_by(subject, L1, switch_noswitch, WordType) %>%
  summarize(GD = mean(FPGD_target, na.rm=T)
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

# GD, targer homograph stimuli
homs.target.sum.gd.subj <- homographdata %>%
  group_by(subject, L1, switch_noswitch, WordType) %>%
  summarize(GD = mean(FPGD_target, na.rm=T)
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

# Combine all data together for plotting
all.target.sum.avg <- bind_rows(cogs.target.sum.avg, homs.target.sum.avg)
all.target.sum.avg$targ_cont <- ifelse(all.target.sum.avg$WordType=="Cognate" | 
                                         all.target.sum.avg$WordType=="Homograph", "Target", "Control")
switch.sum.avg$WordType = NA
switch.sum.avg$dataset = "Switch"
switch.sum.avg$targ_cont = NA

plot_data<-bind_rows(switch.sum.avg, all.target.sum.avg)

# Plot Figure 1 ---------------------------------------------------------------------
# Plots data for each experiment, all condition
plot_data$l1_measure <- paste(plot_data$L1, plot_data$measure, sep="_")
plot_data$l1_measure <- factor(plot_data$l1_measure, levels=c("French_GD", "French_TRT", "English_GD",
                                                              "English_TRT"))
plot_data$dataset    <- factor(plot_data$dataset, levels=c("Switch", "Homographs", "Cognates")) 
plot_data$targ_cont  <- factor(plot_data$targ_cont, levels=c("Target", "Control"))

dataset_names <- c(
  Switch = "Switch region",
  Homographs = "Downstream target region:\nHomographs",
  Cognates = "Downstream target region:\nCognates"
)

l1_measure_names <- c(
  French_GD   = "Experiment 1: L2 reading\nGaze duration",
  French_TRT  = "Experiment 1: L2 reading\nTotal reading time",
  English_GD  = "Experiment 2: L1 reading\nGaze duration",
  English_TRT = "Experiment 2: L1 reading\nTotal reading time"
)

ggplot(plot_data,
       aes(x = switch_noswitch, y = meanRT, ymin = meanRT-serr, ymax = meanRT+serr, group = targ_cont)) +
  geom_bar(stat = "identity", position = position_dodge(.9), aes(fill=targ_cont, colour = targ_cont)) +
  geom_errorbar(position = position_dodge(.9), width=.2) +
  facet_grid(l1_measure ~ dataset, labeller = labeller(dataset = dataset_names, l1_measure = l1_measure_names)) +
  xlab("Switching")+
  theme_apa() +
  guides(fill=guide_legend(title = "Word Type"))+
  coord_cartesian(ylim = c(225, 575)) +
  theme(plot.title = element_text(size = rel(1), lineheight = .8, face="bold")) +
  ylab("Mean reading time (in ms)") +
  scale_fill_manual(breaks = c("Target", "Control"),
                    values = c("#808080", "#D3D3D3") , na.value="#202020") +
  scale_colour_manual(breaks = c("Target","Control"),
                      values = c("#808080", "#D3D3D3"), na.value="#202020", guide=F)+
  theme(legend.position = "bottom") +
  theme(strip.text.y = element_text(angle = 360))

ggsave("figures/figure 1.png", width = 12)


# Experiment 1 analyses ----------------------
# Switch region ----------------------
# Core model, GD ----------------------
switch.gd.base <- lmer(lFPGD_switch ~ cswitch_word_len + ctrial * switch_noswitch +
                         (1 + cswitch_word_len + ctrial + switch_noswitch | subject) + (1 + ctrial | switch_word), 
                       data=alldata.f)
summary(switch.gd.base)

# Core model, TRT ----------------------
switch.trt.base <- lmer(lTRT_switch ~ cswitch_word_len  + ctrial * switch_noswitch +
                          (1 + cswitch_word_len + ctrial * switch_noswitch | subject) + 
                          (1 + ctrial | switch_word), 
                        data=alldata.f)
summary(switch.trt.base)

# Interactions with L2 exposure, GD ----------------------
switch.gd.exposure <- lmer(lFPGD_switch ~ cswitch_word_len + ctrial * switch_noswitch * ccurrent_exposure_L2 +
                             (1 | subject) + (1 | switch_word), data = alldata.f)
summary(switch.gd.exposure)

# Interactions with L2 exposure, TRT ----------------------
switch.trt.exposure <- lmer(lTRT_switch ~ cswitch_word_len  + ctrial * switch_noswitch*ccurrent_exposure_L2 +
                              (1 |subject) + (1 |switch_word), data = alldata.f)
summary(switch.trt.exposure)

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
  facet_grid(measure ~ ctrial, labeller = labeller(measure = measure_names,ctrial = trial_names)) +
  theme_apa() +
  scale_colour_grey() + scale_fill_grey() +
  coord_cartesian(ylim = c(5.0,6.5)) +
  ylab("Estimated reading time (in log-ms)")+
  xlab("Z-score current L2 exposure")+
  labs(fill='Word type', colour='Word type') +
  theme(plot.title = element_text(size = rel(1), lineheight=.8, face="bold"))+
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
  labs(fill='Word type', colour='Word type') +
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
                              confidence.level=.6827))
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
  geom_errorbar(position = position_dodge(.9),width=.2)+
  facet_grid(measure ~ ctrial, labeller = labeller(measure = measure_names,ctrial = trial_names)) +
  theme_apa() +
  scale_fill_manual(breaks = c("Homograph", "Nonhomograph"),
                    values = c("#808080", "#D3D3D3"), na.value = "#202020") +
  coord_cartesian(ylim = c(5.0,6.5)) +
  ylab("Estimated reading time (in log-ms)")+
  xlab("Switching")+
  labs(fill='Word type', colour='Word type') +
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
