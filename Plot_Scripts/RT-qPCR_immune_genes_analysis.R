library(ggplot2)
library(ggpubr)
library(tidyverse)
library(patchwork)
library(lme4)
library(emmeans)
library(multcomp)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Bonito Lab/Virus_interactions/virus_endophyte_interactions/")

ct_dat <- read.csv("./Data/RT-qPCR_immune_ct_data_Feb7-8_2021.csv", stringsAsFactors = F)

exp_stat <- function(exper, inoculum){
  exper <- as.numeric(exper)
  if (inoculum == "virus"){
    if (exper %% 2 == 0){
      return(1)
    }
    else{
      return(0)
    }
  }
  else if(inoculum == "fungi"){
    if (exper < 5){
      return(1)
    }
    else{
      return(0)
    }
  }
  else if(inoculum == "bact"){
    if (exper < 3){
      return(1)
    }
    else if(exper > 4 && exper < 7){
      return(1)
    }
    else{
      return(0)
    }
  }
}

exp_stat(3, "virus")

ct_dat %>%
  mutate(Sample = gsub("^.*_", "", Name),
         block = substr(Sample, 1, 1),
         exper = substr(Sample, 2, 2),
         virus = unlist(map(exper, exp_stat, inoculum="virus")),
         bact = unlist(map(exper, exp_stat, inoculum="bact")),
         fungi = unlist(map(exper, exp_stat, inoculum="fungi"))) -> ct_dat_mut
ct_dat_mut %>%
  group_by(Name) %>%
  mutate(ct_mean = mean(Ct_SYBR),
            ct_sd = sd(Ct_SYBR)) -> ct_dat_sum
ct_dat_mut %>%
  filter(Target == "Ubi4",
         RT == 1) %>%
  arrange(Sample) %>%
  group_by(Name) %>%
  summarize(ct_mean = mean(Ct_SYBR),
         ct_sd = sd(Ct_SYBR)) -> ref_gene

ct_dat_mut %>%
  filter(Target == "AOS",
         RT == 1) %>%
  arrange(Sample) %>%
  group_by(Name) %>%
  summarize(ct_mean = mean(Ct_SYBR),
            ct_sd = sd(Ct_SYBR)) -> AOS_target
ct_dat_mut %>%
  filter(Target == "W45L1",
         RT == 1) %>%
  arrange(Sample) %>%
  group_by(Name) %>%
  summarize(ct_mean = mean(Ct_SYBR),
            ct_sd = sd(Ct_SYBR)) -> W45L1_target
ct_dat_mut %>%
  filter(Target == "TAR2",
         RT == 1) %>%
  arrange(Sample) %>%
  group_by(Name) %>%
  summarize(ct_mean = mean(Ct_SYBR),
            ct_sd = sd(Ct_SYBR)) -> TAR2_target
TAR2_target

ct_dat_mut %>%
  filter(Target == "Ubi4",
         RT == 1,
         Rep == 1) %>%
  arrange(Sample) -> factorial_dat
all_genes_wide <- data.frame(
  AOS = 2^(ref_gene$ct_mean - AOS_target$ct_mean),
  W45L1 = 2^(ref_gene$ct_mean - W45L1_target$ct_mean),
  TAR2 = 2^(ref_gene$ct_mean - TAR2_target$ct_mean)
)
all_genes_wide <- cbind(all_genes_wide, factorial_dat[,7:12])
all_genes_wide

hist(log(all_genes_wide$AOS))
hist(log(all_genes_wide$W45L1))
hist(log(all_genes_wide$TAR2))

all_genes_wide %>%
  mutate(l_AOS = log(AOS),
         l_W45L1 = log(W45L1),
         l_TAR2 = log(TAR2)) -> all_genes_wide

all_genes_wide
m_AOS <- lmer(l_AOS ~ virus + bact + fungi + (1|block), all_genes_wide)
summary(m_AOS)
m_AOS.emm <- emmeans(m_AOS, ~ virus + bact + fungi)
cld.AOS <- cld(m_AOS.emm, alpha = 0.05, Letters = LETTERS)
cld.AOS

m_W45L1 <- lmer(l_W45L1 ~ virus + bact + fungi + (1|block), all_genes_wide)
summary(m_W45L1)
m_W45L1.emm <- emmeans(m_W45L1, ~ virus + bact + fungi)
cld.W45L1 <- cld(m_W45L1.emm, alpha = 0.05, Letters = LETTERS)
cld.W45L1

m_TAR2 <- lmer(l_TAR2 ~ virus + bact + fungi + (1|block), all_genes_wide)
summary(m_TAR2)
m_TAR2.emm <- emmeans(m_TAR2, ~ virus + bact + fungi)
cld.TAR2 <- cld(m_TAR2.emm, alpha = 0.05, Letters = LETTERS)
cld.TAR2$emmean

all_genes_wide

all_genes_wide %>%
  pivot_longer(cols = AOS:TAR2,
              names_to = "target") -> all_genes_long
all_genes_long %>%
  ggplot(aes(x = block, y = value)) + 
  geom_boxplot() + geom_point() + 
  facet_grid(target~., scales = "free") +
  scale_y_log10()
  
cld.AOS %>%
  ggplot(aes(x = factor(virus),
             color = factor(bact),
             shape = factor(fungi),
             y = exp(emmean))) +
  geom_point(position = position_dodge(width = 0.9),
             size = 3) +
  geom_errorbar(aes(ymin = exp(emmean - SE),
                ymax = exp(emmean + SE)),
                width = 0,
                position = position_dodge(width = 0.9)) +
  geom_point(data = all_genes_wide,
             aes(x = factor(virus),
                 color = factor(bact),
                 shape = factor(fungi),
                 y = AOS),
             position = position_dodge(width = 0.9),
             alpha = 0.7) +
  geom_text(aes(label=.group, y = exp(emmean + 5*SE)),
            # colour = "black",
            position = position_dodge(width = 0.9)) +
  scale_y_log10() +
  theme_pubr() +
  scale_x_discrete(labels=c("-", "+")) +
  scale_color_discrete(labels=c("-", "+")) +
  scale_fill_discrete(labels=c("-", "+")) +
  labs(x = "Virus",
       color = "Bacteria",
       shape = "Fungi",
       y = "Normalized AOS expression") -> AOS_plot
cld.W45L1 %>%
  ggplot(aes(x = factor(virus),
             color = factor(bact),
             shape = factor(fungi),
             y = exp(emmean))) +
  geom_point(position = position_dodge(width = 0.9),
             size = 3) +
  geom_errorbar(aes(ymin = exp(emmean - SE),
                    ymax = exp(emmean + SE)),
                width = 0,
                position = position_dodge(width = 0.9)) +
  geom_point(data = all_genes_wide,
             aes(x = factor(virus),
                 color = factor(bact),
                 shape = factor(fungi),
                 y = W45L1),
             position = position_dodge(width = 0.9),
             alpha = 0.7) +
  geom_text(aes(label=.group, y = exp(emmean + 5*SE)),
            # colour = "black",
            position = position_dodge(width = 0.9)) +
  scale_y_log10() +
  theme_pubr() +
  scale_x_discrete(labels=c("-", "+")) +
  scale_color_discrete(labels=c("-", "+")) +
  scale_fill_discrete(labels=c("-", "+")) +
  labs(x = "Virus",
       color = "Bacteria",
       shape = "Fungi",
       y = "Normalized W45L1 expression") -> W45L1_plot
cld.TAR2 %>%
  ggplot(aes(x = factor(virus),
             color = factor(bact),
             shape = factor(fungi),
             y = exp(emmean))) +
  geom_point(position = position_dodge(width = 0.9),
             size = 3) +
  geom_errorbar(aes(ymin = exp(emmean - SE),
                    ymax = exp(emmean + SE)),
                width = 0,
                position = position_dodge(width = 0.9)) +
  geom_point(data = all_genes_wide,
             aes(x = factor(virus),
                 color = factor(bact),
                 shape = factor(fungi),
                 y = TAR2),
             position = position_dodge(width = 0.9),
             alpha = 0.7) +
  geom_text(aes(label=.group, y = exp(emmean + 5*SE)),
            # colour = "black",
            position = position_dodge(width = 0.9)) +
  scale_y_log10() +
  theme_pubr() +
  scale_x_discrete(labels=c("-", "+")) +
  scale_color_discrete(labels=c("-", "+")) +
  scale_fill_discrete(labels=c("-", "+")) +
  labs(x = "Virus",
       color = "Bacteria",
       shape = "Fungi",
       y = "Normalized TAR2 expression") -> TAR2_plot

g <- AOS_plot + W45L1_plot + TAR2_plot + plot_layout(guides = "collect")

ggsave("./Figures/immune_gene_exp_14dpi.png", g, dpi = 400, units="in", height = 8, width = 10)

all_genes_wide %>%
  ggplot(aes(x = l_AOS, y=l_W45L1)) +
  geom_point()
all_genes_wide %>%
  ggplot(aes(x = l_AOS, y=l_TAR2)) +
  geom_point()
all_genes_wide %>%
  ggplot(aes(x = l_W45L1, y=l_TAR2)) +
  geom_point()
