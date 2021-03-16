rm(list=ls())
library(ggpubr)
library(ggplot2)
library(dplyr)
library(GGally)
library(gridExtra)
library(patchwork)
setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Bonito Lab/Virus_interactions")
plants <- read.csv("VxE_plant_data_Nov52019.csv")
plants <- plants[,1:9]
plants$root_shoot_ratio <- plants$Dry_root_mass/plants$Dry_shoot_mass
plants$area_height_ratio <- plants$Leaf_area/plants$Height
plants$area_shoot_ratio <- plants$Leaf_area/plants$Dry_shoot_mass
plants$Virus <- as.factor(plants$Virus)
plants$Endophyte <- as.factor(plants$Endophyte)

hist(log(plants$area_height_ratio))

LA_bp <- ggplot(plants, aes(x=Virus,y=Leaf_area,color=Endophyte)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.4, size = 3) +
  stat_summary(geom = "point", size = 5, position=position_dodge(width=0.75)) +
  theme_pubr() +
  scale_x_discrete(labels=c("-", "+")) +
  scale_color_discrete(labels=c("-", "+")) +
  labs(y=expression(Leaf ~ Area ~ (cm^{2})), x=NULL, tag = "A") + 
  stat_compare_means(aes(group = Endophyte), label = "p.signif", method = "t.test") +
  theme(legend.position = "none")
LA_bp

H_bp <- ggplot(plants, aes(x=Virus,y=Height,color=Endophyte)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.4, size = 3) +
  stat_summary(geom = "point", size = 5, position=position_dodge(width=0.75)) +
  theme_pubr() +
  scale_x_discrete(labels=c("-", "+")) +
  scale_color_discrete(labels=c("-", "+")) +
  labs(y=expression(Height ~ (cm)), x=NULL, tag = "B") + 
  stat_compare_means(aes(group = Endophyte), label = "p.signif", method = "t.test") +
  theme(legend.position = "none")
H_bp

DSM_bp <- ggplot(plants, aes(x=Virus,y=Dry_shoot_mass,color=Endophyte)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.4, size = 3) +
  stat_summary(geom = "point", size = 5, position=position_dodge(width=0.75)) +
  theme_pubr() +
  scale_x_discrete(labels=c("-", "+")) +
  scale_color_discrete(labels=c("-", "+")) +
  labs(y=expression(Dry ~ Shoot ~ Mass ~ (g)), x=NULL, tag = "C") + 
  stat_compare_means(aes(group = Endophyte), label = "p.signif", method = "t.test") +
  theme(legend.position = "none")
DSM_bp

DRM_bp <- ggplot(plants, aes(x=Virus,y=Dry_root_mass,color=Endophyte)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.4, size = 3) +
  stat_summary(geom = "point", size = 5, position=position_dodge(width=0.75)) +
  theme_pubr() +
  scale_x_discrete(labels=c("-", "+")) +
  scale_color_discrete(labels=c("-", "+")) +
  labs(y=expression(Dry ~ Root ~ Mass ~ (g)), tag = "D", x = NULL) +
  stat_compare_means(aes(group = Endophyte), label = "p.signif", method = "t.test") +
  theme(legend.position = "none")
DRM_bp

RSR_bp <- ggplot(plants, aes(x=Virus,y=root_shoot_ratio,color=Endophyte)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.4, size = 3) +
  stat_summary(geom = "point", size = 5, position=position_dodge(width=0.75)) + 
  theme_pubr() +
  scale_x_discrete(labels=c("-", "+")) +
  scale_color_discrete(labels=c("-", "+")) +
  labs(y="Root : shoot Ratio", tag = "E") + 
  stat_compare_means(aes(group = Endophyte), label = "p.signif", method = "t.test") +
  theme(legend.position = "none")
RSR_bp

AHR_bp <- ggplot(plants, aes(x=Virus,y=area_height_ratio,color=Endophyte)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.4, size = 3) +
  stat_summary(geom = "point", size = 5, position=position_dodge(width=0.75)) +
  theme_pubr() +
  scale_x_discrete(labels=c("-", "+")) +
  scale_color_discrete(labels=c("-", "+")) +
  labs(y="Leaf Area : Height Ratio", tag = "F") + 
  stat_compare_means(aes(group = Endophyte), label = "p.signif", method = "t.test") +
  theme(legend.position = "none")
AHR_bp

LSR_bp <- ggplot(plants, aes(x=Virus,y=area_shoot_ratio,color=Endophyte)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.4, size = 3) +
  stat_summary(geom = "point", size = 5, position=position_dodge(width=0.75)) +
  theme_pubr() +
  scale_x_discrete(labels=c("-", "+")) +
  scale_color_discrete(labels=c("-", "+")) +
  labs(y=expression(Leaf ~ Area ~ per ~ Shoot ~ mass ~ (cm^{2})/g), x="Virus", tag = "A") + 
  stat_compare_means(aes(group = Endophyte), label = "p.signif", method = "t.test") #+
  # theme(legend.position = "none")
LSR_bp

LSR_bp <- ggplot(plants, aes(x=Endophyte,y=area_shoot_ratio,color=Virus)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.4, size = 3) +
  stat_summary(geom = "point", size = 5, position=position_dodge(width=0.75)) +
  theme_pubr() +
  scale_x_discrete(labels=c("-", "+")) +
  scale_color_discrete(labels=c("-", "+")) +
  labs(y=expression(Leaf ~ Area ~ per ~ Shoot ~ mass ~ (cm^{2})/g), x="Endophyte", tag = "A") + 
  stat_compare_means(aes(group = Virus), label = "p.signif", method = "t.test") #+
# theme(legend.position = "none")
LSR_bp

g <- (LA_bp | H_bp ) / (DSM_bp | DRM_bp ) / (RSR_bp | AHR_bp) &
  theme(legend.position = "right")
g <- g + plot_layout(guides="collect")
g
# ggsave("Plant_character_boxplots.png", g, width=10, height = 6, units="in")
ggsave("Plant_character_boxplots_vert.png", g, width=6, height = 8, units="in")

model <- aov(sqrt(Leaf_area) ~ Endophyte * Virus + Soil_depth, plants)
summary(model)

model <- aov(log(Dry_shoot_mass) ~ Endophyte * Virus + Soil_depth, plants)
summary(model)

model <- aov(Height ~ Endophyte * Virus + Soil_depth, plants)
summary(model)

model <- aov(sqrt(Dry_root_mass) ~ Endophyte * Virus + Soil_depth, plants)
summary(model)

model <- aov(sqrt(root_shoot_ratio) ~ Endophyte * Virus + Soil_depth, plants)
summary(model)

model <- aov(log(area_height_ratio) ~ Endophyte * Virus + Soil_depth, plants)
summary(model)

model <- lm(sqrt(Leaf_area) ~ Endophyte * Virus + Soil_depth, plants)
summary(model)

model <- lm(log(Dry_shoot_mass) ~ Endophyte * Virus + Soil_depth, plants)
summary(model)

model <- lm(Height ~ Endophyte * Virus + Soil_depth, plants)
summary(model)

model <- lm(sqrt(Dry_root_mass) ~ Endophyte * Virus + Soil_depth, plants)
summary(model)

model <- glm(sqrt(Dry_root_mass) ~ Endophyte * Virus + Soil_depth, data=plants)
summary(model)

boxplot(Height ~ Endophyte, plants)
qqnorm(log(plants$Dry_shoot_mass))
qqline(log(plants$Dry_shoot_mass))
ggpairs(plants[,1:11])

summary(aov(Leaf_area ~ Endophyte * Virus + Soil_depth, plants))
summary(aov(Dry_shoot_mass ~ Endophyte * Virus + Soil_depth, plants))
summary(aov(Height ~ Endophyte * Virus + Soil_depth, plants))
summary(aov(Dry_root_mass ~ Endophyte * Virus + Soil_depth, plants))
summary(aov(root_shoot_ratio ~ Endophyte * Virus + Soil_depth, plants))
summary(aov(area_height_ratio ~ Endophyte * Virus + Soil_depth, plants))

boxplot(Dry_root_mass ~ Virus + Endophyte, plants)

root_qpcr <- read.csv("VxE_Nov5_qPCR_summary.csv")

root_qpcr

root_qpcr$med_bd <- apply(root_qpcr[,2:4], 1, FUN=median, na.rm=T)

root_qpcr$med_el <- apply(root_qpcr[,5:7], 1, FUN=median, na.rm=T)

root_qpcr$med_fs <- apply(root_qpcr[,8:10], 1, FUN=median, na.rm=T)

root_qpcr$ratio_Bd_El <- 2^(root_qpcr$med_bd-root_qpcr$med_el)*(root_qpcr$Dil_Bd/root_qpcr$Dil_El)

root_qpcr$ratio_Bd_Fs <- 2^(root_qpcr$med_bd-root_qpcr$med_fs)*(root_qpcr$Dil_Bd/root_qpcr$Dil_Fs)

root_qpcr$ratio_El_Fs <- root_qpcr$ratio_Bd_El/root_qpcr$ratio_Bd_Fs

root_qpcr$virus <- as.factor(c(rep(1, 10), rep(0, 10)))

root_qpcr_2 <- read.csv("VxE_Nov5_qPCR_summary_2.csv", nrows = 20)

root_qpcr_2$med_bd_el <- apply(root_qpcr_2[,2:4], 1, FUN=median, na.rm=T)

root_qpcr_2$med_bd_fs <- apply(root_qpcr_2[,5:7], 1, FUN=median, na.rm=T)

root_qpcr_2$med_el <- apply(root_qpcr_2[,8:10], 1, FUN=median, na.rm=T)

root_qpcr_2$med_fs <- apply(root_qpcr_2[,11:13], 1, FUN=median, na.rm=T)

root_qpcr_2$ratio_Bd_El <- 2^(root_qpcr_2$med_bd_el-root_qpcr_2$med_el)*(root_qpcr_2$Dil_Bd_El/root_qpcr_2$Dil_El)

root_qpcr_2$ratio_Bd_Fs <- 2^(root_qpcr_2$med_bd_fs-root_qpcr_2$med_fs)*(root_qpcr_2$Dil_Bd_Fs/root_qpcr_2$Dil_Fs)

root_qpcr_2$ratio_El_Fs <- root_qpcr_2$ratio_Bd_El/root_qpcr_2$ratio_Bd_Fs

root_qpcr_2$virus <- as.factor(c(rep(1, 10), rep(0, 10)))

bp_bd_el <- ggplot(root_qpcr, aes(x=virus, y=ratio_Bd_El, fill=virus, label=Number)) +
  geom_boxplot() + theme_pubr() + scale_y_log10() +
  geom_text(nudge_x = 0.05) + 
  labs(x = "Virus", y = "Relative Bacteria Cells per Plant Cell", tag = "A") +
  theme(legend.position = "none") + geom_jitter(aes(shape=virus)) +
  stat_compare_means(method = "wilcox.test", label.y = -1.3)
bp_bd_el

bp_bd_el_2 <- ggplot(root_qpcr_2, aes(x=virus, y=ratio_Bd_El, color=virus, label=Number)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.4, size = 3) +
  scale_x_discrete(labels=c("-", "+")) +
  theme_pubr() + 
  scale_y_log10() +
  labs(x = "Virus", y = "Relative Bacteria Cells per Plant Cell", tag = "A") +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", label.y = -1.5)
bp_bd_el_2

bp_bd_fs <- ggplot(root_qpcr, aes(x=virus, y=ratio_Bd_Fs, fill=virus, label=Number)) +
  geom_boxplot() + theme_pubr() + scale_y_log10() + 
  geom_text(nudge_x = 0.05) +
  labs(x = "Virus", y = "Relative Fungal Cells per Plant Cell", tag = "B") +
  theme(legend.position = "none") + geom_jitter(aes(shape=virus)) +
  stat_compare_means(method = "wilcox.test", label.y = -1)
bp_bd_fs

bp_bd_fs_2 <- ggplot(root_qpcr_2, aes(x=virus, y=ratio_Bd_Fs, color=virus, label=Number)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.4, size = 3) +
  scale_x_discrete(labels=c("-", "+")) +
  theme_pubr() + scale_y_log10() + 
  labs(x = "Virus", y = "Relative Fungal Cells per Plant Cell", tag = "B") +
  theme(legend.position = "none") + 
  stat_compare_means(method = "wilcox.test", label.y = -2.5)
bp_bd_fs_2

bp_el_fs <- ggplot(root_qpcr, aes(x=virus, y=ratio_El_Fs, fill=virus, label=Number)) +
  geom_boxplot() + theme_pubr() + scale_y_log10() +
  geom_text(nudge_x = 0.05) +
  labs(x = "Virus", y = "Relative Bacteria Cells per Fungal Cell", tag = "C") +
  theme(legend.position = "none") + geom_jitter(aes(shape=virus)) +
  stat_compare_means(method = "wilcox.test", label.y = 2.9)
bp_el_fs
bp_el_fs_2 <- ggplot(root_qpcr_2, aes(x=virus, y=ratio_El_Fs, color=virus, label=Number)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.4, size = 3) +
  scale_x_discrete(labels=c("-", "+")) +
  theme_pubr() + scale_y_log10() +
  labs(x = "Virus", y = "Relative Bacteria Cells per Fungal Cell", tag = "C") +
  theme(legend.position = "none") + 
  stat_compare_means(method = "wilcox.test", label.y = 2.7)
bp_el_fs_2

g <- bp_bd_el | bp_bd_fs | bp_el_fs
g
ggsave("bp_ratios_panel.png", g, dpi=400, width = 10, height = 6, units="in")

g <- bp_bd_el_2 | bp_bd_fs_2 | bp_el_fs_2
g
ggsave("bp_ratios_panel_equal_dilution.png", g, dpi=400, width = 10, height = 6, units="in")

boxplot(log(ratio_Bd_El)~virus, root_qpcr)

boxplot(log(ratio_Bd_Fs)~virus, root_qpcr)

boxplot(log(ratio_El_Fs)~virus, root_qpcr)

model <- t.test(ratio_Bd_El~virus, root_qpcr)
model

model <- t.test(ratio_Bd_Fs~virus, root_qpcr)
model

model <- t.test(ratio_El_Fs~virus, root_qpcr)
model

model <- kruskal.test(ratio_Bd_El~virus, root_qpcr)
model

model <- kruskal.test(ratio_Bd_Fs~virus, root_qpcr)
model

model <- kruskal.test(ratio_El_Fs~virus, root_qpcr)
model

model <- t.test(log(ratio_Bd_El)~virus, root_qpcr_2)
model

model <- lm(log(ratio_Bd_El)~virus, root_qpcr_2)
summary(model)

model <- t.test(log(ratio_Bd_Fs)~virus, root_qpcr_2)
model

model <- t.test(log(ratio_El_Fs)~virus, root_qpcr_2)
model

model <- wilcox.test(ratio_Bd_El~virus, root_qpcr_2)
model

model <- kruskal.test(ratio_Bd_Fs~virus, root_qpcr_2)
model

model <- kruskal.test(ratio_El_Fs~virus, root_qpcr_2)
model

plants_endo <- cbind(root_qpcr_2, plants[root_qpcr_2$Number,])
plants_endo

colnames(plants_endo)

model <- lm(Height ~ ratio_Bd_El + ratio_Bd_Fs + ratio_El_Fs, plants_endo)
summary(model)

model <- lm(Dry_shoot_mass ~ ratio_Bd_El + ratio_Bd_Fs + ratio_El_Fs, plants_endo)
summary(model)

model <- lm(Dry_root_mass ~ ratio_Bd_El + ratio_Bd_Fs + ratio_El_Fs, plants_endo)
summary(model)

model <- lm(Leaf_area ~ ratio_Bd_El + ratio_Bd_Fs + ratio_El_Fs, plants_endo)
summary(model)

model <- lm(root_shoot_ratio ~ ratio_Bd_El + ratio_Bd_Fs + ratio_El_Fs, plants_endo)
summary(model)

model <- lm(area_height_ratio ~ ratio_Bd_El + ratio_Bd_Fs + ratio_El_Fs, plants_endo)
summary(model)

hist(log(plants_endo$ratio_Bd_Fs))
plot(area_height_ratio ~ log(ratio_Bd_Fs), plants_endo, col = virus)
model <- lm(area_height_ratio ~ log(ratio_Bd_Fs), plants_endo)
summary(model)

model <- lm(log(ratio_Bd_El) ~ log(ratio_Bd_Fs), plants_endo)
summary(model)

plot(log(ratio_Bd_El) ~ log(ratio_Bd_Fs), plants_endo[plants_endo$Virus == 0,])
plot(log(ratio_Bd_El) ~ log(ratio_Bd_Fs), plants_endo, col=Virus)
ggplot(plants_endo, aes(x=log(ratio_Bd_Fs), y= log(ratio_Bd_El), col=Virus, label=Number)) +
  geom_point() + geom_text(nudge_x = 0.1)

bartlett.test(ratio_Bd_El ~ virus, plants_endo)
var.test(ratio_Bd_El ~ virus, plants_endo)
bartlett.test(ratio_Bd_Fs ~ virus, plants_endo)
bartlett.test(ratio_El_Fs ~ virus, plants_endo)

plants_endo

