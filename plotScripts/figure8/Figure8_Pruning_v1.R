
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)
library("ggsci")

home1Directory <- "data/cellCycle/clusterSizes/"
plotDirectory <- "plots/figure8/"


## EduData for bioavailability
p03_EdU <- tibble(read.table(file = paste0(home1Directory, "P03_clusterSizes_noBorder.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

p03_EdU <- p03_EdU %>% mutate(EdUpercent = 100*fraction)
                            
p03_EdU.summary <- p03_EdU %>% group_by(age, duration,clusterSize) %>% summarise(meanEdUPercent = mean(EdUpercent),
                                                                                 sdEdUPercent = sd(EdUpercent))

##singletsOnly
p03_EdUSinglets = p03_EdU %>% filter(clusterSize == 1)
p03_EdUSinglets.summary = p03_EdU.summary %>% filter(clusterSize == 1)

plot1 <- ggplot(data=p03_EdUSinglets, aes(x= duration, y=EdUpercent)) +
  geom_jitter(width = 0.25) +
  geom_line(data = p03_EdUSinglets.summary, aes(x= duration, y = meanEdUPercent)) +
  geom_pointrange(data = p03_EdUSinglets.summary, aes(x= duration, y = meanEdUPercent, ymin = meanEdUPercent-sdEdUPercent, ymax = meanEdUPercent+sdEdUPercent)) +
  theme_classic(base_size = 18) +
  xlab("duration (hours)") +
  ylab ("singlet percentage")

ggsave(plot1, file = paste0(plotDirectory, 'singletons_p03_EdU.svg'), width = 5, height = 4)

####P19
p19_EdU <- tibble(read.table(file = paste0(home1Directory, "P19_allClusterSizes_noBorder.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

p19_EdU <- p19_EdU %>% mutate(EdUpercent = 100*fraction)

p19_EdU.summary <- p19_EdU %>% group_by(age, duration,clusterSize) %>% summarise(meanEdUPercent = mean(EdUpercent),
                                                                                 sdEdUPercent = sd(EdUpercent))

##singletsOnly
p19_EdUSinglets = p19_EdU %>% filter(clusterSize == 1)
p19_EdUSinglets.summary = p19_EdU.summary %>% filter(clusterSize == 1)

plot2 <- ggplot(data=p19_EdUSinglets, aes(x= duration, y=EdUpercent)) +
  geom_jitter(width = 0.25) +
  geom_line(data = p19_EdUSinglets.summary, aes(x= duration, y = meanEdUPercent)) +
  geom_pointrange(data = p19_EdUSinglets.summary, aes(x= duration, y = meanEdUPercent, ymin = meanEdUPercent-sdEdUPercent, ymax = meanEdUPercent+sdEdUPercent)) +
  theme_classic(base_size = 18) +
  xlab("duration (hours)") +
  ylab ("singlet percentage")

ggsave(plot2, file = paste0(plotDirectory, 'singletons_p19_EdU.svg'), width = 5, height = 4)


###################################################################
############# Singletons Cell Size, Dec 20, 2024 #####################
###################################################################

home1Directory <- "data/cellCycle/pruningSmallCells_Area/"
plotDirectory <- "plots/figure8/"


## EduData for bioavailability
sinlgetonArea <- tibble(read.table(file = paste0(home1Directory, "areas_P3-15_pruning.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

sinlgetonArea.summary = sinlgetonArea %>% group_by(singleton) %>% summarise(mean = mean(areaMicro),
                                                                            total_N = length(singleton))


############################################################################################################
############# Taken from legacy Figure7_20231220_normalizedDistance_SingletonPruning.R #####################
############################################################################################################

dataDirectory <- "data/cellCycle/distanceSingletons_pruning/"

plotDirectory = "plots/"

df <- tibble(read.table(file = paste0(dataDirectory,"normalizedDistances_48hr.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

df.summary = df %>% group_by(type, aorta) %>% summarise(sd = sd(normalizedDistance), mean = mean(normalizedDistance))

aorta_mapping <- setNames(seq_along(unique(df.summary$aorta)), unique(df.summary$aorta))
df.summary$aorta_numeric <- aorta_mapping[df.summary$aorta]

merged_df <- inner_join(df, df.summary, by = c("type", "aorta"))
plot <- ggplot(merged_df, aes(x = aorta, y = normalizedDistance, fill = type, color = aorta)) + 
  geom_point(shape = 16, size = 1.5, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  facet_wrap(~ type) +
  scale_color_igv() + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .2, position = position_dodge(.9), color = "black") +
  geom_segment(data = df.summary, aes(x = aorta_numeric - 0.45, xend = aorta_numeric + 0.45, y = mean, yend = mean), color = "black") + 
  theme_classic((base_size = 18))+
  theme(strip.background = element_blank())

ggsave(plot, file = paste0(plotDirectory, 'distanceSingleton48hrs.svg'), width = 15, height = 8)

############################################################################################################
################################## KI67 for OSX CRE #########################################
############################################################################################################

dataDirectory <- "data/diseaseMice/"

ki67_wtMut <- tibble(read.table(file = paste0(dataDirectory, "osxcre_P5672130_ki67_density.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE)) 

ki67_wtMutFirstWave = ki67_wtMut %>% filter(age_combined %in% c("5_6", "7"))

ki67_wtMutSecondWave = ki67_wtMut %>% filter(age_combined %in% c("21", "30"))


plot1 <- ggplot(ki67_wtMutFirstWave, aes(x = age_combined, y = 100*ki67Active, color = treatment)) + 
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 1) +
  ylim (0,26) +
  theme_classic() +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 5,
    shape = 18,
    position = position_dodge(width = 0.5)
  ) +
  # Add error bars (95% confidence interval) for each subgroup
  stat_summary(
    fun.data = calc_mean_sem,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"  # Removed legend
  )
ggsave(plot1, file = paste0(plotDirectory, 'fig8_Ki67FirstWave_SEM.svg'), width = 3.5, height = 5)


plot2 <- ggplot(ki67_wtMutSecondWave, aes(x = age_combined, y = 100*ki67Active, color = treatment)) + 
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 1) +
  ylim (0,20) +
  theme_classic() +
  # Add error bars (95% confidence interval) for each subgroup
  stat_summary(
    fun = mean,
    geom = "point",
    size = 5,
    shape = 18,
    position = position_dodge(width = 0.5)
  ) +
  # Add error bars (95% confidence interval) for each subgroup
  stat_summary(
    fun.data = calc_mean_sem,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"  # Removed legend
  )
ggsave(plot2, file = paste0(plotDirectory, 'fig8_Ki67SecondWave_SEM.svg'), width = 3.5, height = 5)

# Function to calculate mean and confidence interval
calc_mean_sem <- function(x) {
  mean_x <- mean(x)
  sem <- sd(x) / sqrt(length(x))    # Standard Error of the Mean
  data.frame(
    y = mean_x,
    ymin = mean_x - sem,            # Mean minus SEM
    ymax = mean_x + sem             # Mean plus SEM
  )
}
calc_mean_ci <- function(x) {
  mean_x <- mean(x)
  se <- sd(x) / sqrt(length(x))
  ci <- qt(0.975, df = length(x) - 1) * se
  data.frame(
    y = mean_x,
    ymin = mean_x - ci,
    ymax = mean_x + ci
  )
}
calc_mean_sd <- function(x) {
  mean_x <- mean(x)
  sd_x <- sd(x)           # Standard Deviation
  data.frame(
    y = mean_x,
    ymin = mean_x - sd_x, # Mean minus SD
    ymax = mean_x + sd_x  # Mean plus SD
  )
}

############################################################################################################
################################## Density (from KI67) for OSX CRE #########################################
############################################################################################################

dataDirectory <- "data/diseaseMice/"

DensityKi67_wtMut <- tibble(read.table(file = paste0(dataDirectory, "osxcre_P5672130_ki67_density.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE)) 

DensityKi67_wtMutFirstWave = ki67_wtMut %>% filter(age_combined %in% c("5_6", "7"))

ki67_wtMutSecondWave = ki67_wtMut %>% filter(age_combined %in% c("21", "30"))

plot3 <- ggplot(ki67_wtMutFirstWave, aes(x = age_combined, y = 1000*density, color = treatment)) + 
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 1) +
  theme_classic() +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 5,
    shape = 18,
    position = position_dodge(width = 0.5)
  ) +
  # Add error bars (95% confidence interval) for each subgroup
  stat_summary(
    fun.data = calc_mean_sem,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"  # Removed legend
  ) +
  ylim(0, 4)
ggsave(plot3, file = paste0(plotDirectory, 'fig8_DensityKi67FirstWave_SEM.svg'), width = 3.5, height = 5)


plot4 <- ggplot(ki67_wtMutSecondWave, aes(x = age_combined, y = 1000*density, color = treatment)) + 
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 1) +
  theme_classic() +
  # Add error bars (95% confidence interval) for each subgroup
  stat_summary(
    fun = mean,
    geom = "point",
    size = 5,
    shape = 18,
    position = position_dodge(width = 0.5)
  ) +
  # Add error bars (95% confidence interval) for each subgroup
  stat_summary(
    fun.data = calc_mean_sem,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"  # Removed legend
  )+
  ylim(0, 2.25)
ggsave(plot4, file = paste0(plotDirectory, 'fig8_DensityKi67SecondWave_SEM.svg'), width = 3.5, height = 5)

############################################################################################################
################################## EdU for OSX CRE #########################################################
############################################################################################################

dataDirectory <- "data/diseaseMice/"

EdU_wtMut <- tibble(read.table(file = paste0(dataDirectory, "osxcre_P19-24_48_72hrs_EdUpositive_density.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE)) 

EdU_wtMut24h = EdU_wtMut %>% filter(age %in% c("24hrs"))

plot5 <- ggplot(EdU_wtMut24h, aes(x = age, y = 100*percentagePositive, color = treatment)) + 
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 1) +
  ylim (0,5) +
  theme_classic() +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 5,
    shape = 18,
    position = position_dodge(width = 0.5)
  ) +
  # Add error bars (95% confidence interval) for each subgroup
  stat_summary(
    fun.data = calc_mean_sem,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"  # Removed legend
  )
ggsave(plot5, file = paste0(plotDirectory, 'fig7_EdU_24h.svg'), width = 2, height = 5)

plot6 <- ggplot(EdU_wtMut24h, aes(x = age, y = 1000*density, color = treatment)) + 
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 1) +
  ylim (0,3.5) +
  theme_classic() +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 5,
    shape = 18,
    position = position_dodge(width = 0.5)
  ) +
  # Add error bars (95% confidence interval) for each subgroup
  stat_summary(
    fun.data = calc_mean_sem,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"  # Removed legend
  )
ggsave(plot6, file = paste0(plotDirectory, 'fig8_DensityEdU_24h.svg'), width = 2, height = 5)

plot8 <- ggplot(EdU_wtMut, aes(x = age, y = 100*density, color = treatment)) + 
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 1) +
  ylim (0,3.5) +
  theme_classic() +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 5,
    shape = 18,
    position = position_dodge(width = 0.5)
  ) +
  # Add error bars (95% confidence interval) for each subgroup
  stat_summary(
    fun.data = calc_mean_sem,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"  # Removed legend
  )


plot8 <- ggplot(EdU_wtMut, aes(x = age, y = 100*percentagePositive, color = treatment)) + 
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 1) +
  ylim (0,7) +
  theme_classic() +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 5,
    shape = 18,
    position = position_dodge(width = 0.5)
  ) +
  # Add error bars (95% confidence interval) for each subgroup
  stat_summary(
    fun.data = calc_mean_sem,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"  # Removed legend
  )

#################################################################################################################
################################## EdU SINGLETS for OSX CRE #####################################################
#################################################################################################################
dataDirectory <- "data/diseaseMice/"
plotDirectory <- "plots/figure8/"

EdU_singlets <- tibble(read.table(file = paste0(dataDirectory, "clusterSizes_FractionSingletsNEW.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE)) 

EdU_singletsSelect = EdU_singlets %>% select(-resolution, -Size, -Count) %>% filter(duration!= "48hrs")

plot4 <- ggplot(EdU_singlets, aes(x = duration, y = 100*PercentageSinglets, color = experiment, group = experiment)) + 
  theme_classic() +
  # Add lines connecting means
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 0.5,
  ) +
  # Add points for means
  stat_summary(
    fun = mean,
    geom = "point",
    size = 5,
    shape = 18,
  ) +
  # Add error bars
  stat_summary(
    fun.data = calc_mean_sem,
    geom = "errorbar",
    width = 0.2,
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  )
ggsave(plot4, file = paste0(plotDirectory, 'fig8_percentSinglets.svg'), width = 4, height = 5)