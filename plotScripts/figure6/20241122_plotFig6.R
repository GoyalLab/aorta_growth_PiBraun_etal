library(ggplot2)
library(dplyr)
library(tidyr)

dataDirectory <- "ki67Rainbow/first_wave/"
plotDirectory <- "ki67Rainbow/first_wave/plots/"


df_expData_real <- tibble(read.table(file = paste0(dataDirectory, "realDistance_ExponentialData_belowCount10.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
df_expData_negControl <- tibble(read.table(file = paste0(dataDirectory, "negativeControl_ExponentialData_belowCount10.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

df_markDiv_real <- tibble(read.table(file = paste0(dataDirectory, "realDistance_MarkDiv_belowCount10.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
df_markDiv_negControl <- tibble(read.table(file = paste0(dataDirectory, "negativeControl_MarkDiv_belowCount10.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

df_markCl_real <- tibble(read.table(file = paste0(dataDirectory, "realDistance_MarkCl_belowCount10.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
df_markCl_negControl <- tibble(read.table(file = paste0(dataDirectory, "negativeControl_MarkCl_belowCount10.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

#sample neg control and merge
df_expData_negControl_sample <- df_expData_negControl[sample(nrow(df_expData_negControl), 1000), ]
dfEXP = rbind(df_expData_real, df_expData_negControl_sample)

df_markDiv_negControl_sample <- df_markDiv_negControl[sample(nrow(df_markDiv_negControl), 1000), ]
df_markDiv = rbind(df_markDiv_real, df_markDiv_negControl_sample)

df_markCl_negControl_sample <- df_markCl_negControl[sample(nrow(df_markCl_negControl), 1000), ]
df_markCl = rbind(df_markCl_real, df_markCl_negControl_sample)





###### plot real distribution
df_distribution <- tibble(read.table(file = paste0(dataDirectory, "ki67Rainbow_clusterSize_ki67_distribution_expData.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

#0-7; filter cluster sizes which appear less than 10 times
df_distribution_filtered <- df_distribution[df_distribution$age == "0-7" & df_distribution$count >= 10, ]
unique_ticks <- unique(df_distribution_filtered$clusterSize)

plot <- ggplot(df_distribution_filtered, aes(x = clusterSize, y = ki67PerCluster)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.7, fill = "black") + 
  geom_line(color = "black") +
  scale_x_continuous(breaks = unique_ticks) +
  labs(title = "", y = "ki67 active cells per cluster") +
  theme(
    legend.position = "none", # Remove the legend
    axis.line = element_line(colour = "black"), 
    panel.border = element_blank(), 
    panel.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

ggsave(plot, file = paste0(plotDirectory, 'distribution_experimental_0-7_removeClusterSizesBelow10Counts.svg'), width = 5, height = 5)



#1-21; filter cluster sizes which appear less than 10 times

df_distribution_filtered <- df_distribution[df_distribution$age == "1-21" & df_distribution$count >= 10, ]
unique_ticks <- unique(df_distribution_filtered$clusterSize)

plot <- ggplot(df_distribution_filtered, aes(x = clusterSize, y = ki67PerCluster)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.7, fill = "black") + 
  geom_line(color = "black") + 
  scale_x_continuous(breaks = unique_ticks) +
  labs(title = "", y = "ki67 active cells per cluster") +
  theme(
    legend.position = "none", # Remove the legend
    axis.line = element_line(colour = "black"), 
    panel.border = element_blank(), 
    panel.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

ggsave(plot, file = paste0(plotDirectory, 'distribution_experimental_1-21_removeClusterSizesBelow10Counts.svg'), width = 5, height = 5)


## 1-21 all cluster sizes

df_distribution_filtered <- df_distribution[df_distribution$age == "1-21", ]

unique_ticks <- unique(df_distribution_filtered$clusterSize)

plot <- ggplot(df_distribution_filtered, aes(x = clusterSize, y = ki67PerCluster)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.7, fill = "black") + 
  geom_line(color = "black") + 
  scale_x_continuous(breaks = unique_ticks) +
  labs(title = "", y = "ki67 active cells per cluster") +
  theme(
    legend.position = "none", # Remove the legend
    axis.line = element_line(colour = "black"), 
    panel.border = element_blank(), 
    panel.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

ggsave(plot, file = paste0(plotDirectory, 'distribution_experimental_1-21.svg'), width = 5, height = 5)



#############################DISTANCE MEASURESE###########################

###plot exponential
plot_expFunction_Wasserstein = ggplot(data = dfEXP, aes(x=simulation, y=Wasserstein)) +
  geom_boxplot(width = .5, shape =16, outlier.size = 0.1) +
  theme_classic() +
  ylab ("Wasserstein distance from exponential distribution")

ggsave(plot_expFunction_Wasserstein, file = paste0(plotDirectory, 'WassersteinDistanceFromDistribution_exponential.svg'), width = 5, height = 5)


plot_expFunction_JensenShannon = ggplot(data = dfEXP, aes(x=simulation, y=JensenShannon)) +
  geom_boxplot(width = .5, shape =16, outlier.size = 0.1) +
  theme_classic() +
  ylab ("Jensen Shannon divergence from exponential distribution")

ggsave(plot_expFunction_JensenShannon, file = paste0(plotDirectory, 'JensenShannonDivergenceFromDistribution_exponential.svg'), width = 5, height = 5)

plot_expFunction_KL = ggplot(data = dfEXP, aes(x=simulation, y=KL)) +
  geom_boxplot(width = .5, shape =16, outlier.size = 0.1) +
  theme_classic() +
  ylab ("KL divergence from exponential distribution")

ggsave(plot_expFunction_KL, file = paste0(plotDirectory, 'KLdivergenceFromDistribution_exponential.svg'), width = 5, height = 5)


###plot markovian div

plot_markCluster_Wasserstein = ggplot(data = df_markDiv, aes(x=simulation, y=Wasserstein)) +
  geom_boxplot(width = .5, shape =16, outlier.size = 0.1) +
  theme_classic() +
  ylab ("Wasserstein distance from markovian division distribution")

ggsave(plot_markCluster_Wasserstein, file = paste0(plotDirectory, 'WassersteinDistanceFromDistribution_MarkovianDivision.svg'), width = 5, height = 5)

plot_markCluster_JensenShannon = ggplot(data = df_markDiv, aes(x=simulation, y=JensenShannon)) +
  geom_boxplot(width = .5, shape =16, outlier.size = 0.1) +
  theme_classic() +
  ylab ("Jensen Shannon divergence from markovian division distribution")

ggsave(plot_markCluster_JensenShannon, file = paste0(plotDirectory, 'JensenShannonDivergenceFromDistribution_MarkovianDivision.svg'), width = 5, height = 5)

plot_markCluster_KL = ggplot(data = df_markDiv, aes(x=simulation, y=KL)) +
  geom_boxplot(width = .5, shape =16, outlier.size = 0.1) +
  theme_classic() +
  ylab ("KL divergence from markovian division distribution")

ggsave(plot_markCluster_KL, file = paste0(plotDirectory, 'KLdivergenceFromDistribution_MarkovianDivision.svg'), width = 5, height = 5)



### plot mark CLuster
plot_markCluster_Wasserstein = ggplot(data = df_markCl, aes(x=simulation, y=Wasserstein)) +
  geom_boxplot(width = .5, shape =16, outlier.size = 0.1) +
  theme_classic() +
  ylab ("Wasserstein distance from markovian cluster distribution")

ggsave(plot_markCluster_Wasserstein, file = paste0(plotDirectory, 'WassersteinDistanceFromDistribution_MarkovianCluster.svg'), width = 5, height = 5)

plot_markCluster_JensenShannon = ggplot(data = df_markCl, aes(x=simulation, y=JensenShannon)) +
  geom_boxplot(width = .5, shape =16, outlier.size = 0.1) +
  theme_classic() +
  ylab ("Jensen Shannon divergence from markovian cluster distribution")

ggsave(plot_markCluster_JensenShannon, file = paste0(plotDirectory, 'JensenShannonDivergenceFromDistribution_MarkovianCluster.svg'), width = 5, height = 5)

plot_markCluster_KL = ggplot(data = df_markCl, aes(x=simulation, y=KL)) +
  geom_boxplot(width = .5, shape =16, outlier.size = 0.1) +
  theme_classic() +
  ylab ("KL divergence from markovian cluster distribution")

ggsave(plot_markCluster_KL, file = paste0(plotDirectory, 'KLdivergenceFromDistribution_MarkovianCluster.svg'), width = 5, height = 5)

