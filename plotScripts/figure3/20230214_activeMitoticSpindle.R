#plot active mitotic spindle promille

library(ggplot2)
library(dplyr)
library(tidyr)

dataDirectory <- "data/spindleAngles/"
plotDirectory <- 'plots/'

allSpindle <- tibble(read.table(file = paste0(dataDirectory, "actmitotSpindleWOP19.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))




dfMean <- df %>% group_by(age) %>% dplyr::summarise(mean = mean(fractionActivePromille)) 
dfSD <- df %>% group_by(age) %>% dplyr::summarise(sd = sd(fractionActivePromille))
df_join <- dfMean %>% left_join(dfSD)

plot <- ggplot(data = df_join, aes(x = age)) +  
  geom_line(aes(y = mean), size = 1, color = "black") + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd), alpha = .2) +
  geom_point(data = df, mapping = aes(x=age, y=fractionActivePromille),size=2, color = "black") +
  theme_bw() + ylab(expression(active~mitotic~nuclei~(per~mille))) + xlab("age (days)") +  
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  scale_x_continuous(limits = c(0,60.5), breaks = c(0,5,7,10, 20, 30, 60)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot, file = paste0(plotDirectory, 'activeMitoticSpindle.svg'), width = 10, height = 5)