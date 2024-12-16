library(ggplot2)
library(dplyr)
library(tidyr)

dataDirectory <- "data/ki67/"
plotDirectory <- "plots/"


df <- tibble(read.table(file = paste0(dataDirectory, "ki67ActiveFinal.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

dfSummary = df %>% group_by(age) %>% dplyr::summarise(mean = mean(percentageKi67Active),
                                                      sd = sd(percentageKi67Active))


plot <- ggplot(data = dfSummary, aes(x = age, y= mean)) +  
  geom_point(data = df, aes(x=age, y=percentageKi67Active, color = aorta),size=2, shape = 16)+ 
  geom_line(size = 0.5, color = "black")+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd), alpha = .2) +
  theme_classic() +  ylab(expression(percentage~KI67~active)) + xlab("age (days)") +  theme(axis.text=element_text(size=22), axis.title=element_text(size=25)) + 
  scale_x_continuous(limits = c(0,60.5), breaks = c(0,3, 5,7,10, 14, 21, 30, 50, 60)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot, file = paste0(plotDirectory, 'ki67Fraction.svg'), width = 6, height = 3)