library(dplyr)
library(ggplot2)

viral_load <- read.table('EBV_all_populations_report', sep ='\t')

#histogram for viral load for all population at one time
hist(viral_load$V4,  main="EBV viral load frequency in 1000 Genomes data (peak of viral load)", 
     xlab = "Viral load", xlim = c(0,30), cex.main=1.5, cex.lab=1.3, breaks=500, las = 1)
hist(viral_load$V4,  main="EBV viral load frequency in 1000 Genomes data", 
     xlab = "Viral load", cex.main=1.5, cex.lab=1.3, breaks=200, las = 1)

#histogram of viral load of different population at one graph
density_plot <- ggplot(viral_load, aes(viral_load$V4, fill = viral_load$V1), labels = T)+
  scale_fill_brewer(palette="Set1") + coord_cartesian(xlim = c(0,300)) + 
  geom_density(alpha = 0.5) + ggtitle('The viral load of Epstein-Barr virus in 1000 Genomes data')

print(density_plot + labs(x = 'Viral load', fill = 'Population') + 
        theme(legend.text = element_text(size=15),
              plot.title = element_text(size=20),
              axis.title.x = element_text(size=20),
              axis.title.y = element_text(size=20),
              axis.text.x = element_text(size=20),
              panel.background = element_rect(fill = "white", colour = "grey50")))

#boxpot for determing differences between populations
boxplot(viral_load$V4 ~ viral_load$V1, col = c("#CC6666", "#9999CC", "#66CC99", 'red', 'blue'), ylab = 'viral load', xlab = 'population', ylim = c(0,300))

######################creating a table for GWAS##########################

EBV_GWAS <- viral_load
EBV_GWAS$group <- rep('EBV', 522)
EBV_GWAS <- EBV_GWAS[c('group','V1','V2','V4')]
colnames(EBV_GWAS) <- c('Group', 'Population', 'ID', 'Viral_load')

#creating groups according to the viral load value
EBV_GWAS$Group <-  ifelse(EBV_GWAS$Viral_load < 24, 'EBV_low', 'EBV_high')

#Checking the number of samples in two groups
table(EBV_GWAS$Group)

#exporting the table
write.table(EBV_GWAS, file = 'EBV_GWAS', row.names = F, col.names = T, sep = '\t', quote=FALSE)
