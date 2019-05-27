library(dplyr)
library(ggplot2)

adeno <- read.table('adeno_all_populations_report', sep ='\t')

#merging data of all adenoviruses viral load for every sample
adeno_new <- adeno %>% 
  group_by(V1,V2) %>%
  summarise(V4 = sum(V4),V5 = sum(V5))


#histogram of viral load of different population at one graph
density_plot_adeno <- ggplot(adeno_new, aes(adeno_new$V5, fill = adeno_new$V1))+
  scale_fill_brewer(palette="Set1") + coord_cartesian(xlim = c(0,40)) + 
  geom_density(alpha = 0.5) + ggtitle('The viral load of Mastadenoviruses in 1000 Genomes data')

print(density_plot_adeno + labs(x = 'Viral load', fill = 'Population') + 
        theme(legend.text = element_text(size=15),
              plot.title = element_text(size=20),
              axis.title.x = element_text(size=20),
              axis.title.y = element_text(size=20),
              axis.text.x = element_text(size=20),
              panel.background = element_rect(fill = "white", colour = "grey50")))

#boxpot for determing differences between populations
boxplot(adeno_new$V5 ~ adeno_new$V1, col = c("#CC6666", "#9999CC", "#66CC99", 'red', 'blue'), ylab = 'viral load', xlab = 'population', ylim = c(0,300))


###################creating a table for GWAS############################

adeno_GWAS <- adeno_new
adeno_GWAS$group <- rep('adenovirus', 123)
adeno_GWAS <- adeno_GWAS[c('group','V1','V2','V5')]
colnames(adeno_GWAS) <- c('Group', 'Population', 'ID', 'Viral_load')

#Uploading samples without adenoviruses and adding them to GWAS table as control
adeno_free <- read.table('adeno_free_control_samples', sep ='\t')

#adding 0 in viral load for control
adeno_free$Viral_load <- rep(0, 402)

#filter only unique samples
adeno_free <- adeno_free %>% 
  group_by(V1,V2) %>%
  summarise(Viral_load = sum(Viral_load))
adeno_free$Group <- rep('Control', 401)

adeno_free <- adeno_free[c('Group','V1','V2','Viral_load')]
colnames(adeno_free) <- c('Group', 'Population', 'ID', 'Viral_load')

#merging two groups in one data frame
adeno_GWAS <- rbind.data.frame(adeno_GWAS,adeno_free)


#Checking the number of samples in two groups
table(adeno_GWAS$Group)


#exporting the table
write.table(adeno_GWAS, file = 'adeno_GWAS', row.names = F, col.names = T, sep = '\t', quote=FALSE)
