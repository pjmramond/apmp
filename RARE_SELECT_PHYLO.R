########################################################################################
# Multivariate selection of rare species
########################################################################################

# Load package
library(vegan)
library(ade4)
library(ggplot2)
library(reshape2)

# Data import and formating
########################################################################################

# otu (swarm2) table
data<-read.table("POHEM/DATA/otu_table.txt", header = TRUE )
rownames(data)<-data$cid

# ID samples
tab<-read.table("POHEM/DATA/id_table.csv", header=TRUE, sep = ";")
tab<-tab[grep("Temoin",tab$ST, invert = TRUE ), ]
tab$Filter<-factor(tab$Filtre, levels=c(20, 10, 3, 0.2),labels = c("Microplankton","Microplankton", "Nanoplankton" ,"Picoplankton"), ordered = TRUE)
tab$Filter<-factor(tab$Filter, levels = c("Microplankton", "Nanoplankton" ,"Picoplankton"), ordered = TRUE)
tab$DATE<-as.POSIXct(tab$DATE,tz="GMT","%d/%m/%Y")
tab<-tab[grep("DA|SEN", tab$CAMP, invert = TRUE),]

# size-fraction subsets 
########################################################################################

# one subdataset for each size-fraction
micro<-data[, colnames(data) %in% tab[tab$Filter == "Microplankton", "TAB"] ];micro<-micro[rowSums(micro) != 0, ]
nano<-data[, colnames(data) %in% tab[tab$Filter == "Nanoplankton", "TAB"] ];nano<-nano[rowSums(nano) != 0, ]
pico<-data[, colnames(data) %in% tab[tab$Filter == "Picoplankton", "TAB"] ];pico<-pico[rowSums(pico) != 0, ]

# rank abundance files for each sub dataset
########################################################################################

# Microplankton otus rank sorted 
######################################@
rmicro<-data.frame(nreads = rowSums(micro), noccur= rowSums(micro !=0), cid = names(rowSums(micro)), SF = "Micro")
# otu rank sorting
rmicro<-rmicro[order(rmicro$nreads, decreasing = TRUE),]
rmicro$rank<-c(1:nrow(rmicro))
rmicro$rela<-rev(rmicro$rank)/nrow(rmicro)*100
# reads ranking and cumulation percentage
rmicro<-rmicro[order(rmicro$rank, decreasing = TRUE),]
rmicro$re_rank<-(rmicro$nreads/sum(rmicro$nreads))*100 
rmicro$cumul<-NA
rmicro$cumul[1]<-rmicro$re_rank[1]
for (i in 2:nrow(rmicro)) { rmicro[i,"cumul"]<-rmicro[i,"re_rank"] + rmicro[i-1,"cumul"] }
# occurence ranking
rmicro<-rmicro[order(rmicro$noccur, decreasing = TRUE),]
rmicro$oc_rank<-c(1:nrow(rmicro))
rmicro$oc_rela<-rmicro$noccur/max(rmicro$noccur)*100

# Nanoplankton otus rank sorted 
######################################@
rnano<-data.frame(nreads = rowSums(nano), noccur= rowSums(nano !=0), cid = names(rowSums(nano)), SF = "nano")
# otu rank sorting
rnano<-rnano[order(rnano$nreads, decreasing = TRUE),]
rnano$rank<-c(1:nrow(rnano))
rnano$rela<-rev(rnano$rank)/nrow(rnano)*100
# reads ranking and cumulation percentage
rnano<-rnano[order(rnano$rank, decreasing = TRUE),]
rnano$re_rank<-(rnano$nreads/sum(rnano$nreads))*100 
rnano$cumul<-NA
rnano$cumul[1]<-rnano$re_rank[1]
for (i in 2:nrow(rnano)) { rnano[i,"cumul"]<-rnano[i,"re_rank"] + rnano[i-1,"cumul"] }
# occurence ranking
rnano<-rnano[order(rnano$noccur, decreasing = TRUE),]
rnano$oc_rank<-c(1:nrow(rnano))
rnano$oc_rela<-rnano$noccur/max(rnano$noccur)*100

# Picoplankton otus rank sorted 
######################################@
rpico<-data.frame(nreads = rowSums(pico), noccur= rowSums(pico !=0), cid = names(rowSums(pico)), SF = "pico")
# otu rank sorting
rpico<-rpico[order(rpico$nreads, decreasing = TRUE),]
rpico$rank<-c(1:nrow(rpico))
rpico$rela<-rev(rpico$rank)/nrow(rpico)*100
# reads ranking and cumulation percentage
rpico<-rpico[order(rpico$rank, decreasing = TRUE),]
rpico$re_rank<-(rpico$nreads/sum(rpico$nreads))*100 
rpico$cumul<-NA
rpico$cumul[1]<-rpico$re_rank[1]
for (i in 2:nrow(rpico)) { rpico[i,"cumul"]<-rpico[i,"re_rank"] + rpico[i-1,"cumul"] }
# occurence ranking
rpico<-rpico[order(rpico$noccur, decreasing = TRUE),]
rpico$oc_rank<-c(1:nrow(rpico))
rpico$oc_rela<-rpico$noccur/max(rpico$noccur)*100

# Rank abundance Gobet selection plots
########################################################################################

load("POHEM/ANALYSIS/RANK_OCC_WODA.RData")
load("POHEM/ANALYSIS/RANK_OCC_WODASEN.RData")
load("POHEM/ANALYSIS/RANK_RANK_WODA.RData")
load("POHEM/ANALYSIS/RANK_RANK_WODASEN.RData")
load("POHEM/ANALYSIS/RANK_READS_WODA.RData")
load("POHEM/ANALYSIS/RANK_READS_WODASEN.RData")

ra<-rbind(ra_micro, ra_nano, ra_pico)
(ra<-melt(ra, measure.vars = c("Mantel's test" ,"Spearman's test")))

ggra<-ggplot(data = ra, aes(x = `% rare reads removed`, y = value ))+
  geom_hline(yintercept = 0.9, color = "grey25", size = 0.75 )+
  geom_vline(xintercept = 50, color = "coral2", size = 1.5)+
  geom_vline(xintercept = 60, color = "coral3", size = 1.5)+
  geom_vline(xintercept = 65, color = "coral4", size = 1.5)+
  geom_point(aes(color = variable, shape = variable), size = 2)+
  geom_line(aes(color = variable ))+
  scale_color_manual(name = "",values = c("seagreen4","dodgerblue4"))+
  scale_shape_manual(name = "",values = c(16,15))+
  scale_y_continuous(name = "Correlation with Original Dataset",limits = c(0,1),breaks = c(0,0.3, 0.6, 0.9) )+
  facet_grid(SF~.)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid")) +
  theme(axis.title = element_text(size = 18, face = "plain"),plot.title = element_text(size = 18),legend.text=element_text(size=18))+
  theme(axis.text = element_text(size = 14))+
  theme(strip.text = element_text(size = 14));ggra

micro_rare<-micro[rownames(micro) %in% rmicro[rmicro$cumul < 30, "cid"], ]; micro_rare<-micro_rare[rowSums(micro_rare) != 0, ]
micro_ab<-micro[rownames(micro) %in% rmicro[rmicro$cumul > 30, "cid"], ]; micro_ab<-micro_ab[rowSums(micro_ab) != 0, ]

nano_rare<-nano[rownames(nano) %in% rnano[rnano$cumul < 30, "cid"], ]; nano_rare<-nano_rare[rowSums(nano_rare) != 0, ]
nano_ab<-nano[rownames(nano) %in% rnano[rnano$cumul > 30, "cid"], ]; nano_ab<-nano_ab[rowSums(nano_ab) != 0, ]

pico_rare<-pico[rownames(pico) %in% rpico[rpico$cumul < 30, "cid"], ]; pico_rare<-pico_rare[rowSums(pico_rare) != 0, ]
pico_ab<-pico[rownames(pico) %in% rpico[rpico$cumul > 30, "cid"], ]; pico_ab<-pico_ab[rowSums(pico_ab) != 0, ]




