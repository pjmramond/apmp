########################################################################################
# Script to compute the phylogenetic distance between all OTUs in POHEM
########################################################################################

# Load package
library(DECIPHER)
library(ggplot2)
library(gridExtra)
library(knitr)
library(dada2)
library(phyloseq)
library(phangorn)
library(reshape2)
library(data.table)
library(picante)
library(cluster)
library(seqinr)
library(hexbin)
library(ggpubr)
# terminal: > module load gsl/2.5
library(energy)
library(ggpubr)
library(dplyr)

# Data import and formating
########################################################################################

# otu (swarm2) table
data<-read.table("POHEM/DATA/otu_table.txt", header = TRUE )
rownames(data)<-data$cid

# representative sequence of each otus
seq<-read.table("POHEM/DATA/sequence.txt", header = TRUE)
seq<-seq[as.character(seq$cid) %in% as.character(data$cid),]
rownames(seq)<-seq$cid

# import otu id data
otu_id<-read.table("POHEM/DATA/otu_corresp_table.csv", header=TRUE, sep =";")

# import trait table table and formating
trait<-read.table("POHEM/DATA/trait_table.csv", header=TRUE, sep =";")
trait<-merge(otu_id, trait, by.x = "lineage", by.y = "Lineage")
trait_NA<-na.omit(trait)

#trait_NA[,"SizeMin"]<-log10(trait_NA[,"SizeMin"])
#trait_NA[,"SizeMax"]<-log10(trait_NA[,"SizeMax"])
trait_NA[,"Cyst_Spore"]<-factor(trait_NA[,"Cyst_Spore"], ordered = TRUE)
trait_NA[,"Spicule"]<-factor(trait_NA[,"Spicule"], ordered = TRUE)
trait_NA[,"Cover"]<-factor(trait_NA[,"Cover"], levels=c("Naked", "Organic", "Siliceous","Calcareous","StrontiumSulphate" ), ordered = TRUE )
trait_NA[,"Shape"]<-factor(trait_NA[,"Shape"], levels=c( "Amoeboid","Round", "Elongated"), ordered = TRUE )
trait_NA[,"Symmetry"]<-factor(trait_NA[,"Symmetry"], levels=c( "Asymetrical", "Spherical", "Bilateral", "Radial"), ordered = TRUE )
trait_NA[,"Polarity"]<-factor(trait_NA[,"Polarity"], levels = c("Heteropolar", "Isopolar"), ordered = TRUE)
trait_NA$Colony<-ifelse(trait_NA$Colony!="None", "Colonial","NonColonial")
trait_NA[,"Colony"]<-factor(trait_NA[,"Colony"], ordered = TRUE )
trait_NA[,"Motility"]<-factor(trait_NA[,"Motility"], levels=c( "Attached", "Floater" , "Gliding" , "Swimmer"  ), ordered = TRUE )
trait_NA[,"Plast_Origin"]<-factor(trait_NA[,"Plast_Origin"], levels=c("None" ,"Kleptoplastidic", "Endosymbiotic" ,"Constitutive"  ), ordered = TRUE )
trait_NA[,"Ingestion"]<-factor(trait_NA[,"Ingestion"],levels=c( "No", "Osmotrophic",  "Saprotrophic", "Intern", "Extern" ),labels = c("No", "Osmotrophic",  "Saprotrophic", "Phagotrophic", "Myzocytotic") ,ordered = TRUE )
trait_NA[,"Symbiontic"]<-factor(trait_NA[,"Symbiontic"], levels=c("No","Commensalist",      "MutualistNonPhotosynthetic", "MutualistPhotosynthetic",  "Parasite" ), ordered = TRUE )
rownames(trait_NA)<-trait_NA$cid
traitf<-trait_NA[,c(8:16, 18:20, 22)]; remove(trait_NA)

# additional data import and format
# ID samples
################################
tab<-read.table("POHEM/DATA/id_table.csv", header=TRUE, sep = ";")
tab<-tab[grep("Temoin",tab$ST, invert = TRUE ), ]
tab$Filter<-factor(tab$Filtre, levels=c(20, 10, 3, 0.2),labels = c("Microplankton","Microplankton", "Nanoplankton" ,"Picoplankton"), ordered = TRUE)
tab$Filter<-factor(tab$Filter, levels = c("Microplankton", "Nanoplankton" ,"Picoplankton"), ordered = TRUE)
tab$DATE<-as.POSIXct(tab$DATE,tz="GMT","%d/%m/%Y")
tab<-tab[grep("DA|SEN", tab$CAMP, invert = TRUE),]

# one subdataset for each size-fraction
################################
micro<-data[, colnames(data) %in% tab[tab$Filter == "Microplankton", "TAB"] ];micro<-micro[rowSums(micro) != 0, ]
nano<-data[, colnames(data) %in% tab[tab$Filter == "Nanoplankton", "TAB"] ];nano<-nano[rowSums(nano) != 0, ]
pico<-data[, colnames(data) %in% tab[tab$Filter == "Picoplankton", "TAB"] ];pico<-pico[rowSums(pico) != 0, ]

# Ranks
################################

# Microplankton otus rank sorted 
rmicro<-data.frame(nreads = rowSums(micro), noccur= rowSums(micro !=0), cid = names(rowSums(micro)), SF = "Micro")
rmicro<-rmicro[order(rmicro$nreads, decreasing = TRUE),]
rmicro$rank<-c(1:nrow(rmicro))
# reads ranking and cumulation percentage
rmicro<-rmicro[order(rmicro$rank, decreasing = TRUE),]
rmicro$re_rank<-(rmicro$nreads/sum(rmicro$nreads))*100 
rmicro$cumul<-NA
rmicro$cumul[1]<-rmicro$re_rank[1]
for (i in 2:nrow(rmicro)) { rmicro[i,"cumul"]<-rmicro[i,"re_rank"] + rmicro[i-1,"cumul"] }

# Nanoplankton otus rank sorted 
rnano<-data.frame(nreads = rowSums(nano), noccur= rowSums(nano !=0), cid = names(rowSums(nano)), SF = "Nano")
rnano<-rnano[order(rnano$nreads, decreasing = TRUE),]
rnano$rank<-c(1:nrow(rnano))
# reads ranking and cumulation percentage
rnano<-rnano[order(rnano$rank, decreasing = TRUE),]
rnano$re_rank<-(rnano$nreads/sum(rnano$nreads))*100 
rnano$cumul<-NA
rnano$cumul[1]<-rnano$re_rank[1]
for (i in 2:nrow(rnano)) { rnano[i,"cumul"]<-rnano[i,"re_rank"] + rnano[i-1,"cumul"] }

# Picoplankton otus rank sorted 
rpico<-data.frame(nreads = rowSums(pico), noccur= rowSums(pico !=0), cid = names(rowSums(pico)), SF = "Pico")
rpico<-rpico[order(rpico$nreads, decreasing = TRUE),]
rpico$rank<-c(1:nrow(rpico))
# reads ranking and cumulation percentage
rpico<-rpico[order(rpico$rank, decreasing = TRUE),]
rpico$re_rank<-(rpico$nreads/sum(rpico$nreads))*100 
rpico$cumul<-NA
rpico$cumul[1]<-rpico$re_rank[1]
for (i in 2:nrow(rpico)) { rpico[i,"cumul"]<-rpico[i,"re_rank"] + rpico[i-1,"cumul"] }

# merge size-fraction ranked OTUs
mm<-rbind(rmicro, rnano, rpico)
mmr<-as.data.frame(acast(data = mm, cid ~ SF , value.var=c("nreads"), fill = 0))
mmr$Dom_SF<-colnames(mmr)[apply(mmr,1,which.max)]

# trait selection
###########@###########@###########@
summary(rownames(traitf) %in% unique(mm$cid))

# heterotrophs
#traitf<-traitf[!traitf$Ingestion %in% c("No", "Osmotrophic"),]
# phototrophs
#traitf<-traitf[traitf$Plast_Origin %in% c("Endosymbiotic", "Constitutive"),]

# The loop
###########@###########@###########

# create empty result table
res<-data.frame(r2 = NULL,pv = NULL,centered_dcor = NULL, nc_dcor = NULL,centered_dcov = NULL, nc_dcov = NULL , a = NULL,b = NULL)

# create count table
counts<-data.frame(dtrait = rep(seq(0,1, 0.01),100) , dphylo = rep(seq(0,1, 0.01), each=100), count = 0)

# repeat loop for 999 different subsets of 1000 OTUs (1000*1000 = 1M pairwise comparison of OTUs)
for (i in 1:999){
  
  # make a random selection of 1000 OTUs that have traits
  rand_select<-as.character(sample(rownames(traitf),size = 100, replace = FALSE))
  
  # Compute trait distance
  dtrait<-daisy(x = traitf[rand_select,], metric = "gower")
  
  # compute phylogenetic distance
  seqs <- seq[rand_select,"SEQUENCE"];names(seqs)<-seq[rand_select,"cid"]
  alignment <- AlignSeqs(DNAStringSet(seqs, use.names = TRUE), anchor=NA ,verbose = FALSE)
  dphylo <- dist.ml(phyDat(as(alignment, "matrix"), type="DNA") ) 

  # correlation results
  res[i,"r2"]<-cor.test(dtrait, dphylo, method = "pearson")$estimate
  res[i,"pv"]<-cor.test(dtrait, dphylo, method = "pearson")$p.value
  
  # test distance correlation and covariance from package "energy"
  # Centered transformation of mat/dist in energy
  res[i,"centered_dcor"]<-dcor(Dcenter(dtrait), Dcenter(dphylo))
  res[i,"nc_dcor"]<-dcor(dtrait, dphylo)
  res[i,"centered_dcov"]<-dcov(Dcenter(dtrait), Dcenter(dphylo))
  res[i,"nc_dcov"]<-dcov(dtrait, dphylo)
  
  # merge phylo and trait distance
  mdtrait<-reshape2::melt(as.matrix(dtrait))
  mdphylo<-reshape2::melt(as.matrix(dphylo))
  mp<-merge(mdtrait, mdphylo,by= c("Var1","Var2"))
  mp<-mp[!mp$Var1 == mp$Var2,]
  
  # linear model coefficients
  reg<-lm(mp$value.y ~mp$value.x)
  res[i,"a"]<-round(reg$coefficients,2)[[2]]
  res[i,"b"]<-round(reg$coefficients,2)[[1]]

  # gather counts information
  cc<-reshape2::melt(table(round(mp[,c(3,4)], 2)))
  for (j in 1:nrow(counts)){
    if ( nrow(cc[cc$value.x == counts[j,"dtrait"] & cc$value.y == counts[j,"dphylo"],]) == 0 ) next
    counts[j, "count"]<-counts[j,"count"] + cc[cc$value.x == counts[j,"dtrait"] & cc$value.y == counts[j,"dphylo"], "value" ]
    }
  
  print(paste(i, "/999",sep = "" ) )
}
#saveRDS(res, "POHEM/ANALYSIS/RES999_100_OTUS_TRAIT_PHYLO.RData")                                                                                                
#saveRDS(counts, "POHEM/ANALYSIS/COUNTS999_100_OTUS_TRAIT_PHYLO.RData")                                                                                           

# plot results from 999 linear models
####################################################################

# all counts data are merged and formatted for plotting
all<-readRDS("POHEM/ANALYSIS/PHYLO_TRAIT/COUNTS999_100_OTUS_TRAIT_PHYLO.RData");all$type<-"All Protists"
het<-readRDS("POHEM/ANALYSIS/PHYLO_TRAIT/HET_COUNTS999_100_OTUS_TRAIT_PHYLO.RData");het$type<-"Heterotrophs"
phot<-readRDS("POHEM/ANALYSIS/PHYLO_TRAIT/PHOT_COUNTS999_100_OTUS_TRAIT_PHYLO.RData");phot$type<-"Phototrophs"
pt<-as.data.frame(rbind(all, het, phot))
pt$type<-factor(pt$type, levels = c("All Protists", "Phototrophs", "Heterotrophs"), ordered = TRUE)
pt[pt$count == 0, "count"]<-NA
pt<-na.omit(pt)

# all coefficient data
mall<-as.data.frame(t(round(apply(na.omit(readRDS("POHEM/ANALYSIS/PHYLO_TRAIT/RES999_100_OTUS_TRAIT_PHYLO.RData")), MARGIN = 2, mean),2)));mall$type<-"All Protists"
mhet<-as.data.frame(t(round(apply(na.omit(readRDS("POHEM/ANALYSIS/PHYLO_TRAIT/HET_RES999_100_OTUS_TRAIT_PHYLO.RData")), MARGIN = 2, mean),2)));mhet$type<-"Heterotrophs"
mphot<-as.data.frame(t(round(apply(na.omit(readRDS("POHEM/ANALYSIS/PHYLO_TRAIT/PHOT_RES999_100_OTUS_TRAIT_PHYLO.RData")), MARGIN = 2, mean),2)));mphot$type<-"Phototrophs"
ll<-as.data.frame(rbind(mall, mhet, mphot))
ll$type<-factor(ll$type, levels = c("All Protists", "Phototrophs", "Heterotrophs"), ordered = TRUE)
ll$y<-rep(0.9375, 3)
ll$x<-rep(0.55, 3)
ll$label<-paste("Average R2: ", ll$r2, "\np.value: ", 1e-8, sep = "")

# plot 7.5 x 12
gpt<-ggplot(pt) +
  aes(x = dtrait, y = dphylo, z = count, group = -1) +
  stat_summary_hex(size = 0.1, fun = function(x) log10(sum(x)), na.rm = TRUE ) +
  geom_abline(data = ll, aes(slope= a, intercept=b), color = "gray5", size = 1.1 ) +
  geom_text(data = ll, aes(x =x, y = y , label = label),hjust = 0, inherit.aes = FALSE, size = 4, fontface = "bold")+
  facet_grid(.~type)+
  scale_fill_continuous(name = "Nb. of Observations\nin 999 subsets\nof 100 OTUs" , type = "viridis", breaks = c(1,2,3,4,5), labels = c(10^1,10^2, 10^3, 10^4,10^5))+
  xlab("Trait-based\nOTUxOTU Distance")+
  ylab("Phylogenetic-based\nOTUxOTU Distance")+
  ggtitle("")+
  scale_x_continuous(limits = c(0,1), expand = c(0,0.05))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0.05))+
  theme(panel.grid.minor = element_line(size = 0.1, color = "gray45"))+
  theme(panel.grid.major = element_line(size = 0.1,color = "gray45"))+
  theme(panel.background = element_rect(color = 'gray15', fill = "#FF000000", size = 3))+
  theme(strip.text =  element_text(size = 18, face = "bold"))+
  theme(strip.background = element_blank())+
  theme(axis.ticks = element_line(color = 'gray15', size = 1))+
  theme(legend.title = element_text(size = 20, face = "bold", color = "gray5", hjust = 0))+
  theme(legend.text = element_text(size = 16, color = "gray5"))+
  theme(legend.justification = "top")+
  theme(plot.title = element_text(size = 18, face = "bold", color = "gray5"))+
  theme(axis.title = element_text(size = 18, face = "bold", color = "gray5", hjust = 0.5))+
  theme(axis.text = element_text(size = 16, color = "gray5", hjust = 0.5))+
  theme(aspect.ratio = 1)+
  theme(panel.spacing = unit(1.5, "lines"));gpt

########################################################################################
# phylogenetic signal of rarity
########################################################################################

# empty result files
list_sf_pdrank<-list()
list_sf_pdab<-list()

# Loop to study and save the correlation betwen PD and abundance rank, per size-fraction
for (i in unique(mm$SF) ) {

  # Empty file to receive output
  res<-data.frame(r2 = NULL,pv = NULL, a = NULL,b = NULL)
  counts<-matrix(data = 0,nrow = length(seq(0,628227, 100)), ncol =  length(seq(0,1, 0.01)) );rownames(counts)<-seq(0,628227, 100) ;colnames(counts)<-seq(0,1, 0.01)
  
  # repeat loop for 999 different subsets of 100 OTUs (1000*1000 = 1M pairwise comparison of OTUs)
  for (j in 1:999){
    
    # make a random selection of 100 OTUs that have traits
    rand_select<-as.character(sample(mm[mm$SF == i, "cid"],size = 100, replace = FALSE))
    
    # Compute OTU rank distance
    rank<-as.data.frame(mm[mm$SF == i & mm$cid %in% rand_select,c("nreads")]); colnames(rank)<-"rank";rownames(rank)<-mm[mm$SF == i & mm$cid %in% rand_select,c("cid")]
    drank<-dist(rank, method = "euclidean")
    
    # compute phylogenetic distance
    seqs<-seq[seq$cid %in% rand_select,"SEQUENCE"];names(seqs)<-seq[seq$cid %in% rand_select,"cid"]
    alignment <- AlignSeqs(DNAStringSet(seqs, use.names = TRUE), anchor=NA ,verbose = FALSE)
    dphylo <- dist.ml(phyDat(as(alignment, "matrix"), type="DNA") ) 
    
    # correlation results
    res[j,"r2"]<-cor.test(drank, dphylo, method = "pearson")$estimate
    res[j,"pv"]<-cor.test(drank, dphylo, method = "pearson")$p.value
    
    # merge phylo and trait distance
    mdrank<-reshape2::melt(as.matrix(drank))
    mdphylo<-reshape2::melt(as.matrix(dphylo))
    mp<-merge(mdrank, mdphylo,by= c("Var1","Var2"))
    mp<-mp[!mp$Var1 == mp$Var2,]
    
    # linear model coefficients
    reg<-lm(mp$value.y ~mp$value.x)
    res[j,"a"]<-round(reg$coefficients,2)[[2]]
    res[j,"b"]<-round(reg$coefficients,2)[[1]]
    
    # gather counts information
    cc<-reshape2::melt(table(round(mp$value.x,-2),round(mp$value.y, 2)))
    cc<-cc[cc$value != 0,]
    for (z in 1:nrow(cc)){
      counts[as.character(cc[z, "Var1"]) , as.character(cc[z, "Var2"])] <- counts[as.character(cc[z, "Var1"]) , as.character(cc[z, "Var2"])] + cc[z,"value"]
    }
    # timeline 
    print(paste(i, " / ",j, " / 999",sep = "" ) )
  }
  #save outputs
  list_sf_pdab[[i]]$counts<-counts
  list_sf_pdab[[i]]$res<-res
  
}
#saveRDS(list_sf_pdrank, "POHEM/ANALYSIS/RANK/RANK_PHYLO.RData")
#saveRDS(list_sf_pdab, "POHEM/ANALYSIS/RANK/AB_PHYLO.RData")

# test phylo vs rank data import and format
list_sf_pdrank<-readRDS("POHEM/ANALYSIS/RANK/RANK_PHYLO.RData")
count_sf_pdrank=NULL;res_rank=NULL
for (i in names(list_sf_pdrank)){
  # counts
  c<-reshape2::melt(list_sf_pdrank[[i]]$counts)
  c$SF<-paste(i,"plankton", sep="");c<-c[c$value >0,]
  count_sf_pdrank=rbind(count_sf_pdrank,c)
  
  # res
  z<-list_sf_pdrank[[i]]$res;z$SF<-paste(i,"plankton", sep="")
  res_rank=rbind(res_rank,z)
}
colnames(count_sf_pdrank)<-c("drank",'dphylo','count','SF')

# results of the linear model
llrank<-aggregate(.~SF,data = res_rank, FUN = function(x) mean(abs(x)) )
llrank$y<-rep(0.4, 3)
llrank$x<-rep(10000, 3)
llrank$label<-paste("Average R2: ", round(llrank$r2,2), "\np.value: ", round(llrank$pv,2), sep = "")

gran<-ggplot(count_sf_pdrank) +
  aes(x = drank, y = dphylo, z = count, group = -1) +
  stat_summary_hex(size = 0.1, fun = function(x) log10(sum(x)), na.rm = TRUE ) +
  facet_grid(.~SF)+
  geom_abline(data = llrank, aes(slope= a, intercept=b), color = "gray5", size = 1.1 ) +
  geom_text(data = llrank, aes(x =x, y = y , label = label),hjust = 0, inherit.aes = FALSE, size = 4, fontface = "bold")+
  scale_fill_continuous(name = "Nb. of Observations\nin 999 subsets\nof 100 OTUs" , type = "viridis", breaks = c(1,2,3,4,5), labels = c(10^1,10^2, 10^3, 10^4,10^5))+
  xlab("Rank-based\nOTUxOTU Distance")+
  ylab("OTUxOTU\nPhylogenetic Distance")+
  ggtitle("")+
  #scale_x_continuous(limits = c(0,1), expand = c(0,0.05))+
  #scale_y_continuous(limits = c(0,1), expand = c(0,0.05))+
  theme(panel.grid.minor = element_line(size = 0.1, color = "gray15"))+
  theme(panel.grid.major = element_line(size = 0.1,color = "gray15"))+
  theme(panel.background = element_rect(color = 'gray15', fill = "#FF000000", size = 3))+
  theme(strip.text =  element_text(size = 18, face = "bold"))+
  theme(strip.background = element_blank())+
  theme(axis.ticks = element_line(color = 'gray15', size = 1))+
  theme(legend.title = element_text(size = 20, face = "bold", color = "gray5", hjust = 0))+
  theme(legend.text = element_text(size = 16, color = "gray5"))+
  theme(plot.title = element_text(size = 18, face = "bold", color = "gray5"))+
  theme(axis.title = element_text(size = 18, face = "bold", color = "gray5", hjust = 0.5))+
  theme(axis.text = element_text(size = 16, color = "gray5", hjust = 0.5))+
  theme(aspect.ratio = 1);gran

# test phylo vs ab data import and format
list_sf_pdab<-readRDS("POHEM/ANALYSIS/RANK/AB_PHYLO.RData")
count_sf_pdab=NULL;res_ab=NULL
for (i in names(list_sf_pdab)){
  #counts
  c<-reshape2::melt(list_sf_pdab[[i]]$counts)
  c$SF<-paste(i,"plankton", sep="");c<-c[c$value >0,]
  count_sf_pdab=rbind(count_sf_pdab,c)
  
  # res 
  z<-list_sf_pdrank[[i]]$res;z$SF<-paste(i,"plankton", sep="")
  res_ab=rbind(res_ab,z)
}
colnames(count_sf_pdab)<-c("dab",'dphylo','count','SF')

# Results from the linear model
llab<-aggregate(.~SF,data = res_ab, FUN = function(x) mean(abs(x)))
llab$y<-rep(0.4, 3)
llab$x<-rep(2.5, 3)
llab$label<-paste("Average R2: ", round(llab$r2,2), "\np.value: ", round(llab$pv,2), sep = "")

# plot
gab<-ggplot(count_sf_pdab) +
  aes(x = log10(dab+1), y = dphylo, z = count, group = -1) +
  stat_summary_hex(size = 0.1, fun = function(x) log10(sum(x)), na.rm = TRUE ) +
  facet_grid(.~SF)+
  geom_abline(data = llab, aes(slope= a, intercept=b), color = "gray5", size = 1.1 ) +
  geom_text(data = llab, aes(x =x, y = y , label = label),hjust = 0, inherit.aes = FALSE, size = 4, fontface = "bold")+
  scale_fill_continuous(name = "Nb. of Observations\nin 999 subsets\nof 100 OTUs" , type = "viridis", breaks = c(1,2,3,4,5), labels = c(10^1,10^2, 10^3, 10^4,10^5))+
  xlab("OTUxOTU\nAbundance Difference")+
  ylab("OTUxOTU\nPhylogenetic Distance")+
  scale_x_continuous(limits = c(0,6.05),breaks = c(0,2,4,6), labels = c("0","1e2", "1e4", "1e6"))+
  ggtitle("")+
  theme(panel.grid.minor = element_line(size = 0.1, color = "gray15"))+
  theme(panel.grid.major = element_line(size = 0.1,color = "gray15"))+
  theme(panel.background = element_rect(color = 'gray15', fill = "#FF000000", size = 3))+
  theme(strip.text =  element_text(size = 18, face = "bold"))+
  theme(strip.background = element_blank())+
  theme(axis.ticks = element_line(color = 'gray15', size = 1))+
  theme(legend.title = element_text(size = 20, face = "bold", color = "gray5", hjust = 0))+
  theme(legend.text = element_text(size = 16, color = "gray5"))+
  theme(plot.title = element_text(size = 18, face = "bold", color = "gray5"))+
  theme(axis.title = element_text(size = 18, face = "bold", color = "gray5", hjust = 0.5))+
  theme(axis.text = element_text(size = 16, color = "gray5", hjust = 0.5))+
  theme(aspect.ratio = 1);gab

# 15 X 10
ggpubr::ggarrange(gran, gab, common.legend = TRUE, 
                  align = "hv", ncol = 1, legend = c("right"),
                  labels = c("A", "B"),font.label = list(face = "bold", size = 18))

