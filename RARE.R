########################################################################################
# Multivariate selection of rare species
########################################################################################

# Load package
library(vegan)
library(ade4)
library(ggplot2)
library(reshape2)
library(patchwork)
library(stringr)
#library("energy")

# fade plot background
library(ggplot2) 
library(grid)
library(RColorBrewer)

# function
make_gradient <- function(deg = 45, n = 100, cols = blues9) {
  cols <- colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(
    data = rep(seq(0, 1, length.out = n) * cos(rad), n),
    byrow = TRUE,
    ncol = n) +  
    matrix(data = rep(seq(0, 1, length.out = n) * sin(rad), n),
           byrow = FALSE,
           ncol = n)
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(
    image = mat,
    width = unit(1, "npc"),
    height = unit(1, "npc"), 
    interpolate = TRUE
  )
}

# pick colors and angle
g <- make_gradient(deg = -45, n = 500, cols = c("#f7f5f0","#fae3a2"))

# add this to any plot
annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

# Data import and formating
########################################################################################

# otu (swarm2) table
data<-read.table("POHEM/DATA/otu_table.txt", header = TRUE )
rownames(data)<-data$cid

# this is a reannotation with latest version of PR2 and mothur, newfound metazoans were discarded
taxo<-read.csv("POHEM/DATA/SEQ.0_18S_mothur.wang.taxonomy", sep = "\t", header = FALSE)
taxo<-cbind(taxo[,1],str_split_fixed( taxo[,2], ";", n=8)) # get a column for each taxonomic level
taxo[grep("_unclassified", taxo)]<-NA # replace mothur's "_unclassified" as NAs 
taxo[grep("_X", taxo)]<-NA # replace mothur's "_unclassified" as NAs 
#taxo[grep("Metazoa", taxo)]<-NA # replace mothur's "_unclassified" as NAs 
colnames(taxo)<-c("cid","Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species") # name the columns with PR2 taxonomic levels
taxo<-as.data.frame(taxo) # df
for (j in 1:ncol(taxo)){taxo[,j]<-gsub("\\s*\\([^\\)]+\\)","",as.character(taxo[,j]))} # remove the bootstrap values in parenthesis at the end of each cells in the taxo table
taxo[,-1]<-lapply(taxo[,-1], gsub, pattern=';', replacement='') # remove ";"
taxo$level<-colnames(taxo)[apply(taxo, 1, function(x) length(which(!is.na(x))))] # get taxonomic level at which each sequence is annotated
taxo[is.na(taxo)]<-"Unclassified" # replace NAs by Unclassified
data<-data[!data$cid %in% taxo[taxo$Division == "Metazoa", "cid"],]
taxo<-taxo[taxo$Division != "Metazoa",]

# ID samples
tab<-read.table("POHEM/DATA/id_table.csv", header=TRUE, sep = ";")
tab<-tab[grep("Temoin",tab$ST, invert = TRUE ), ]
tab$Filter<-factor(tab$Filtre, levels=c(20, 10, 3, 0.2),labels = c("Microplankton","Microplankton", "Nanoplankton" ,"Picoplankton"), ordered = TRUE)
tab$Filter<-factor(tab$Filter, levels = c("Microplankton", "Nanoplankton" ,"Picoplankton"), ordered = TRUE)
tab$DATE<-as.POSIXct(tab$DATE,tz="GMT","%d/%m/%Y")
tab<-tab[grep("DA|SEN", tab$CAMP, invert = TRUE),]

# import otu id data and taxonomy
otu_id<-read.table("/export/lv1/user/pramond/POHEM/DATA/otu_corresp_table.csv", header=TRUE, sep =";")

# import trait table table and formating
trait<-read.table("/export/lv1/user/pramond/POHEM/DATA/trait_table.csv", header=TRUE, sep =";")
trait<-merge(otu_id, trait, by.x = "lineage", by.y = "Lineage")
trait_NA<-trait
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

########################################################################################
# size-fraction subsets 
########################################################################################

# one subdataset for each size-fraction
micro<-data[, colnames(data) %in% tab[tab$Filter == "Microplankton", "TAB"] ];micro<-micro[rowSums(micro) != 0, ]
nano<-data[, colnames(data) %in% tab[tab$Filter == "Nanoplankton", "TAB"] ];nano<-nano[rowSums(nano) != 0, ]
pico<-data[, colnames(data) %in% tab[tab$Filter == "Picoplankton", "TAB"] ];pico<-pico[rowSums(pico) != 0, ]

########################################################################################
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
rnano<-data.frame(nreads = rowSums(nano), noccur= rowSums(nano !=0), cid = names(rowSums(nano)), SF = "Nano")
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
rpico<-data.frame(nreads = rowSums(pico), noccur= rowSums(pico !=0), cid = names(rowSums(pico)), SF = "Pico")
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

# merge all size-fraction ranked files
mm<-rbind(rmicro, rnano, rpico)

########################################################################################
# Community comparison of the original dataset to subsets with progressive removal of rare OTUs
###############################################################################################

# bray-curtis distance of the original dataset
dmicro<-vegdist(t(decostand(micro, MARGIN = 2 , method = "total")), method = "bray" )
dnano<-vegdist(t(decostand(nano, MARGIN = 2 , method = "total")), method = "bray" )
dpico<-vegdist(t(decostand(pico, MARGIN = 2 , method = "total")), method = "bray")

# reads abundance cut-off
########################################################################################

# all sf
rar<-as.data.frame(matrix(NA, nrow = 60 , ncol = 5)); colnames(rar)<-c("% rare reads removed", "Mantel's test", "Spearman's test", "nrare OTUs","SF")
rar$`% rare reads removed`<-rep(seq(0,95,5), 3)
rar$SF<-c(rep("Micro", 20), rep("Nano",20), rep("Pico", 20))
for (i in 1:nrow(rar) ) {
  if(rar[i,"SF"] == "Micro"){dsf<-dmicro ; dsb<-micro; rsb<-rmicro}
  if(rar[i,"SF"] == "Nano"){dsf<-dnano ; dsb<-nano; rsb<-rnano}
  if(rar[i,"SF"] == "Pico"){dsf<-dpico; dsb<-pico; rsb<-rpico}
  if (nrow(rsb[rsb$cumul > rar[i,1] ,]) < 2) next
  
  rar[ i , "Mantel's test"]<-mantel.rtest(dsf , br<-vegdist(t(decostand(dsb[rownames(dsb) %in% rsb[rsb$cumul > rar[i,1] ,"cid" ], ], MARGIN = 2,method = "total")), method = "bray"), nrepet = 999)$obs
  rar[i , "Spearman's test"]<-cor(dsf, br, method = "spearman")
  rar[i , "nrare OTUs"]<-nrow(rsb[rsb$cumul > rar[i,1] ,])
  
  print(na.omit(rar));remove(br)
}
saveRDS(object = rar,file = "POHEM/ANALYSIS/RANK/FIN_RANK_READS.RData")
rar<-readRDS("POHEM/ANALYSIS/RANK/FIN_RANK_READS.RData")

# OTUs rank abundance cut-off
############################################

# all sf
raf<-as.data.frame(matrix(NA, nrow = 60 , ncol = 5)); colnames(raf)<-c("% rare reads removed", "Mantel's test", "Spearman's test", "nrare OTUs","SF")
raf$`% rare reads removed`<-rep(seq(0,95,5), 3)
raf$SF<-c(rep("Micro", 20), rep("Nano",20), rep("Pico", 20))
for (i in 1:nrow(raf) ) {
  if(raf[i,"SF"] == "Micro"){dsf<-dmicro ; dsb<-micro; rsb<-rmicro}
  if(raf[i,"SF"] == "Nano"){dsf<-dnano ; dsb<-nano; rsb<-rnano}
  if(raf[i,"SF"] == "Pico"){dsf<-dpico; dsb<-pico; rsb<-rpico}
  if (nrow(rsb[rsb$rela > rar[i,1] ,]) < 2) next
  
  raf[ i , "Mantel's test"]<-mantel.rtest(dsf , br<-vegdist(t(decostand(dsb[rownames(dsb) %in% rsb[rsb$rela > raf[i,1] ,"cid" ], ], MARGIN = 2,method = "total")), method = "bray"), nrepet = 999)$obs
  raf[i , "Spearman's test"]<-cor(dsf, br, method = "spearman")
  raf[i , "nrare OTUs"]<-nrow(rsb[rsb$rela > raf[i,1] ,])
  
  print(na.omit(raf));remove(br)
}
saveRDS(object = raf,file = "POHEM/ANALYSIS/RANK/RANK_OTUS.RData")
raf<-readRDS("POHEM/ANALYSIS/RANK/RANK_OTUS.RData")

# OTUs samples occurence cut-off
############################################

# all sf
rao<-as.data.frame(matrix(NA, nrow = 60 , ncol = 5)); colnames(rao)<-c("% rare reads removed", "Mantel's test", "Spearman's test", "nrare OTUs","SF")
rao$`% rare reads removed`<-rep(seq(0,95,5), 3)
rao$SF<-c(rep("Micro", 20), rep("Nano",20), rep("Pico", 20))
for (i in 1:nrow(rao) ) {
  if(rao[i,"SF"] == "Micro"){dsf<-dmicro ; dsb<-micro; rsb<-rmicro}
  if(rao[i,"SF"] == "Nano"){dsf<-dnano ; dsb<-nano; rsb<-rnano}
  if(rao[i,"SF"] == "Pico"){dsf<-dpico; dsb<-pico; rsb<-rpico}
  if (nrow(rsb[rsb$oc_rela > rar[i,1] ,]) < 2) next
  
  rao[ i , "Mantel's test"]<-mantel.rtest(dsf , br<-vegdist(t(decostand(dsb[rownames(dsb) %in% rsb[rsb$oc_rela > rao[i,1] ,"cid" ], ], MARGIN = 2,method = "total")), method = "bray"), nrepet = 999)$obs
  rao[i , "Spearman's test"]<-cor(dsf, br, method = "spearman")
  rao[i , "nrare OTUs"]<-nrow(rsb[rsb$oc_rela > rao[i,1] ,])
  
  print(na.omit(rao));remove(br)
}
saveRDS(object = rao,file = "POHEM/ANALYSIS/RANK/RANK_OCCURENCE.RData")
rao<-readRDS("POHEM/ANALYSIS/RANK/RANK_OCCURENCE.RData")

##############################################
# plots: 
# ggra = effects of the removal of an increasing number of rare reads, otus, and low occuring OTUs
# ggrank = rank abundance curves and removal
##############################################

ra<-na.omit(reshape2::melt(rar, id.vars = c("% rare reads removed", "SF" , "nrare OTUs")  ))
ra90<-ra[ra$value %in% aggregate(value ~SF ,data = ra[ra$value > 0.90,], FUN = min)$value,]
ra95<-ra[ra$value %in% aggregate(value ~SF ,data = ra[ra$value > 0.95,], FUN = min)$value,]

thres<-data.frame(SF = c("Micro","Nano", "Pico"),
                  int = c(40, 55, 55), 
                  rank = c(nrow(mm[mm$SF == "Micro" & mm$cumul > 40,]),
                           nrow(mm[mm$SF == "Nano" & mm$cumul > 55,]),
                           nrow(mm[mm$SF == "Pico" & mm$cumul > 55,]) ),
                  notus= c(nrow(micro),nrow(nano),nrow(pico)),
                  nreads = c(sum(rowSums(micro)), sum(rowSums(nano)), sum(rowSums(pico))))

ggra<-ggplot(data = ra, aes(x = `% rare reads removed`, y = value ))+
  facet_grid(SF~.)+
  geom_hline(yintercept = 0.9, color = "grey25", size = 0.75 )+
  geom_vline(data = thres, aes(xintercept = int), color = "coral3", size = 1.5)+
  geom_point(aes(color = variable, shape = variable), size = 2)+
  geom_line(aes(color = variable ))+
  geom_rect(data = thres, aes(xmin=-Inf, xmax=int, ymin = -Inf, ymax = Inf, NULL, NULL), alpha=0.1, fill = "gray5")+
  scale_color_manual(name = "",values = c("seagreen4","dodgerblue4"))+
  scale_shape_manual(name = "",values = c(16,15))+
  scale_y_continuous(name = "Correlation with\nOriginal Dataset",limits = c(0,1),breaks = c(0,0.3, 0.6, 0.9) )+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(strip.text.y =  element_text(size = 18, face = "bold", angle =0, hjust = 0))+
  theme(strip.background = element_blank())+
  #theme(strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid")) +
  theme(axis.title = element_text(size = 18, face = "plain"),plot.title = element_text(size = 18),legend.text=element_text(size=18))+
  theme(axis.text = element_text(size = 14))+
  theme(strip.text = element_text(size = 14))

ggrank<-ggplot(data = mm, aes(x = rank, y = log10(nreads) ))+
  geom_point(pch = 1, color = "grey25", alpha = 0.4)+
  #geom_bar(stat = "identity", color = "grey25")+
  facet_grid(SF~.)+
  scale_y_continuous(expand = c(0,0), limits = c(-0.5,7),breaks = c(0,2,4,6), labels = c("0", "100", "10,000", "1,000,000"))+
  scale_x_continuous(expand = c(0,0), limits = c(0,79000) )+
  labs(x = "OTUs Rank Abundance", y = "Number of Reads\n(log10 transformed)")+
  geom_vline(data = thres,aes(xintercept = rank), color = "coral3", size = 1.5)+
  geom_rect(data = thres, aes(xmin= rank, xmax=notus, ymin = -Inf, ymax = Inf, NULL, NULL), alpha=0.1, fill = "gray5")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  #theme(strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid")) +
  theme(axis.title = element_text(size = 18, face = "plain"),plot.title = element_text(size = 18),legend.text=element_text(size=18))+
  theme(axis.text = element_text(size = 14))+
  theme(strip.text = element_text(size = 14))

lab<-data.frame(x=c(10,10,10), y=c(0,0,0), SF = c("Micro", "Nano", "Pico"), ra=c("Abundant","Abundant","Abundant"), ab = c("Rare", "Rare", "Rare"))
ggzoom<-ggplot(data = mm, aes(x = rank, y = log10(nreads)))+
  geom_point(pch = 1, color = "grey25", alpha = 0.4)+
  #geom_bar(stat = "identity", color = "grey25")+
  facet_grid(SF~.)+
  scale_y_continuous(expand = c(0,0), limits = c(-0.5,7),breaks = c(0,2,4,6), labels = c("0", "100", "10,000", "1,000,000"))+
  scale_x_continuous(expand = c(0,0), limits = c(0,105) )+
  labs(x = "OTUs Rank Abundance", y = "Number of Reads\n(log10 transformed)")+
  geom_vline(data = thres,aes(xintercept = rank), color = "coral3", size = 1.5)+
  geom_rect(data = thres, aes(xmin= rank, xmax=Inf, ymin = -Inf, ymax = Inf, NULL, NULL), alpha=0.1, fill = "gray5")+
  geom_text(data = lab, aes(x=x, y=y, label=ra), hjust = 0, vjust = 0, size = 5, angle = 90, fontface = "bold")+
  geom_text(data = lab, aes(x=thres$rank+10, y=y, label=ab), hjust = 0, vjust = 0, size = 5, angle = 90, fontface = "bold")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid")) +
  theme(axis.title = element_text(size = 18, face = "plain"),plot.title = element_text(size = 18),legend.text=element_text(size=18))+
  theme(axis.text = element_text(size = 14))+
  theme(strip.text = element_text(size = 14))

ggzoom+theme(strip.text = element_blank()) + 
  ggrank+theme(strip.text = element_blank(), axis.title.y = element_blank()) + 
  ggra + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold'))

################################
# illustration plot
################################

# explanation of multicola
ggplot(data = mm, aes(x = rank, y = log10(nreads)))+
  geom_point(pch = 1, color = "grey25", alpha = 0.4)+
  #geom_bar(stat = "identity", color = "grey25")+
  facet_grid(SF~.)+
  scale_y_continuous(expand = c(0,0), limits = c(-0.5,7),breaks = c(0,2,4,6), labels = c("0", "100", "10,000", "1,000,000"))+
  scale_x_continuous(expand = c(0,0), limits = c(0,105) )+
  labs(x = "OTUs Rank Abundance", y = "Number of Reads\n(log10 transformed)")+
  geom_vline(data = thres,aes(xintercept = rank), color = "coral3", size = 1.5)+
  geom_rect(data = ra, aes(xmin= `nrare OTUs`, xmax=Inf, ymin = -Inf, ymax = Inf, NULL, NULL), alpha=0.1, fill = "gray55")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid")) +
  theme(axis.title = element_text(size = 18, face = "plain"),plot.title = element_text(size = 18),legend.text=element_text(size=18))+
  theme(axis.text = element_text(size = 14))+
  theme(strip.text = element_text(size = 14))

# trait distribution plots
mmtrait<-merge(mm, traitf,by.x = "cid", by.y = "row.names") 
head(mmtrait)

ggplot(data = mmtrait, aes(x = rank, y = log10(nreads) ))+
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_point(pch = 1, aes(color = Ingestion), alpha = 0.7)+
  scale_color_viridis_d(direction = -1,na.value="gray75")+
  #scale_color_manual(values = pal,guide = "none")+
  #geom_point(pch = 1, color = "gray25",alpha = 0.4)+
  #geom_vline(xintercept = seq(0, 70000, 5000), color = "coral3", size = 0.75, alpha = 0.5)+
  facet_grid(SF~.)+
  scale_y_continuous(expand = c(0,0), limits = c(-0.5,7),breaks = c(0,2,4,6), labels = c("0", "100", "10,000", "1,000,000"))+
  scale_x_continuous(expand = c(0,0), limits = c(-1000,79000) )+
  labs(x = "Rank", y = "Abundance")+
  theme(panel.grid = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text=element_text(size=16),
        legend.key  =element_blank())+
  theme(axis.text = element_text(size = 16))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(strip.text.y = element_text(size = 16, face = "bold", angle = 0),
        strip.background = element_blank())+
  guides(color = guide_legend(override.aes = list(stroke = 2)))

# trait distribution plots
mmplot<-mm[mm$SF == "Nano",]
mmplot$trait_sorted<-cut(mmplot$rank, breaks = c(-Inf,seq(0,70000, 20000),Inf))
mmplot$trait_rand<-sample(unique(mmplot$trait_sorted), size = nrow(mmplot), replace = TRUE)
mmplot$trait_sorted_rand<-mmplot$trait_sorted

for (i in levels(mmplot$trait_sorted)[2:5]){
  lines<-rownames(mmplot[mmplot$trait_sorted %in% i, ])
  spp<-sample(x = lines, size = 10000, replace = FALSE)
  mmplot[spp, "trait_sorted_rand"]<-sample(levels(mmplot$trait_sorted)[2:5], size = 10000, replace = TRUE)
}

pal<-sample(c("#ff6961","#ffb480","#f8f38d","#42d6a4","#08cad1","#59adf6","#9d94ff","#c780e8"), size = 8, replace = TRUE)
pal<-c("coral","royalblue","seagreen","lightblue")

ggplot(data = mmplot, aes(x = rank, y = log10(nreads) ))+
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_point(pch = 1, aes(color = trait_sorted_rand), alpha = 0.4)+
  scale_color_manual(values = pal,guide = "none")+
  #geom_point(pch = 1, color = "gray25",alpha = 0.4)+
  #geom_vline(xintercept = seq(0, 70000, 5000), color = "coral3", size = 0.75, alpha = 0.5)+
  facet_grid(SF~.)+
  scale_y_continuous(expand = c(0,0), limits = c(-0.5,7),breaks = c(0,2,4,6), labels = c("0", "100", "10,000", "1,000,000"))+
  scale_x_continuous(expand = c(0,0), limits = c(-1000,79000) )+
  labs(x = "Rank", y = "Abundance")+
  theme(panel.grid = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(axis.title = element_text(size = 18, face = "plain"),plot.title = element_text(size = 18),legend.text=element_text(size=18))+
  theme(axis.text = element_blank())+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(strip.text = element_blank())

(g2+g1)

################################@################################@
# End
################################@################################@

# plot the occurence frequency curves
ggfreq<-ggplot(data = mm, aes(x = oc_rank, y = log10(noccur)+1 ))+
  geom_bar(stat = "identity", color = "grey25")+
  facet_grid(SF~.)+
  scale_y_continuous(expand = c(0,0), limits = c(0,4),breaks = c(0,1,2,3), labels = c("0", "1", "10", "100"))+
  labs(x = "OTUs Rank Sample Frequency", y = "Number of Samples\n(log10 transformed)")+
  #geom_vline(data = df,aes(xintercept = p35), color = "coral3", size = 1.5)+
  #geom_vline(data = df,aes(xintercept = p55), color = "coral4", size = 1.5)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid")) +
  theme(axis.title = element_text(size = 18, face = "plain"),plot.title = element_text(size = 18),legend.text=element_text(size=18))+
  theme(axis.text = element_text(size = 14))+
  theme(strip.text = element_text(size = 14))



