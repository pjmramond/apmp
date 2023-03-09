##################################################################################
# Script for computing Assembly processes on subsets of ranked 95% OTUs in POHEM
##################################################################################

# Load package
library(sf)
library(broom)
library(vegan)
library(ade4)
library(ggplot2)
library(reshape2)
library(picante)
library(DECIPHER)
library(phangorn)
library(phyloseq)
library(Biostrings)
library(ape)
library(iCAMP)
library(NST)
library(microeco)
library(tidyverse)
library(pspearman)
library(viridis)
library(dplyr)
library(iNEXT)
library(patchwork)
library(limma)
library(ggforce)
library(stringr)
library(FD)
library(mapdata)
library(castor)

# location of the script
rstudioapi::getSourceEditorContext()$path

########################################################################################
# Data import and formating
########################################################################################

# otu (swarm2) table
data<-read.table("/export/lv1/user/pramond/POHEM/DATA/otu_table.txt", header = TRUE )

# ID samples
tab<-read.table("/export/lv1/user/pramond/POHEM/DATA/id_table.csv", header=TRUE, sep = ";")
tab<-tab[grep("Temoin",tab$ST, invert = TRUE ), ]
tab$Filter<-factor(tab$Filtre, levels=c(20, 10, 3, 0.2),labels = c("Microplankton","Microplankton", "Nanoplankton" ,"Picoplankton"), ordered = TRUE)
tab$Filter<-factor(tab$Filter, levels = c("Microplankton", "Nanoplankton" ,"Picoplankton"), ordered = TRUE)
tab$DATE<-as.POSIXct(tab$DATE,tz="GMT","%d/%m/%Y")
tab<-tab[grep("DA|SEN", tab$CAMP, invert = TRUE),]
tab<-tab[tab$TAB %in% colnames(data),]

# env. variables
env<-read.csv("/export/lv1/user/pramond/POHEM/DATA/env_table.csv")
rownames(env)<-env$TAB;env<-na.omit(env[,-1])

# aggregate at 95%
vs95<-read.table("/export/lv1/user/pramond/POHEM/DATA/SEQ_0.95.uc")
vs95<-vs95[vs95$V1 != "C" ,c(1:2, 9)]
colnames(vs95)<-c("type", "ncid", "cid" )
for (i in 1:ncol(vs95)){vs95[,i]<-as.character(vs95[,i])}

# this is a reannotation with latest version of PR2 and mothur, newfound metazoans were discarded
taxo<-read.csv("/export/lv1/user/pramond/POHEM/DATA/SEQ.0_18S_mothur.wang.taxonomy", sep = "\t", header = FALSE)
taxo<-cbind(taxo[,1],stringr::str_split_fixed( taxo[,2], ";", n=8)) # get a column for each taxonomic level
taxo[grep("_unclassified", taxo)]<-NA # replace mothur's "_unclassified" as NAs 
taxo[grep("_X", taxo)]<-NA # replace mothur's "_unclassified" as NAs 
#taxo[grep("Metazoa", taxo)]<-NA # replace mothur's "_unclassified" as NAs 
colnames(taxo)<-c("cid","Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species") # name the columns with PR2 taxonomic levels
taxo<-as.data.frame(taxo) # df
for (j in 1:ncol(taxo)){taxo[,j]<-gsub("\\s*\\([^\\)]+\\)","",as.character(taxo[,j]))} # remove the bootstrap values in parenthesis at the end of each cells in the taxo table
taxo[,-1]<-lapply(taxo[,-1], gsub, pattern=';', replacement='') # remove ";"
taxo$level<-colnames(taxo)[apply(taxo, 1, function(x) length(which(!is.na(x))))] # get taxonomic level at which each sequence is annotated
taxo[is.na(taxo)]<-"Unclassified" # replace NAs by Unclassified

# import otu id data and taxonomy
otu_id<-read.table("/export/lv1/user/pramond/POHEM/DATA/otu_corresp_table.csv", header=TRUE, sep =";")

# import trait table table and formating
trait<-read.table("/export/lv1/user/pramond/POHEM/DATA/trait_table.csv", header=TRUE, sep =";")
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

# import tree build with MAFT and Fasttree
phylo<-read.tree(file = "/export/lv1/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/seqs.tre")

# representative sequence of each otus
seq<-read.table("/export/lv1/user/pramond/POHEM/DATA/sequence.txt", header = TRUE)
seq<-seq[as.character(seq$cid) %in% as.character(data$cid),]
rownames(seq)<-seq$cid

########################################################################################
# Sub-datasets per size fractiion additional information
########################################################################################

# 95% OTUs abundance clustering
data95<-merge(data, vs95, by = "cid")
data95<-aggregate(data95[, grep( "POHEM|RA", colnames(data95))], by = list(data95$ncid),FUN = sum )
data95<-merge(vs95[vs95$type == "S",],data95, by.x = "ncid", by.y = "Group.1")
rownames(data95)<-data95$cid;data95<-data95[, -c(1:3)]
data95<-data95[!rownames(data95) %in% taxo[taxo$Division == "Metazoa", "cid"],]

# Size-fraction subsets
micro95<-data95[, colnames(data95) %in% tab[tab$Filter == "Microplankton", "TAB"] ];micro95<-micro95[rowSums(micro95) != 0, ]
nano95<-data95[, colnames(data95) %in% tab[tab$Filter == "Nanoplankton", "TAB"] ];nano95<-nano95[rowSums(nano95) != 0, ]
pico95<-data95[, colnames(data95) %in% tab[tab$Filter == "Picoplankton", "TAB"] ];pico95<-pico95[rowSums(pico95) != 0, ]

# phyloseq object for each size-fraction
ps_micro<-phyloseq(otu_table(micro95, taxa_are_rows = TRUE), phy_tree(phylo) )
ps_nano<-phyloseq(otu_table(nano95, taxa_are_rows = TRUE), phy_tree(phylo) )
ps_pico<-phyloseq(otu_table(pico95, taxa_are_rows = TRUE), phy_tree(phylo) )

########################################################################################
# Sub-datasets per size fractiion additional information
########################################################################################

# Microplankton otus rank sorted 
rmicro<-data.frame(nreads = rowSums(micro95), noccur= rowSums(micro95 !=0), cid = names(rowSums(micro95)), SF = "Micro")
rmicro<-rmicro[order(rmicro$nreads, decreasing = TRUE),]
rmicro$rank<-c(1:nrow(rmicro))
# reads ranking and cumulation percentage
rmicro<-rmicro[order(rmicro$rank, decreasing = TRUE),]
rmicro$re_rank<-(rmicro$nreads/sum(rmicro$nreads))*100 
rmicro$cumul<-NA
rmicro$cumul[1]<-rmicro$re_rank[1]
for (i in 2:nrow(rmicro)) { rmicro[i,"cumul"]<-rmicro[i,"re_rank"] + rmicro[i-1,"cumul"] }

# Nanoplankton otus rank sorted 
rnano<-data.frame(nreads = rowSums(nano95), noccur= rowSums(nano95 !=0), cid = names(rowSums(nano95)), SF = "Nano")
rnano<-rnano[order(rnano$nreads, decreasing = TRUE),]
rnano$rank<-c(1:nrow(rnano))
# reads ranking and cumulation percentage
rnano<-rnano[order(rnano$rank, decreasing = TRUE),]
rnano$re_rank<-(rnano$nreads/sum(rnano$nreads))*100 
rnano$cumul<-NA
rnano$cumul[1]<-rnano$re_rank[1]
for (i in 2:nrow(rnano)) { rnano[i,"cumul"]<-rnano[i,"re_rank"] + rnano[i-1,"cumul"] }

# Picoplankton otus rank sorted 
rpico<-data.frame(nreads = rowSums(pico95), noccur= rowSums(pico95 !=0), cid = names(rowSums(pico95)), SF = "Pico")
rpico<-rpico[order(rpico$nreads, decreasing = TRUE),]
rpico$rank<-c(1:nrow(rpico))
# reads ranking and cumulation percentage
rpico<-rpico[order(rpico$rank, decreasing = TRUE),]
rpico$re_rank<-(rpico$nreads/sum(rpico$nreads))*100 
rpico$cumul<-NA
rpico$cumul[1]<-rpico$re_rank[1]
for (i in 2:nrow(rpico)) { rpico[i,"cumul"]<-rpico[i,"re_rank"] + rpico[i-1,"cumul"] }

# merge size fraction files
mm<-rbind(rmicro, rnano, rpico)
# create file where OTUs are only present in their main size-fraction
mmr<-as.data.frame(acast(data = mm, cid ~ SF , value.var=c("nreads"), fill = 0))
mmr$Dom_SF<-colnames(mmr)[apply(mmr,1,which.max)]

######################################@######################################
# Rarefaction curves
######################################@######################################

## Diversity estimates with iNEXT
all<-data[,c(colnames(micro95),colnames(nano95),colnames(pico95))];all<-all[rowSums(all) != 0, ]
lsf<-data.frame(All = rowSums(all), 
                Microplankton=rowSums(all[,colnames(micro95)]),
                Nanoplankton=rowSums(all[,colnames(nano95)]),
                Picoplankton=rowSums(all[,colnames(pico95)]) ,
                a2009=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2009,"TAB"  ] ]),
                a2010=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2010,"TAB"  ] ]),
                a2011=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2011,"TAB"  ] ]),
                a2012=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2012,"TAB"  ] ]),
                a2013=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2013,"TAB"  ] ]),
                a2014=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2014,"TAB"  ] ]),
                a2015=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2015,"TAB"  ] ]) )
#div_esti<-iNEXT(x = lsf, datatype="abundance", endpoint = sum(lsf[,"All"])*3)
div_esti<-readRDS("/export/lv6/user/pramond/POHEM/DATA/OTU_95/iNEXT.RData")

# data for plt
df <- fortify(div_esti, type=1)
curves_labels<-aggregate( cbind(x,y) ~ site,  data = df, FUN = max)
curves_labels$type<-ifelse(curves_labels$site %in% c("All", "Microplankton", "Nanoplankton", "Picoplankton"), "SF", "Years" )
curves_labels$site<-gsub("a20","20",curves_labels$site)
df$type<-ifelse(df$site %in% c("All", "Microplankton", "Nanoplankton", "Picoplankton"), "SF", "Years" )
df$site<-gsub("a20","20",df$site)

# plot
library(RColorBrewer)

grac<-ggplot(df, aes(x=x, y=y, colour=site)) + 
  facet_grid(type ~.)+
  geom_line(aes(linetype = method), size = 1)+
  geom_point(data = df[df$method == "observed",], size = 2.5)+
  geom_label(data = curves_labels, aes(label = site, x = Inf, y = y, color = site), hjust = 1, size = 4, fontface = "bold", alpha = 0.75)+
  scale_linetype_manual(values = c("dashed", "solid", "blank"), guide = FALSE)+
  scale_color_manual(guide =FALSE , values = c(brewer.pal(name = "YlOrRd", n = 7),rev(c("#00fff7", "#494fc1", "#fd084a", "coral3")) ) )+
  scale_y_continuous(limits = c(-2500,105000), expand = c(0,0) )+
  scale_x_continuous(limits = c(-2e6,4.9e7), expand = c(0,0) )+
  xlab("# reads")+ylab("# OTUs")+
  #theme(aspect.ratio = 1)+
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16, color = "gray5", hjust = 0.5),
        legend.position="none",
        legend.text=element_text(size=16), 
        legend.title = element_text(size=18, face = "bold"),
        axis.line = element_line(size = 1.25),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"), # get rid of legend panel bg
        axis.ticks = element_line(linewidth = 1));grac

######################################@######################################
# Sliding window
######################################@######################################

# Create Ranges
range<-as.data.frame(rbind(cbind(c(seq(2500,nrow(rmicro), 2500),nrow(rmicro)),"Micro"),
                           cbind(c(seq(2500,nrow(rnano), 2500),nrow(rnano)),"Nano"),
                           cbind(c(seq(2500,nrow(rpico), 2500),nrow(rpico)),"Pico") ))
range$names<-paste(range[,2], range[,1], sep = "_")
colnames(range)<-c("range", "SF", "name")
range$range<-as.numeric(range$range)

# make a list of OTUs within ranges of the sliding window
list_range<-list()
for (i in 1:nrow(range)){
  list_range[[i]]<-subset(mm[mm$SF ==range[i,"SF"],], rank > range$range[i]-2500 & rank <= range$range[i]+2500)$cid
  names(list_range)[[i]]<-paste(range[i,2], range[i,1], sep = "_")
}

# list all samples per subset
list_all=NULL
for (i in names(list_range)) {
  if(range[range$name == i,"SF"] == "Micro" ){SF="Micro"; ps_dsb<-ps_micro}
  if(range[range$name == i,"SF"] == "Nano" ){SF="Nano";ps_dsb<-ps_nano}
  if(range[range$name == i,"SF"] == "Pico" ){SF="Pico"; ps_dsb<-ps_pico}

  # data subset
  sub_dsb<-phyloseq(otu_table(subset(otu_table(ps_dsb), rownames(otu_table(ps_dsb)) %in% list_range[[i]] )), phy_tree(phylo))
  sub_dsb<-prune_samples(colSums(sub_dsb@otu_table) > 0, sub_dsb)
  sub_dsb_env<-prune_samples(sample_names(sub_dsb) %in% rownames(env), sub_dsb)
  list_all$subset[i]<-i
  list_all$nb_samples[i]<-nsamples(sub_dsb)
  list_all$nb_env_samples[i]<-nsamples(sub_dsb_env)
  list_all$nb_otus[i]<-ntaxa(sub_dsb)
  print(i)
}
list_all<-data.frame(nb_samples = list_all$nb_samples, nb_env_samples = list_all$nb_env_samples, nb_otus = list_all$nb_otus)

######################################@######################################
# Taxonomic plots in sliding window
######################################@######################################

# format data of otus per subset
sub=NULL
for (i in names(list_range)){
  c<-data.frame(cid = list_range[[i]], id = i)
  sub<-rbind(sub,c)
}

# get taxonomy info of otus in each subset
sub_taxo<-merge(sub, taxo, by = "cid")
sub_taxo$SF<-gsub("_.*", "",sub_taxo$id)
sub_taxo$subset<-gsub(".*_", "",sub_taxo$id)
sub_taxo$subset<-as.numeric(sub_taxo$subset)

# change value of last subset
last<-unique(sub_taxo[,c("SF","subset")]) %>% group_by(SF) %>% arrange(dplyr::desc(subset)) %>% dplyr::slice(2) 
last$last<-last$subset+2500
last<-merge(last, aggregate(subset ~ SF,data = sub_taxo, FUN = max), by = "SF")
for (i in unique(sub_taxo$SF) ){
  sub_taxo[sub_taxo$SF == i & sub_taxo$subset == last[last$SF==i, 4] , "subset"]<-last[last$SF==i, 3]
}
sub_taxo$ab<-1

# Organize taxonomy by supergroup 
unique(sub_taxo$Division)
labtax<-unique(sub_taxo[order(sub_taxo$Supergroup),3:5]);rownames(labtax)<-NULL
labtax<-labtax[labtax$Division != "Unclassified",]
labtax<-labtax[order(labtax$Supergroup, labtax$Division),];rownames(labtax)<-NULL
labtax$colors<-c(colorRampPalette(c("aquamarine4", "aquamarine1"))(4),
                 colorRampPalette(c("cornsilk4", "cornsilk2"))(3),
                 colorRampPalette(c("darkorange4", "darkorange2"))(3),
                 colorRampPalette(c("palegreen4", "palegreen1"))(4),
                 "violetred",
                 colorRampPalette(c("dodgerblue4", "lightsteelblue3"))(6),
                 colorRampPalette(c("goldenrod4", "goldenrod1"))(3),
                 colorRampPalette(c("slateblue4", "slateblue3"))(2),
                 colorRampPalette(c("red4", "pink1"))(4)
)
labtax<-rbind(labtax, data.frame(Kingdom = "Unclassified", Supergroup = "Unclassified", Division = "Unclassified", colors= "gray25"))
sub_taxo$Division<-factor(sub_taxo$Division, levels = c(labtax$Division))
pal<-labtax$colors;names(pal)<-labtax$Division

# tax plot per subset
ag_sub_taxo<-aggregate(ab ~id + SF + subset + Division , data = sub_taxo, FUN =sum)
ana_tax<-merge(ag_sub_taxo, list_all, by.x = "id", by.y = "row.names")
ana_tax$per_otus<-(ana_tax$ab/ana_tax$nb_otus)*100
gstaxo<-ggplot(data = ag_sub_taxo , aes(x = subset, y = ab, fill = Division))+
  facet_grid(SF~., scale = "free_x", space = "free_x")+
  geom_bar(stat = "identity", width = 2200)+
  scale_y_continuous(expand =c(0,0), limits = c(0, 5500), breaks = seq(0,5000, 1000))+
  scale_x_continuous(expand =c(0,0), limits =c(0,15000), breaks = seq(0,12500, 2500))+
  scale_fill_manual(values = c(pal))+
  xlab("Subsets of Rank Abundance Curves (5000 OTUs)")+
  ylab("# OTUs")+
  ggtitle("Protistan Taxonomy across\nRank Abundance Curves")+
  theme(axis.line.x = element_line(colour = 'gray15', linewidth =1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', linewidth =1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text.x =  element_blank())+
  theme(strip.text.y =  element_text(size = 18, face = "bold", angle =0, hjust = 0))+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0))+
  guides(fill=guide_legend(ncol=2));gstaxo

######################################@#######################################
# Compute ßNTI and RC for subsets
######################################@######################################

# bNTI.big and RC.pc {iCAMP}: Loop to compute ßNTI and RC on subdatasets as phyloseq objects with a phylogenetic tree
list_sub_bnti<-list() # prepare list to receive ßNTI per subdataset
list_sub_rc<-list() # prepare list to receive RC per subdataset
list_sub_bmntd<-list() # prepare list to receive bMNTD per subdataset

# LOOP
for (i in 1:nrow(range) ) {
  
  # dataset selection according to size-fraction
  if(range[i,"SF"] == "Micro"){ps_dsb<-ps_micro}
  if(range[i,"SF"] == "Nano"){ps_dsb<-ps_nano}
  if(range[i,"SF"] == "Pico"){ps_dsb<-ps_pico}
  
  # data subset
  sub_dsb<-phyloseq(otu_table(subset(otu_table(ps_dsb), rownames(otu_table(ps_dsb)) %in% list_range[[i]] )), phy_tree(phylo))
  sub_dsb<-prune_samples(colSums(sub_dsb@otu_table) > 0, sub_dsb)
  
  # make sure the names on the phylogeny are ordered the same as the names in otu table
  match.phylo.otu = match.phylo.data(sub_dsb@phy_tree, sub_dsb@otu_table)
  
  # New command to compute BNTI
  save.wd="/export/lv6/user/pramond/POHEM/DATA/OTU_95/TEMP/"
  pd.big=pdist.big(tree = match.phylo.otu$phy, wd=save.wd, nworker = 40 )
  bNTI=bNTI.big(comm=t(match.phylo.otu$data), pd.desc=pd.big$pd.file,pd.spname=pd.big$tip.label,pd.wd=pd.big$pd.wd,
                spname.check=TRUE, nworker=40, weighted=TRUE, exclude.consp=FALSE,rand=999,
                output.dtail=TRUE, RC=FALSE, trace=FALSE)
  
  # New command to compute Raup-Crick
  RC=RC.pc(comm=t(match.phylo.otu$data), rand = 999,
           nworker = 40, weighted = TRUE,
           sig.index="RC", silent = TRUE)
  
  # save files in list
  list_sub_bnti[[i]]<-bNTI$bNTI
  names(list_sub_bnti)[[i]]<-range[i,3]
  list_sub_bmntd[[i]]<-bNTI$bMNTD
  names(list_sub_bmntd)[[i]]<-range[i,3]
  list_sub_rc[[i]]<-RC
  names(list_sub_rc)[[i]]<-range[i,3]
  
  # remove by-products
  unlink("/export/lv6/user/pramond/POHEM/DATA/OTU_95/TEMP/", recursive = TRUE)
  
  # progress
  print(paste("done with = ", i, "/", nrow(range),sep = ""))
}

# save files
#saveRDS(object = list_sub_rc, file = "/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_SUBSET_RC_2500.RData")
#saveRDS(object = list_sub_bnti, file = "/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_SUBSET_BNTI_2500.RData")
#saveRDS(object = list_sub_bmntd, file = "/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_SUBSET_BMNTD_2500.RData")

# import, format and plot output of assembly processes inference
######################################@######################################

list_sub_rc<-readRDS("/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_SUBSET_RC_2500.RData")
list_sub_bnti<-readRDS("/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_SUBSET_BNTI_2500.RData")
list_sub_bmntd<-readRDS("/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_SUBSET_BMNTD_2500.RData")

list_sub_rc<-readRDS("/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_SUBSET_RC.RData")
list_sub_bnti<-readRDS("/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_SUBSET_BNTI.RData")
list_sub_bmntd<-readRDS("/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_SUBSET_BMNTD.RData")


# Loop Assembly process inference across subsets 
list_sub_ass<-list() # prepare list to receive the proportions of assembly processes per subdataset
for (i in 1:nrow(range) ) {
  
  # extract and merge ßNTI and ßRC
  eco<-merge(melt(list_sub_bnti[[i]]),melt(list_sub_rc[[i]]$index), by = c("Var1", "Var2")); colnames(eco)<-c("sp1", "sp2", "bNTI", "bRCBC")
  
  # infer assembly process per community comparison pairs
  eco$`Dominant Ecological Process`<-ifelse(eco$bNTI >= 2, "Variable Selection", 
                                            ifelse(eco$bNTI <= -2, 'Homogeneous Selection',
                                                   ifelse(abs(eco$bNTI) <= 2 & eco$bRCBC >= 0.95, "Dispersal Limitation",
                                                          ifelse(abs(eco$bNTI) <= 2 & eco$bRCBC <= -0.95, "Homogenizing Dispersal",
                                                                 ifelse(abs(eco$bNTI) <= 2 & abs(eco$bRCBC) <= 0.95, "Undominated Scenario", NA)))))
  
  # save files in list
  list_sub_ass[[i]]<-eco
  names(list_sub_ass)[[i]]<-range[i,3]
  remove(eco)
  
  # progress
  print(paste("done with = ", i, "/", nrow(range),sep = ""))
}

# merge and unlist all datasets
sub_ass<-bind_rows(list_sub_ass, .id = "id") # unlist every dataset into a merged object
sub_ass[is.na(sub_ass$`Dominant Ecological Process`),]
sub_ass<-sub_ass[sub_ass$sp1 != sub_ass$sp2,]

# count assembly processes per subset
sub_ass_count<-sub_ass %>% dplyr::count(id, `Dominant Ecological Process`)
sub_ass_count$SF<-gsub("_.*", "",sub_ass_count$id)
sub_ass_count$subset<-as.numeric(gsub(".*_", "",sub_ass_count$id))
sub_ass_count<-na.omit(sub_ass_count) 
sub_ass_count$`Dominant Ecological Process`<-factor(sub_ass_count$`Dominant Ecological Process`, levels = rev(c("Undominated Scenario","Homogenizing Dispersal", "Dispersal Limitation","Homogeneous Selection", "Variable Selection")), ordered = TRUE)

# compute proportions %
sub_ass_count$tot_sub<-NA
sub_ass_count$perc<-NA
for (i in 1:nrow(sub_ass_count) ){
  sub_ass_count[i, "tot_sub"]<-sum(sub_ass_count[sub_ass_count$id %in% sub_ass_count[i,"id"],"n"])
  sub_ass_count[i, "perc"]<-(sub_ass_count[i, "n"]/sub_ass_count[i, "tot_sub"])*100
}

# change value of last subset
last<-unique(sub_ass_count[,c("SF","subset")]) %>% group_by(SF) %>% arrange(dplyr::desc(subset)) %>% dplyr::slice(2) 
last$last<-last$subset+2500
last<-merge(last, aggregate(subset ~ SF,data = sub_ass_count, FUN = max), by = "SF")
for (i in 1:nrow(last) ){
  sub_ass_count[sub_ass_count$subset == last[i, 4], "subset"]<-last[i, 3]
}

sub_ass_count<-merge(sub_ass_count, list_all, by.x = "id", by.y = "row.names")
levels(sub_ass_count$`Dominant Ecological Process`)

# format assembly processes name and types of rarity
sub_ass_count$`Assembly Processes and Types of Rarity`<-ifelse(sub_ass_count$`Dominant Ecological Process` == "Variable Selection", "Variable Selection\nConditionally Rare",
                                                               ifelse(sub_ass_count$`Dominant Ecological Process` == "Homogeneous Selection", "Homogeneous Selection\nPermanently Rare",
                                                                      ifelse(sub_ass_count$`Dominant Ecological Process` == "Dispersal Limitation", "Dispersal Limitation\nTransiently Rare",
                                                                             ifelse(sub_ass_count$`Dominant Ecological Process` == "Homogenizing Dispersal", "Homogenizing Dispersal\nPermanently Rare (with periodic distribution)",
                                                                                    ifelse(sub_ass_count$`Dominant Ecological Process` == "Undominated Scenario", "Undominated Scenario\nTransiently Rare",NA)))))
sub_ass_count$`Assembly Processes and Types of Rarity`<-factor(sub_ass_count$`Assembly Processes and Types of Rarity`, levels = c("Variable Selection\nConditionally Rare", "Homogeneous Selection\nPermanently Rare", "Dispersal Limitation\nTransiently Rare", "Homogenizing Dispersal\nPermanently Rare (with periodic distribution)","Undominated Scenario\nTransiently Rare"), ordered = TRUE)
ag_sub_ass<-aggregate(perc ~ subset+ id+SF+nb_samples,data = sub_ass_count, FUN = sum)

# plot the assembly processes per subdataset and SF
sub_ass_count[sub_ass_count$`Dominant Ecological Process` == "Homogeneous Selection",] %>% 
  arrange(desc(perc)) %>%
  group_by(SF) %>%
  slice(1)

gsass<-ggplot(data = sub_ass_count, aes(x = subset, y = perc, fill = `Assembly Processes and Types of Rarity`))+
  facet_grid(SF~., scale = "free_x", space = "free")+
  geom_bar(stat = "identity", width = 2200)+
  #geom_line()+
  geom_text(data = ag_sub_ass, aes(x =subset, y = 105 , group =id , label = nb_samples), inherit.aes = FALSE, size = 3.00, vjust =0)+
  #geom_text(aes(label = nb_samples, group = id, y = 1), angle = 0, hjust = 0.5, vjust = 0 , color = "gray15", size= 3.00)+
  scale_y_continuous(expand =c(0,0), limits = c(0, 115), breaks = seq(0,100, 25))+
  scale_x_continuous(expand =c(0,0), limits =c(0,15000), breaks = seq(0,12500, 2500))+
  scale_fill_manual(values = c("#063175","#7aadff", "#a60202", "#ff491c", "#cccccc"))+
  xlab("Subsets of Rank Abundance Curves (1000 OTUs)")+ylab("Frequency (%)")+
  ggtitle("Assembly Processes across\nRank Abundance Curves")+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text =  element_text(size = 18, face = "bold", angle =0, hjust = 0))+
  theme(strip.text.y =  element_text(size = 18, face = "bold", angle =0, hjust = 0))+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gsass

############################################################################################################
# NMDS per subsets of rank abundance curves
############################################################################################################

# Loop PERMANOVA across subsets (only for samples that have env. conditions)
list_sub_perma<-list()
list_sub_perma_bray<-list()
for (i in names(list_sub_bmntd) ) {
  
  # ßMNTD distance
  common<-colnames(list_sub_bmntd[[i]])[colnames(list_sub_bmntd[[i]]) %in% rownames(na.omit(env))]
  comm<-list_sub_bmntd[[i]][common,common]
  
  # Bray curtis distance
  if(range[range$name == i,"SF"] == "Micro" ){SF="Micro"; dsb<-micro95}
  if(range[range$name == i,"SF"] == "Nano" ){SF="Nano";dsb<-nano95}
  if(range[range$name == i,"SF"] == "Pico" ){SF="Pico";dsb<-pico95}
  cbray<- dsb[rownames(dsb) %in% list_range[[i]], common]
  bray<-vegdist(t(cbray), method = "bray")
  
  # Env. conditions
  env_cond<-env[rownames(env) %in% common,]
  
  # permanova
  perma<-as.data.frame(adonis2(comm ~ ., data = env_cond, permutations = 999))
  perma_bray<-as.data.frame(adonis2(bray ~ ., data = env_cond, permutations = 999))
  
  # save files in list
  list_sub_perma[[i]]<-perma
  list_sub_perma_bray[[i]]<-perma_bray
  
  # progress
  print(paste("done with = ", match(i, names(list_sub_bmntd)), "/", length(list_sub_bmntd),sep = ""))
}

# save files
#saveRDS(object = list_sub_perma, file = "/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_PERMANOVA_2500.RData")
#saveRDS(object = list_sub_perma_bray, file = "/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_PERMANOVA_BC_2500.RData")
# read files
list_sub_perma<-readRDS("/export/lv6/user/pramond/POHEM/DATA/OTU_95/FIN_PERMANOVA_2500.RData")

# unlist all datasets
sub_perma<-bind_rows(list_sub_perma, .id = "id") # unlist every dataset into a merged object
sub_perma$env<-gsub("\\..*","",rownames(sub_perma))
sub_perma$SF<-gsub("_.*", "",sub_perma$id)
sub_perma$subset<-as.numeric(gsub(".*_", "",sub_perma$id))

# Format for better plotting
last<-unique(sub_perma[,c("SF","subset")]) %>% group_by(SF) %>% arrange(dplyr::desc(subset)) %>% dplyr::slice(2) 
last$last<-last$subset+2500
last<-merge(last, aggregate(subset ~ SF,data = sub_perma, FUN = max), by = "SF")
for (i in unique(sub_perma$SF) ){
  sub_perma[sub_perma$SF == i & sub_perma$subset == last[last$SF==i, 4] , "subset"]<-last[last$SF==i, 3]
}
sub_perma<-sub_perma[!sub_perma$env %in% c("Total", "Residual"),]
sub_perma$env_cond<-ifelse(sub_perma$`Pr(>F)` < 0.05, as.character(sub_perma$env), "N.S.")
sub_perma$env_cond<-factor(sub_perma$env_cond, 
                           levels = c( "N.S.","NOX", "PO4", "SiOH4", "NH4", "Temp", "Sal", "DEPTH", "Lat", "Long"),
                           labels = c( "N.S.","NOX", "PO4", "SiOH4", "NH4", "Temperature", "Salinity", "Depth", "Latitude", "Longitude"), ordered = TRUE)
sub_perma<-merge(sub_perma, list_all, by.x = "id", by.y = "row.names")
ag_perma<-aggregate(abs(R2) ~ subset + nb_env_samples + SF + id, data =  sub_perma, FUN = sum)

# plot PERMANOVA
pal_p<-c("gray75","#9ef7cf","#71bf9b","#46705d", "#152920", "#f0624f", "#913e33", "#4378e0", "#224282","#0c1830")
names(pal_p)<-levels(sub_perma$env_cond)

aggregate(R2~ id, sub_perma, sum)
gsperma<-ggplot(data = sub_perma, aes(x = subset, y = abs(R2), fill = env_cond))+
  facet_grid(SF~.)+
  geom_bar(stat = "identity", width = 2200)+
  scale_fill_manual(name = "Environmental\nVariable",values = pal_p)+
  scale_y_continuous(expand =c(0,0), limits = c(0, 0.49))+
  scale_x_continuous(expand =c(0,0), limits =c(0,15000), breaks = seq(0,12500, 2500))+
  geom_text(data = ag_perma[ag_perma$subset != 10,], aes(x =subset, y = `abs(R2)`+0.025 , group =id , label = nb_env_samples), inherit.aes = FALSE, size = 3.00)+
  #geom_text(aes(label = nb_env_samples, group = id, y = 0.32), angle = 0, hjust = 0.5, vjust = 0 , color = "gray15", size= 3.00)+
  xlab("Subsets of Rank Abundance Curves (1000 OTUs)")+ylab("PERMANOVA's R2")+
  ggtitle("PERMANOVA across\nRank Abundance Curves")+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text.y =  element_text(size = 18, face = "bold", angle =0, hjust = 0))+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gsperma

library(patchwork)
#gstaxo + gsass + gsperma + plot_layout(ncol = 1) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold'))

p1<-gstaxo+ggtitle("") + plot_annotation(title = 'A: Protistan Taxonomy along Rank Abundance Curves') & theme(title = element_text(size = 16, face = "bold"))
p2<-gsass+ggtitle("") + plot_annotation(title = 'B: Assembly Processes along Rank Abundance Curves') & theme(title = element_text(size = 16, face = "bold"))
p3<-gsperma+ggtitle("") + plot_annotation(title = 'C: PERMANOVA along Rank Abundance Curves') & theme(title = element_text(size = 16, face = "bold"))

ggpubr::ggarrange(p1,p2,p3, ncol = 1, align = "hv")

############################################################################################################
# Signal between traits and abundance
############################################################################################################

# file with information on OTUs traits and abundance
mmt<-merge(mm, traitf, by.x = "cid", by.y = "row.names" )
mmt<-merge(mmt, mmr, by.x = 'cid', by.y = "row.names") # add these two lines to sort OTUs into a single SF
mmt<-mmt[mmt$SF == mmt$Dom_SF,]

# create and format output files and list
tt<-as.data.frame(matrix(nrow = 13, ncol = 3 ))
colnames(tt) = colnames(mmt)[c(2,3,5)]
rownames(tt) = colnames(mmt)[8:20]
tt2<-tt
list_spear_R2<-list()
list_spear_pv<-list()

# loop spearman test over SF and metrics
for (sf in c("Micro", "Nano", "Pico")){
  for (i in colnames(tt))  {
    for (j in rownames(tt) ) {
      tt[j , i] = spearman.test(mmt[mmt$SF %in% sf , j], mmt[mmt$SF %in% sf, i])$estimate
      tt2[j , i] = spearman.test(mmt[mmt$SF %in% sf , j], mmt[mmt$SF %in% sf, i])$p.value
    }}
  list_spear_R2[[sf]]<-tt
  list_spear_pv[[sf]]<-tt2
}
list_spear_R2<-bind_rows(list_spear_R2, .id = "SF");list_spear_R2$traits<-gsub("[.].*","",rownames(list_spear_R2))
list_spear_pv<-bind_rows(list_spear_pv, .id = "SF");list_spear_pv$traits<-gsub("[.].*","",rownames(list_spear_pv))

# format for ggplot
R2<-melt(list_spear_R2, id = c("SF","traits"))
pv<-melt(list_spear_pv, id = c("SF","traits"))
co<-merge(R2,pv, by = c("SF","traits", "variable") );colnames(co)[4:5]<-c("R2","pv")
co$sign<-ifelse(co$R2 >0, "+", "-")
co$traits<-factor(co$traits, levels = rev(colnames(mmt)[8:20]), labels = c("Resting Stage", "Symbiosis", "Ingestion", "Plast_Origin", "Motility","Colony Forming", "Polarity","Symmetry", "Spicule", "Shape", "Cover", "SizeMax", "SizeMin"), ordered = TRUE)
co$variable<-factor(co$variable, levels = c("rank", "nreads", "noccur"), labels = c("Rank", "Abundance", "Occurrence") , ordered = TRUE)
#co<-co[co$variable != "Rank",]

# plot correlation 6x5
summary(abs(co$R2))
gcort<-ggplot(data = co, aes(x = variable, y = traits))+  geom_point(aes(size = abs(R2), alpha = abs(R2),color = sign ))+
  geom_point(data = co[co$pv < 0.001, ], aes(size = abs(R2)), shape = 1 , color = "gray15")+
  facet_grid(.~SF)+
  scale_x_discrete(position = "bottom")+
  scale_alpha(guide = 'none')+
  scale_color_manual(values = c("#871717", "#174687"))+
  scale_size(range = c(0,7))+
  labs(size = bquote(R^2), color = bquote("Sign"~R^2 ))+
  theme_bw()+xlab("")+ylab("Traits")+
  theme(strip.text =  element_text(size = 18, face = "bold", angle =0, hjust = 0.5))+
  theme(strip.background = element_blank())+
  theme(panel.border = element_blank(), panel.grid = element_blank() )+
  theme(axis.ticks = element_blank(), axis.text = element_text(size = 16, hjust = 1))+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 16, face = "bold") )+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1) )+
  theme(axis.title = element_text(size = 18, face = "bold"));gcort

############################################################################################################
# Signal between phylogeny and abundance
############################################################################################################

# empty result files
list_sf_pdrank<-list()
list_sf_pdab<-list()

# Loop to study and save the correlation betwen PD, abundance, and rank, per size-fraction
for (i in unique(mm$SF) ) {
  
  # Empty file to receive output
  res.rank<-data.frame(r2 = NULL,pv = NULL, a = NULL,b = NULL)
  res.ab<-data.frame(r2 = NULL,pv = NULL, a = NULL,b = NULL)
  counts.rank<-matrix(data = 0,nrow = length(seq(0,max(mm$rank), 100)), ncol =  length(seq(0,1, 0.01)) );rownames(counts.rank)<-seq(0,max(mm$rank), 100) ;colnames(counts.rank)<-seq(0,1, 0.01)
  counts.ab<-matrix(data = 0,nrow = length(seq(0,max(mm$nreads), 100)), ncol =  length(seq(0,1, 0.01)) );rownames(counts.ab)<-seq(0,max(mm$nreads), 100) ;colnames(counts.ab)<-seq(0,1, 0.01)
  
  # repeat loop for 999 different subsets of 100 OTUs (1000*1000 = 1M pairwise comparison of OTUs)
  for (j in 1:999){
    
    # make a random selection of 100 OTUs 
    rand_select<-as.character(sample(mm[mm$SF == i, "cid"],size = 100, replace = TRUE))
    rank<-as.data.frame(mm[mm$SF == i & mm$cid %in% rand_select,])
    
    # Compute OTU rank distance
    drank<-dist(rank$rank, method = "euclidean")
    drank<-as.matrix(drank); colnames(drank)<-rank$cid; rownames(drank)<-rank$cid
    
    # Compute OTU abundance distance
    dab<-dist(rank$nreads, method = "euclidean")
    dab<-as.matrix(dab); colnames(dab)<-rank$cid; rownames(dab)<-rank$cid
    
    # compute phylogenetic distance
    seqs<-seq[seq$cid %in% rand_select,]#;names(seqs)<-seq[seq$cid %in% names(drank),"cid"]
    seqs2<-seqs$SEQUENCE;names(seqs2)<-seqs$cid
    alignment <- AlignSeqs(DNAStringSet(seqs2, use.names = TRUE), anchor=NA ,verbose = FALSE)
    dphylo <- dist.ml(phyDat(as(alignment, "matrix"), type="DNA") ) 
    dphylo<-as.matrix(dphylo)[rownames(drank),colnames(drank)]
    
    # correlation results
    res.rank[j,"r2"]<-cor.test(drank, dphylo, method = "pearson")$estimate
    res.rank[j,"pv"]<-cor.test(drank, dphylo, method = "pearson")$p.value
    
    res.ab[j,"r2"]<-cor.test(dab, dphylo, method = "pearson")$estimate
    res.ab[j,"pv"]<-cor.test(dab, dphylo, method = "pearson")$p.value
    
    # merge phylo and trait distance
    mdrank<-reshape2::melt(as.matrix(drank))
    mdab<-reshape2::melt(as.matrix(dab))
    mdphylo<-reshape2::melt(as.matrix(dphylo))
    mp<-merge(mdrank, mdab, by= c("Var1","Var2"), suffixes = c(".rank", ".ab"))
    mp<-merge(mp, mdphylo, by= c("Var1","Var2"));colnames(mp)[ncol(mp)]<-"value.phylo"
    mp<-mp[!mp$Var1 == mp$Var2,]
    
    ## abundance
    # linear model coefficients
    reg<-lm(mp$value.phylo ~mp$value.ab)
    res.ab[j,"a"]<-round(reg$coefficients,2)[[2]]
    res.ab[j,"b"]<-round(reg$coefficients,2)[[1]]
    # gather counts information
    cc<-reshape2::melt(table(round(mp$value.ab,-2),round(mp$value.phylo, 2)))
    cc<-cc[cc$value != 0,]
    cc[cc$Var2 > 1, "Var2"]<-1
    for (z in 1:nrow(cc)){
      counts.ab[as.character(cc[z, "Var1"]) , as.character(cc[z, "Var2"])] <- counts.ab[as.character(cc[z, "Var1"]) , as.character(cc[z, "Var2"])] + cc[z,"value"]
    }
    
    ## rank
    # linear model coefficients
    reg<-lm(mp$value.phylo ~mp$value.rank)
    res.rank[j,"a"]<-round(reg$coefficients,2)[[2]]
    res.rank[j,"b"]<-round(reg$coefficients,2)[[1]]
    # gather counts information
    cc<-reshape2::melt(table(round(mp$value.rank,-2),round(mp$value.phylo, 2)))
    cc<-cc[cc$value != 0,]
    cc[cc$Var2 > 1, "Var2"]<-1
    for (z in 1:nrow(cc)){
      counts.rank[as.character(cc[z, "Var1"]) , as.character(cc[z, "Var2"])] <- counts.rank[as.character(cc[z, "Var1"]) , as.character(cc[z, "Var2"])] + cc[z,"value"]
    }
    # timeline 
    print(paste(i, " / ",j, " / 999",sep = "" ) )
  }
  #save outputs
  list_sf_pdab[[i]]$counts<-counts.ab
  list_sf_pdab[[i]]$res<-res.ab
  list_sf_pdrank[[i]]$counts<-counts.rank
  list_sf_pdrank[[i]]$res<-res.rank
}
saveRDS(list_sf_pdrank, "/export/lv6/user/pramond/POHEM/DATA/OTU_95/RANK_PHYLO.RData")
saveRDS(list_sf_pdab, "/export/lv6/user/pramond/POHEM/DATA/OTU_95/AB_PHYLO.RData")

list_sf_pdrank<-readRDS("/export/lv6/user/pramond/POHEM/DATA/OTU_95/RANK_PHYLO.RData")
list_sf_pdab<-readRDS( "/export/lv6/user/pramond/POHEM/DATA/OTU_95/AB_PHYLO.RData")

# test phylo vs rank data import and format
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
llrank$y<-rep(0.45, 3)
llrank$x<-rep(3000, 3)
llrank$label<-paste("Average R2: ", round(llrank$r2,2), "\np.value: ", round(llrank$pv,2), sep = "")

gran<-ggplot(count_sf_pdrank) +
  aes(x = drank, y = dphylo, z = count, group = -1) +
  stat_summary_hex(size = 0.01, fun = function(x) log10(sum(x)), na.rm = TRUE ) +
  facet_grid(.~SF)+
  geom_abline(data = llrank, aes(slope= a, intercept=b), color = "gray5", size = 1.1 ) +
  geom_text(data = llrank, aes(x =x, y = y , label = label),hjust = 0, inherit.aes = FALSE, size = 4, fontface = "bold")+
  scale_fill_continuous(name = "Nb. of Observations\nin 999 subsets\nof 100 OTUs" , type = "viridis", breaks = c(1,2,3,4,5), labels = c(10^1,10^2, 10^3, 10^4,10^5))+
  xlab("Rank-based\nOTUxOTU Distance")+
  ylab("OTUxOTU\nPhylogenetic Distance")+
  ggtitle("")+
  scale_x_continuous(breaks =c(0,2500,5000,7500,10000))+
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
llab$y<-rep(0.45, 3)
llab$x<-rep(2.25, 3)
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

