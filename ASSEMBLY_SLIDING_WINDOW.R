##############################################################################
# Script for computing Assembly processes on subsets of ranked OTUs in POHEM
##############################################################################

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

# location of the script
rstudioapi::getSourceEditorContext()$path

########################################################################################
# Data import and formating
########################################################################################

# otu (swarm2) table
data<-read.table("/export/lv1/user/pramond/POHEM/DATA/otu_table.txt", header = TRUE )
data$cid<-as.character(data$cid)
rownames(data)<-data$cid

# import otu id data and taxonomy
otu_id<-read.table("/export/lv1/user/pramond/POHEM/DATA/otu_corresp_table.csv", header=TRUE, sep =";")

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
data<-data[!data$cid %in% taxo[taxo$Division == "Metazoa", "cid"],]
taxo<-taxo[taxo$Division != "Metazoa",]
taxo[grep("Proteobacteria", taxo[,2]),2]<-taxo[grep("Proteobacteria", taxo[,2]),3]

# representative sequence of each otus
seq<-read.table("/export/lv1/user/pramond/POHEM/DATA/sequence.txt", header = TRUE)
seq<-seq[as.character(seq$cid) %in% as.character(data$cid),]
rownames(seq)<-seq$cid
#microseq::writeFasta( data.frame(Header = seq$cid, Sequence = seq$SEQUENCE), out.file = "POHEM/DATA/SEQ.fasta")

# import tree build with MAFT and Fasttree
phylo<-read.tree(file = "/export/lv1/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/seqs.tre")

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

# ID samples
tab<-read.table("/export/lv1/user/pramond/POHEM/DATA/id_table.csv", header=TRUE, sep = ";")
tab<-tab[grep("Temoin",tab$ST, invert = TRUE ), ]
tab$Filter<-factor(tab$Filtre, levels=c(20, 10, 3, 0.2),labels = c("Microplankton","Microplankton", "Nanoplankton" ,"Picoplankton"), ordered = TRUE)
tab$Filter<-factor(tab$Filter, levels = c("Microplankton", "Nanoplankton" ,"Picoplankton"), ordered = TRUE)
tab$DATE<-as.POSIXct(tab$DATE,tz="GMT","%d/%m/%Y")
tab<-tab[grep("DA|SEN", tab$CAMP, invert = TRUE),]
tab<-tab[tab$TAB %in% colnames(data),]
tab$DEPTH

# depth at each location
library(marmap)
atl<-getNOAA.bathy(lon1 = -10, lon2 = 5, lat1 = 40, lat2 = 55)
coord<-unique(na.omit(tab[,c("Long", "Lat")]))
mean(get.depth(atl, x=coord$Long, y=coord$Lat, locator=FALSE)$depth[-nrow(coord)])

# Map of locations
geo<-tab[,c("Long", "Lat","CAMP")]
geo[is.na(geo$Long) ==TRUE,"Long"]<--5.150000
geo[is.na(geo$Lat) ==TRUE,"Lat"]<-48.31667
geo$CAMP<-gsub('[[:digit:]]+', '', geo$CAMP)
geo<-reshape2::melt(table(geo[,c("Long", "Lat","CAMP")]))
geo<-geo[geo$value > 0,]
geo$point<-as.character(geo$CAMP)
geo[geo$CAMP %in% c("PHYTEC", "IPARO") & geo$Lat > 48,"point"]<-"rade"
geo[geo$point %in% c("DYNAPSE","PHYTEC", "IPARO") & geo$Long < -3.75 ,"point"]<-"DY"
geo[geo$point %in% c("PHYTEC", "IPARO") & geo$Long > -3.75 ,"point"]<-"PHY"
geo<-merge(aggregate(cbind(Long,Lat)~point ,data = geo, FUN = mean),aggregate(value~point ,data = geo, FUN = sum), by = "point")

# Cours d'eau, issue de SANDRE.EAUFRANCE.FR
loire <- st_read("/export/lv1/user/pramond/POHEM/DATA/CARTO/CoursEau_04_Loire-Bretagne.shp")
loire<-st_transform(loire, 4326)
loire<-loire[loire$CdEntiteHy == "----0000",]
loire = as.data.frame(st_coordinates(loire))

giro <- st_read("/export/lv1/user/pramond/POHEM/DATA/CARTO/CoursEau_05_Adour-Garonne.shp")
giro<-st_transform(giro, 4326)
garonne<-as.data.frame(st_coordinates(giro[giro$CdEntiteHy %in% c("O---0000"),]))
dordogne<-as.data.frame(st_coordinates(giro[giro$CdEntiteHy %in% c("P---0000"),]))

# importing and formating map details
france.map <- map_data('worldHires','France')
country_shapes <- geom_polygon(data = france.map, aes(x = long, y = lat, group = group),fill = "#CECECE", color = "#CECECE", size = 0.15)
ewbrks <- seq(-6.5,0.5,2)
nsbrks <- seq(45.5,52,2)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(abs(x), "°W"), ifelse(x > 0, paste(abs(x), "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "°S"), ifelse(x > 0, paste(x, "°N"),x))))

#plot map
gc<-ggplot()+ country_shapes +
  geom_path(data = loire, aes(x = X,Y), color = "#596673")+
  geom_path(data = garonne, aes(x = X,Y), color = "#596673")+
  geom_path(data = dordogne, aes(x = X,Y), color = "#596673")+
  geom_jitter(data = geo, aes(x = Long, y = Lat,size = value), color = "coral3", alpha = 0.85)+
  geom_text(data = geo, aes(x = Long, y = Lat, label = value, size = value/25), color = "gray5", fontface = "bold")+
  scale_size(range = c(5,20))+
  #scale_color_manual(values= c("#ebe417", "#99952e", "#4d4b23", "#f22929","#ad2f2f","#4d2323" ))+
  coord_fixed(xlim = c(-7,0), ylim = c(44.5,50.5))+
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0))+
  labs(title = "Spatial Span", y = '', x = '')+
  annotate("text",x=-5, y=47, label = 'atop(italic("Bay of"),italic("        Biscay"))',color = "grey78", parse = TRUE, size = 5)+
  annotate("text",x=-4, y=49.5, label = 'italic("English Channel")',color = "grey78", parse = TRUE, size = 5)+
  annotate("text",x=-6.25, y=48.5, label = 'italic("Iroise \n    Sea")',color = "grey78", parse = TRUE, size = 5)+
  annotate("text",x=-1, y=48, label = 'bold("France")',size = 6,color = "grey18", parse = TRUE)+
  theme(panel.background = element_rect(fill = "#596673"),
        plot.background = element_rect(fill = "transparent", color = "#00000000"),
        panel.grid = element_blank(),
        legend.key  = element_blank(),
        legend.text = element_text(size = 16, color = "gray5"),
        legend.title = element_text(size = 18, color = "gray5", face = "bold"),
        axis.text = element_text(size = 16, color = "gray5"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        title = element_text(size = 16, color = "gray5", face = "bold"))+
  annotate("text",x=-1.5, y=47.5, label = 'italic("Loire River")', color ="#596673",size = 4, parse = TRUE )+
  annotate("text",x=-0.5, y=45.5, label = 'atop(italic("Gironde"),"     River")',color ="#596673" , size = 4, parse = TRUE)+
  guides(size = "none");gc

# seasonality summary
year<-melt(table(format(tab$DATE,"%b-%Y") )); colnames(year)<-c("month","value")
year<-as.data.frame(cbind(str_split_fixed(year$month,"-", 2),year$value)); colnames(year)<-c("month", "year", "value")
year$month<-factor(year$month, levels = substr(month.name,1,3) )
year$Year<-as.numeric(year$year); year$value<-as.numeric(year$value); 

# plot number of samples per month
aggregate(value~Year,year, sum)
gy<-ggplot(year ,aes(x = month, y = value, fill = Year) )+
  geom_bar(stat = "identity", alpha = 0.85)+
  scale_fill_distiller(palette = "YlOrRd", direction = 1)+
  labs(title = "Temporal Span", y = '# of samples', x = '')+
  scale_y_continuous(expand = c(0,0), limits = c(0,325))+
  theme(panel.background = element_rect(fill = "grey88"),
        panel.grid.major.y = element_line(color = "grey78"),
        panel.grid.minor.y = element_line(color = "grey78"),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.key  = element_blank(),
        legend.text = element_text(size = 16, color = "gray5"),
        legend.title = element_text(size = 18, color = "gray5", face = "bold"),
        axis.text = element_text(size = 16, color = "gray5"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1) ,
        title = element_text(size = 16, color = "gray5", face = "bold"));gy
  
# Environmental table
env<-read.csv("/export/lv1/user/pramond/POHEM/DATA/env_table.csv")
rownames(env)<-env$TAB;env<-env[,-1]
env.tab<-na.omit(env[rownames(env) %in% tab$TAB,])
env.tab$DEPTH.y<-ifelse(env.tab$DEPTH.y == "Surface",0,1)
pair<-reshape2::melt(cor(env.tab))
pair[order(abs(pair$value)),]

# format data for plot
env.p<-melt(env.tab[,4:ncol(env.tab)])
env.p$variable<-factor(env.p$variable, 
                           levels = rev(c("NOX", "PO4", "SiOH4", "NH4", "Temp", "Sal", "DEPTH", "Lat", "Long")),
                           labels = rev(c( "NOX", "PO4", "SiOH4", "NH4", "Temperature", "Salinity", "Depth", "Latitude", "Longitude")), ordered = TRUE)
pal_p<-c("#9ef7cf","#71bf9b","#46705d", "#152920", "#f0624f", "#913e33", "#4378e0", "#224282","#0c1830")
names(pal_p)<-c( "NOX", "PO4", "SiOH4", "NH4", "Temperature", "Salinity", "Depth", "Latitude", "Longitude")
env.p$type<-ifelse(env.p$variable %in% c( "NOX", "PO4", "SiOH4", "NH4"), "Chemical",ifelse(env.p$variable %in% c("Temperature", "Salinity"), "Physical","Distance"))
env.p$type<-factor(env.p$type, levels = c("Chemical", "Physical","Distance"), ordered = TRUE)

# plot env. conditions
ge<-ggplot(env.p ,aes(x =variable , y = value, color = variable) )+
  facet_wrap(~type,scales ="free",  drop=TRUE)+
  geom_jitter(alpha = 0.75)+
  scale_fill_manual(values = rev(pal_p))+
  scale_color_manual(values = rev(pal_p))+
  labs(title = "Conditions", y = 'Value', x = '')+
  theme(panel.background = element_rect(fill = "grey88"),
        panel.grid.major.y = element_line(color = "grey78"),
        panel.grid.minor.y = element_line(color = "grey78"),
        strip.background = element_blank(),
        strip.text= element_text(size = 16, color = "gray5"),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.key  = element_blank(),
        legend.text = element_text(size = 16, color = "gray5"),
        legend.title = element_text(size = 18, color = "gray5", face = "bold"),
        axis.text = element_text(size = 16, color = "gray5"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1) ,
        title = element_text(size = 16, color = "gray5", face = "bold"))+
  guides(fill = "none", color = "none");ge

# 7.5 X 18
gc+gy+ge + plot_layout(ncol = 1) +  plot_annotation(tag_levels = 'A')
gc+gy+ge + plot_layout(nrow = 1) +  plot_annotation(tag_levels = 'A')

########################################################################################
# size-fraction subsets 
########################################################################################

# one subdataset for each size-fraction
micro<-data[, colnames(data) %in% tab[tab$Filter == "Microplankton", "TAB"] ];micro<-micro[rowSums(micro) != 0, ]
nano<-data[, colnames(data) %in% tab[tab$Filter == "Nanoplankton", "TAB"] ];nano<-nano[rowSums(nano) != 0, ]
pico<-data[, colnames(data) %in% tab[tab$Filter == "Picoplankton", "TAB"] ];pico<-pico[rowSums(pico) != 0, ]

# info on the data
#sum(rowSums(micro[rownames(micro) %in% rownames(traitf),]))+sum(rowSums(nano[rownames(nano) %in% rownames(traitf),]))+sum(rowSums(pico[rownames(pico) %in% rownames(traitf),]))
#sum(colSums(micro))+sum(colSums(nano))+sum(colSums(pico))
# all<-data[, c(colnames(micro), colnames(nano), colnames(pico))];all<-all[rowSums(all) != 0, ]
# nrow(all[rownames(all) %in% rownames(traitf),])
# sum(rowSums(all[rownames(all) %in% rownames(traitf),]))
# one subdataset for each size-fraction, transformed in proportions
#micro<-decostand(data[, colnames(data) %in% tab[tab$Filter == "Microplankton", "TAB"] ], MARGIN = 2, method = "total" );micro<-micro[rowSums(micro) != 0, ]
#nano<-decostand(data[, colnames(data) %in% tab[tab$Filter == "Nanoplankton", "TAB"] ], MARGIN = 2, method = "total" );nano<-nano[rowSums(nano) != 0, ]
#pico<-decostand(data[, colnames(data) %in% tab[tab$Filter == "Picoplankton", "TAB"] ], MARGIN = 2, method = "total" );pico<-pico[rowSums(pico) != 0, ]

# trait.data<-trait[trait$cid %in% unique(c(rownames(micro), rownames(nano), rownames(pico))),]
# trait.data<-as.data.frame(cbind(stringr::str_split_fixed( trait.data$lineage, "[|]", n=max(str_count(trait.data$lineage, "[|]"))), trait.data))
# trait.data_NA<-na.omit(trait.data)
# trait.data<-merge(trait.data, taxo, by = "cid")
# trait.data_NA<-merge(trait.data_NA, taxo, by = "cid")
# ntaxoref<-melt(table(trait.data$lineage, trait.data$`3`))
# ntaxoref<-ntaxoref[ntaxoref$value>0,]
# unc<-names(table(trait.data$`3`))[!names(table(trait.data$`3`)) %in% names(table(trait.data_NA$`3`))]
# trait.data$`3`<-ifelse(trait.data$`3` %in% unc,"Unclassified",  trait.data$`3`)
# 
# taxo.ref<-data.frame(Division=names(table(trait.data$`3`)), 
#   notus = table(trait.data$`3`),
#   ntaxoref=table(ntaxoref$Var2),
#   avg.ntaxoref=aggregate(value ~ Var2,data = ntaxoref, FUN = mean) );
# taxo.ref<-taxo.ref[, grep("Var", colnames(taxo.ref), invert = TRUE)]
# View(data.frame(table(trait.data_NA$`3`)))

# phyloseq object for each size-fraction
ps_micro<-phyloseq(otu_table(micro, taxa_are_rows = TRUE), phy_tree(phylo) )
ps_nano<-phyloseq(otu_table(nano, taxa_are_rows = TRUE), phy_tree(phylo) )
ps_pico<-phyloseq(otu_table(pico, taxa_are_rows = TRUE), phy_tree(phylo) )

######################################@######################################@
# Size fraction otus rank sorted 
######################################@######################################@

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

# merge size fraction files
mm<-rbind(rmicro, rnano, rpico)
# create file where OTUs are only present in their main size-fraction
mmr<-as.data.frame(acast(data = mm, cid ~ SF , value.var=c("nreads"), fill = 0))
mmr$Dom_SF<-colnames(mmr)[apply(mmr,1,which.max)]

#summary(rownames(mmr) %in% rownames(traitf[traitf$Plast_Origin %in% c("Endosymbiotic", "Constitutive"),]))
########################################################################################
# Information on the amplicon sequencing dataset and plots
########################################################################################

# plot the relation between occurence and abundance
ocab<-ggplot(data = mm, aes(x = noccur, y = log10(nreads), color = SF))+
  scale_y_continuous(expand = c(0,0), limits = c(-0.25,6.25),breaks = c(0,2,4,6), labels = c("0", "100", "10,000", "1,000,000"))+
  scale_x_continuous(expand = c(0,0), limits = c(-12.5, 350))+
  geom_point(aes(alpha = SF) , shape = 1, size = 1.5 ,stroke = 1.5)+
  scale_alpha_manual(values = c(1,0.5,0.25), guide = FALSE)+
  scale_color_manual(name ="" , values = rev(c("#00fff7", "#494fc1", "#fd084a")))+
  labs(x = "OTU Occurrence\n(# samples)", y = "OTU Abundance\n(Log10 # reads)")+
  theme_bw()+
  theme(panel.grid.minor = element_line(size = 0.1, color = "gray15"))+
  theme(panel.grid.major = element_line(size = 0.1,color = "gray15"))+
  theme(panel.background = element_rect(color = 'gray15', fill = "#FF000000", size = 2))+
  theme(axis.ticks = element_blank())+
  #theme(legend.position = c(0.8, 0.2))+
  theme(legend.justification = "top")+
  theme(legend.background = element_blank(), legend.key = element_blank())+
  theme(legend.title = element_text(size = 20, face = "bold", color = "gray5", hjust = 0))+
  theme(legend.text = element_text(size = 16, color = "gray5", face = "bold"))+
  theme(plot.title = element_text(size = 18, face = "bold", color = "gray5"))+
  theme(axis.ticks = element_line(size = 1))+
  theme(axis.title = element_text(size = 18, face = "bold", color = "gray5", hjust = 0.5))+
  theme(axis.text = element_text(size = 16, color = "gray5", hjust = 0.5))+
  theme(aspect.ratio = 1)#;ocab
  
# venn diagrmas of OTUs shared across size-fractions
df<-as.data.frame(unique(mm$cid)); colnames(df)<-"cid"
df$Micro<-ifelse(df$cid %in% rownames(micro) , 1, 0)
df$Nano<-ifelse(df$cid %in% rownames(nano) , 1, 0)
df$Pico<-ifelse(df$cid %in% rownames(pico) , 1, 0)
rownames(df)<-df$cid; df<-df[,-1]
mydata<-df %>% mutate_all(., as.logical)
vdc <- vennCounts(mydata)
class(vdc) <- 'matrix'
df.vdc <- as.data.frame(vdc)[-1,] %>%mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),y = c(1.2, -0.6, 0.5, -0.6, 0.5, -1, 0))
df.vdc$shared<-rowSums(df.vdc[,1:3])
df.venn <- data.frame(x = c(0, 0.866, -0.866),y = c(1, -0.5, -0.5),
                      labels = c('Micro', 'Nano', 'Pico'))

# venn plot
sf.venn<-ggplot(df.venn) +
  geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels), alpha = .4, size = 1, colour = 'gray55') +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "") +
  labs(x = "# of shared OTUs ", y = "")+
  scale_fill_manual(values = rev(c("#00fff7", "#494fc1", "#fd084a"))) +
  scale_colour_manual(values = rev(c("#00fff7", "#494fc1", "#fd084a")), guide = FALSE) +
  scale_x_continuous(limits = c(-3.2,3.2))+
  annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$Counts, size = 5)+
  annotate("text", x = c(0, 3, -3), y = c(3,-0.6, -0.6), label = c("Micro", "Nano", "Pico"), size = 6,color= rev(c("#00fff7", "#494fc1", "#fd084a")),  fontface = "bold")+
  theme(axis.title = element_text(size = 18, face = "bold", color = "gray5", hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0))#;sf.venn

## Diversity estimates with iNEXT
all<-data[,c(colnames(micro),colnames(nano),colnames(pico))];all<-all[rowSums(all) != 0, ]
lsf<-data.frame(All = rowSums(all), 
                Microplankton=rowSums(all[,colnames(micro)]),
                Nanoplankton=rowSums(all[,colnames(nano)]),
                Picoplankton=rowSums(all[,colnames(pico)]) ,
                a2009=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2009,"TAB"  ] ]),
                a2010=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2010,"TAB"  ] ]),
                a2011=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2011,"TAB"  ] ]),
                a2012=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2012,"TAB"  ] ]),
                a2013=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2013,"TAB"  ] ]),
                a2014=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2014,"TAB"  ] ]),
                a2015=rowSums(all[,colnames(all) %in% tab[format(tab$DATE, format = "%Y")==2015,"TAB"  ] ]) )
#div_esti<-iNEXT(x = lsf, datatype="abundance", endpoint = sum(lsf[,"All"])*4)
# saveRDS(div_esti,"POHEM/ANALYSIS/ASSEMBLY_PROCESS/iNEXT.RData")
div_esti<-readRDS("POHEM/ANALYSIS/ASSEMBLY_PROCESS/iNEXT.RData")

# study output
div_esti$AsyEst
div_esti$iNextEst

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
  scale_x_continuous(limits = c(-2e6,6.9e7), expand = c(0,0) )+
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
        axis.ticks = element_line(size = 1))#;grac

# EXport pdf 15 x 5
ocab<-ocab+ guides(color="none")
#(grac + ocab + sf.venn) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold'))

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
range_all<-rbind(range,data.frame(range = c(10,10,10),SF = c("Micro","Nano", "Pico"), name =c("Micro_AB","Nano_AB", "Pico_AB")))

# make a list of OTUs within ranges of the sliding window
list_range<-list()
for (i in 1:nrow(range)){
  list_range[[i]]<-subset(mm[mm$SF ==range[i,"SF"],], rank > range$range[i]-2500 & rank <= range$range[i]+2500)$cid
  names(list_range)[[i]]<-paste(range[i,2], range[i,1], sep = "_")
}

# Get list of abundant OTUs per SF
list_range_ab<-list()
for (i in c("Micro", "Nano", "Pico")){
  if(i == "Micro"){SF="Micro";thres=40}
  if(i == "Nano"){SF="Nano";thres=55}
  if(i == "Pico"){SF="Pico";thres=55}
  list_range_ab[[i]]<-mm[mm$SF == SF & mm$cumul > thres, "cid"]
};names(list_range_ab)<-paste(names(list_range_ab), "_AB", sep = "")
list_range_all<-c(list_range, list_range_ab) # merge both subset and abundant lists

inf<-bind_rows(lapply(list_range_all, function(x) table(x %in% rownames(traitf) )), .id = "subset")

# Information on the number of samples and env. samples inclueded per subset of OTUs
# list_all=NULL
# for (i in names(list_range_all)) {
#   if(range_all[range_all$name == i,"SF"] == "Micro" ){SF="Micro";thres=40; ps_dsb<-ps_micro}
#   if(range_all[range_all$name == i,"SF"] == "Nano" ){SF="Nano";thres=55;ps_dsb<-ps_nano}
#   if(range_all[range_all$name == i,"SF"] == "Pico" ){SF="Pico";thres=55; ps_dsb<-ps_pico}
# 
#   # data subset
#   sub_dsb<-phyloseq(otu_table(subset(otu_table(ps_dsb), rownames(otu_table(ps_dsb)) %in% list_range_all[[i]] )), phy_tree(phylo))
#   sub_dsb<-prune_samples(colSums(sub_dsb@otu_table) > 0, sub_dsb)
#   sub_dsb_env<-prune_samples(sample_names(sub_dsb) %in% rownames(env), sub_dsb)
#   list_all$subset[i]<-i
#   list_all$nb_samples[i]<-nsamples(sub_dsb)
#   list_all$nb_env_samples[i]<-nsamples(sub_dsb_env)
#   print(i)
# };list_all<-data.frame(nb_samples = list_all$nb_samples, nb_env_samples = list_all$nb_env_samples)
# saveRDS(list_all, "POHEM/ANALYSIS/ASSEMBLY_PROCESS/NB_SAMPLES_ALL.RData")
list_all<-readRDS("/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/NB_SAMPLES_ALL.RData")

######################################@######################################
# Taxonomic plots in sliding window
######################################@######################################

# format data of otus per subset
sub=NULL
for (i in names(list_range_all)){
  c<-data.frame(cid = list_range_all[[i]], id = i)
  sub<-rbind(sub,c)
}

# get taxonomy info of otus in each subset
sub_taxo<-merge(sub, taxo, by = "cid")
sub_taxo$SF<-gsub("_.*", "",sub_taxo$id)
sub_taxo$subset<-gsub(".*_", "",sub_taxo$id)
sub_taxo[sub_taxo$subset == "AB","subset"]<-10
sub_taxo$subset<-as.numeric(sub_taxo$subset)

# change value of last subset
last<-unique(sub_taxo[,c("SF","subset")]) %>% group_by(SF) %>% arrange(desc(subset)) %>% slice(2) 
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

# Exploration
unique(taxo[taxo$Division == "Sagenista","Class"])
#View(taxo[grep("Mamiellophyceae", taxo$Class),])
tot_taxo<-aggregate(ab ~ Division +SF+subset , data = sub_taxo, FUN =sum)
tot_taxo$perc<-round(((tot_taxo$ab)/sum(tot_taxo$ab))*100, 2)
tot_taxo[order(tot_taxo$perc),]
tail(tot_taxo[order(tot_taxo$perc),], 50)
#View(tot_taxo[tot_taxo$Division == "Unclassified",])
tt<-tot_taxo %>%
  group_by(SF, subset) %>%
  mutate(perc = round(ab/sum(ab)*100,2))
#View(tt[tt$Division == "Ochrophyta",])

# explore taxo of most abundant OTUs 
ab<-data.frame(unlist(list_range_ab));colnames(ab)<-"cid";ab$SF<-gsub("_.*", "",rownames(ab))
ab<-merge(ab, taxo, by = "cid")
ab<-merge(ab, mm, by = c("cid", "SF"))
ab$ab<-1
aggregate(ab ~ level +SF, data =ab, FUN = sum)

# plot
ag_sub_taxo<-aggregate(ab ~id + SF + subset + Division , data = sub_taxo, FUN =sum)
ag_sub_taxo$type<-ifelse(ag_sub_taxo$subset == 10, "Abundant", "Subset")
gstaxo<-ggplot(data = ag_sub_taxo[ag_sub_taxo$type == "Subset",] , aes(x = subset, y = ab, fill = Division))+
  facet_grid(SF~type, scale = "free_x", space = "free_x")+
  geom_bar(stat = "identity", width = 2000)+
  scale_y_continuous(expand =c(0,0), limits = c(0, 5500), breaks = seq(0,5000, 1000))+
  scale_x_continuous(expand =c(0,0), limits =c(0,75000), breaks = seq(0,80000, 10000))+
  scale_fill_manual(values = c(pal))+
  xlab("Subsets of Rank Abundance Curves (5000 OTUs)")+
  ylab("# OTUs")+
  ggtitle("Protistan Taxonomy across\nRank Abundance Curves")+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
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

#View(sub_taxo[sub_taxo$subset == 10,])
gstaxo_ab<-ggplot(data = ag_sub_taxo[ag_sub_taxo$type == "Abundant",] , aes(x = "Abundant OTUS", y = ab, fill = Division))+
  facet_grid(SF~type, scale = "free_x", space = "free_x")+
  geom_bar(stat = "identity", width = 2000)+
  scale_y_continuous(expand =c(0,0), limits = c(0, 45), breaks = seq(0,40, 10))+
  #scale_x_continuous(expand =c(0,0))+
  scale_fill_manual(values = c(pal), guide = "none")+
  xlab("")+
  ylab("# OTUs")+
  ggtitle("Protistan Taxonomy of\nmost Abundant OTUs")+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text =  element_blank())+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gstaxo_ab

######################################@######################################
# Taxonomic subsets at Division Level in each size-fraction
######################################@######################################

range

div<-data.frame(Division = rep(labtax$Division,3), 
           SF= c(rep("Micro", length(labtax$Division)),rep("Nano", length(labtax$Division)),rep("Pico", length(labtax$Division))))
div$name<-paste(div$SF,div$Division, sep = "_")
div$notus<-0

list_taxo<-list()
for (i in 1:nrow(div) ){
  list_taxo[[i]]<-mm[mm$SF %in% div[i,"SF"] & mm$cid %in% rownames(taxo[taxo$Division == div[i,"Division"],]) ,"cid"]
  div[i,"notus"]<-length(list_taxo[[i]])
  names(list_taxo)[[i]]<-div[i,3]
}
div<-div[div$notus > 20,]

# Create Ranges + division
range_div<-merge(range[,1:2], div[,1:2], by = "SF")
range_div$name<-paste(range_div$SF,range_div$range, range_div$Division, sep = "_" )

# make a list of OTUs within ranges of the sliding window
list_range_div<-list()
for (i in 1:nrow(range_div)){
  list_range_div[[i]]<-subset(mm[mm$SF ==range_div[i,"SF"]  & mm$cid %in% rownames(taxo[taxo$Division == range_div[i,"Division"],] ), ], rank > range_div$range[i]-2500 & rank <= range_div$range[i]+2500)$cid
  names(list_range_div)[[i]]<-range_div[i,4]
}
df_range_div <- data.frame(matrix(ncol = 2, nrow = length(list_range_div), data = NA));colnames(df_range_div)<-c("name", "notus")
for (i in 1:length(list_range_div)){
  df_range_div$name[i]<-names(list_range_div)[i]
  df_range_div$notus[i]<-length(list_range_div[[i]])
}
df_range_div$SF<-gsub("\\_.*","",df_range_div$name)
df_range_div$Division<-sub(".*_","",df_range_div$name)
df_range_div$range<-gsub(".*_(.+)_.*", "\\1", df_range_div$name)

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
  save.wd="/export/lv1/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/TEMP/"
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
  unlink("POHEM/ANALYSIS/ASSEMBLY_PROCESS/TEMP/", recursive = TRUE)

  # progress
  print(paste("done with = ", i, "/", nrow(range),sep = ""))
}

# save files
#saveRDS(object = list_sub_rc, file = "POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_SUBSET_RC.RData")
#saveRDS(object = list_sub_bnti, file = "POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_SUBSET_BNTI.RData")
#saveRDS(object = list_sub_bmntd, file = "POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_SUBSET_BMNTD.RData")

# read files
list_sub_rc<-readRDS( file = "POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_SUBSET_RC.RData")
list_sub_bnti<-readRDS( file = "POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_SUBSET_BNTI.RData")
list_sub_bmntd<-readRDS( file = "POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_SUBSET_BMNTD.RData")

# Abundant bNTI.big and RC.pc {iCAMP}: Loop to compute ßNTI and RC on abundant OTUS in phyloseq objects with a phylogenetic tree
list_abun_bnti<-list() # prepare list to receive ßNTI per subdataset
list_abun_rc<-list() # prepare list to receive RC per subdataset
list_abun_bmntd<-list() # prepare list to receive bMNTD per subdataset
for (i in c("Micro_AB", "Nano_AB", "Pico_AB") ) {

  if(i == "Micro_AB"){SF="Micro";thres=40; ps_dsb<-ps_micro}
  if(i == "Nano_AB"){SF="Nano";thres=55;ps_dsb<-ps_nano}
  if(i == "Pico_AB"){SF="Pico";thres=55; ps_dsb<-ps_pico}
  
  # data subset
  sub_dsb<-phyloseq(otu_table(subset(otu_table(ps_dsb), rownames(otu_table(ps_dsb)) %in% mm[mm$SF == SF & mm$cumul > thres, "cid"] )), phy_tree(phylo))
  sub_dsb<-prune_samples(colSums(sub_dsb@otu_table) > 0, sub_dsb)
  
  # make sure the names on the phylogeny are ordered the same as the names in otu table
  match.phylo.otu = match.phylo.data(sub_dsb@phy_tree, sub_dsb@otu_table)
  
  # New command to compute BNTI
  save.wd="/export/lv1/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/TEMP/"
  pd.big=pdist.big(tree = match.phylo.otu$phy, wd=save.wd, nworker = 40 )
  bNTI=bNTI.big(comm=t(match.phylo.otu$data), pd.desc=pd.big$pd.file,pd.spname=pd.big$tip.label,pd.wd=pd.big$pd.wd,
                spname.check=TRUE, nworker=40, weighted=TRUE, exclude.consp=FALSE,rand=999,
                output.dtail=TRUE, RC=FALSE, trace=FALSE)
  
  # New command to compute Raup-Crick
  RC=RC.pc(comm=t(match.phylo.otu$data), rand = 999,
           nworker = 40, weighted = TRUE,
           sig.index="RC", silent = TRUE)
  
  # save files in list
  list_abun_bnti[[i]]<-bNTI$bNTI
  list_abun_bmntd[[i]]<-bNTI$bMNTD
  list_abun_rc[[i]]<-RC

  # remove by-products
  unlink("POHEM/ANALYSIS/ASSEMBLY_PROCESS/TEMP/", recursive = TRUE)

  # progress
  print(paste("done with = ", i, "/3",sep = ""))
}
# save files
#saveRDS(object = list_abun_rc, file = "POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_ABUN_RC.RData")
#saveRDS(object = list_abun_bnti, file = "POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_ABUN_BNTI.RData")
#saveRDS(object = list_abun_bmntd, file = "POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_ABUN_BMNTD.RData")

# read files
list_abun_bnti<-readRDS("POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_ABUN_BNTI.RData")
list_abun_bmntd<-readRDS("POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_ABUN_BMNTD.RData")
list_abun_rc<-readRDS("POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_ABUN_RC.RData")

######################################@######################################
# Infer and plot Assembly Processes from phylogenetic turnovers
######################################@######################################

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

# Loop Assembly process inference across subsets 
list_abun_ass<-list() # prepare list to receive the proportions of assembly processes per subdataset
for (i in c("Micro_AB", "Nano_AB", "Pico_AB") ) {
  
  # extract and merge ßNTI and ßRC
  eco<-merge(melt(list_abun_bnti[[i]]),melt(list_abun_rc[[i]]$index), by = c("Var1", "Var2")); colnames(eco)<-c("sp1", "sp2", "bNTI", "bRCBC")
  
  # infer assembly process per community comparison pairs
  eco$`Dominant Ecological Process`<-ifelse(eco$bNTI >= 2, "Variable Selection", 
                                            ifelse(eco$bNTI <= -2, 'Homogeneous Selection',
                                                   ifelse(abs(eco$bNTI) <= 2 & eco$bRCBC >= 0.95, "Dispersal Limitation",
                                                          ifelse(abs(eco$bNTI) <= 2 & eco$bRCBC <= -0.95, "Homogenizing Dispersal",
                                                                 ifelse(abs(eco$bNTI) <= 2 & abs(eco$bRCBC) <= 0.95, "Undominated Scenario", NA)))))
  
  # save files in list
  list_abun_ass[[i]]<-eco
  remove(eco)
  
  # progress
  print(paste("done with = ", i,sep = ""))
}

# merge and unlist all datasets
list_sub_ass<-c(list_sub_ass,list_abun_ass)
sub_ass<-bind_rows(list_sub_ass, .id = "id") # unlist every dataset into a merged object
sub_ass[is.na(sub_ass$`Dominant Ecological Process`),]
sub_ass<-sub_ass[sub_ass$sp1 != sub_ass$sp2,]

# count assembly processes per subset
sub_ass_count<-sub_ass %>% count(id, `Dominant Ecological Process`)
sub_ass_count$SF<-gsub("_.*", "",sub_ass_count$id)
sub_ass_count$subset<-gsub(".*_", "",sub_ass_count$id)
sub_ass_count[sub_ass_count$subset == "AB","subset"]<-"10"
sub_ass_count$subset<-as.numeric(sub_ass_count$subset)
sub_ass_count<-na.omit(sub_ass_count) 
sub_ass_count$`Dominant Ecological Process`<-factor(sub_ass_count$`Dominant Ecological Process`, levels = rev(c("Undominated Scenario","Homogenizing Dispersal", "Dispersal Limitation","Homogeneous Selection", "Variable Selection")), ordered = TRUE)
# NA occurs increasingly towards rare subsets 
# because less samples contain OTUs (a sample without OTUs cannot be compared to other samples)

# compute proportions %
sub_ass_count$tot_sub<-NA
sub_ass_count$perc<-NA
for (i in 1:nrow(sub_ass_count) ){
  sub_ass_count[i, "tot_sub"]<-sum(sub_ass_count[sub_ass_count$id %in% sub_ass_count[i,"id"],"n"])
  sub_ass_count[i, "perc"]<-(sub_ass_count[i, "n"]/sub_ass_count[i, "tot_sub"])*100
}

# change value of last subset
last<-unique(sub_ass_count[,c("SF","subset")]) %>% group_by(SF) %>% arrange(desc(subset)) %>% slice(2) 
last$last<-last$subset+2500
last<-merge(last, aggregate(subset ~ SF,data = sub_ass_count, FUN = max), by = "SF")
for (i in 1:nrow(last) ){
  sub_ass_count[sub_ass_count$subset == last[i, 4], "subset"]<-last[i, 3]
}
sub_ass_count$data<-ifelse(sub_ass_count$subset == 10, "Abundant OTUs", "Rank Abundance Curve")
sub_ass_count$data<-factor(sub_ass_count$data, levels = c("Rank Abundance Curve", "Abundant OTUs"), ordered = TRUE)
sub_ass_count<-merge(sub_ass_count, list_all, by.x = "id", by.y = "row.names")
levels(sub_ass_count$`Dominant Ecological Process`)

# format assembly processes name and types of rarity
sub_ass_count$`Assembly Processes and Types of Rarity`<-ifelse(sub_ass_count$`Dominant Ecological Process` == "Variable Selection", "Variable Selection\nConditionally Rare",
       ifelse(sub_ass_count$`Dominant Ecological Process` == "Homogeneous Selection", "Homogeneous Selection\nPermanently Rare",
              ifelse(sub_ass_count$`Dominant Ecological Process` == "Dispersal Limitation", "Dispersal Limitation\nTransiently Rare",
                     ifelse(sub_ass_count$`Dominant Ecological Process` == "Homogenizing Dispersal", "Homogenizing Dispersal\nPermanently Rare (with periodic distribution)",
                            ifelse(sub_ass_count$`Dominant Ecological Process` == "Undominated Scenario", "Undominated Scenario\nTransiently Rare",NA)))))
sub_ass_count$`Assembly Processes and Types of Rarity`<-factor(sub_ass_count$`Assembly Processes and Types of Rarity`, levels = c("Variable Selection\nConditionally Rare", "Homogeneous Selection\nPermanently Rare", "Dispersal Limitation\nTransiently Rare", "Homogenizing Dispersal\nPermanently Rare (with periodic distribution)","Undominated Scenario\nTransiently Rare"), ordered = TRUE)
levels(sub_ass_count$`Assembly Processes and Type of Rarity`)

# plot the assembly processes per subdataset and SF
ag_sub_ass<-aggregate(perc ~ subset+ id+SF+nb_samples,data = sub_ass_count, FUN = sum)
#View(sub_ass_count[sub_ass_count$`Dominant Ecological Process`=="Homogenous Selection",])
#aggregate(perc ~ SF,data = sub_ass_count[sub_ass_count$`Dominant Ecological Process`=="Homogenous Selection",], FUN = mean)
gsass<-ggplot(data = sub_ass_count[sub_ass_count$data != "Abundant OTUs",], aes(x = subset, y = perc, fill = `Assembly Processes and Types of Rarity`))+
  facet_grid(SF~., scale = "free_x", space = "free")+
  geom_bar(stat = "identity", width = 2000)+
  #geom_line()+
  geom_text(data = ag_sub_ass[ag_sub_ass$subset != 10,], aes(x =subset, y = 105 , group =id , label = nb_samples), inherit.aes = FALSE, size = 3.00, vjust =0)+
  #geom_text(aes(label = nb_samples, group = id, y = 1), angle = 0, hjust = 0.5, vjust = 0 , color = "gray15", size= 3.00)+
  scale_y_continuous(expand =c(0,0), limits = c(0, 115), breaks = seq(0,100, 25))+
  scale_x_continuous(expand =c(0,0), limits =c(0,75000), breaks = seq(0,80000, 10000))+
  scale_fill_manual(values = c("#063175","#7aadff", "#a60202", "#ff491c", "#cccccc"))+
  xlab("Subsets of Rank Abundance Curves (5000 OTUs)")+ylab("Frequency (%)")+
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

gsass_ab<-ggplot(data = sub_ass_count[sub_ass_count$data == "Abundant OTUs",], aes(x = "Abundant OTUs", y = perc, fill = `Assembly Processes and Types of Rarity`))+
  facet_grid(SF~., scale = "free_x", space = "free")+
  geom_bar(stat = "identity", width = 2000)+
  scale_y_continuous(expand =c(0,0), limits = c(0, 115), breaks = seq(0,100, 25))+
  geom_text(data = ag_sub_ass[ag_sub_ass$subset == 10,], aes(x =subset, y = 105 , group =id , label = nb_samples), inherit.aes = FALSE, size = 3.00, vjust =0)+
  #scale_x_continuous(expand =c(0,0), breaks = seq(0,80000, 10000),limits = c(0,79000))+
  scale_fill_manual(values = c("#063175","#7aadff", "#a60202", "#ff491c", "#cccccc"), guide ="none")+
  xlab("")+ylab("Frequency (%)")+
  ggtitle("Assembly Processes of\nmost Abundant OTUs")+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))+
  #theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 1))+
  theme(axis.text = element_text(size = 16))+
  theme(strip.text =  element_blank())+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gsass_ab

######################################@######################################
# Get Permanova for each sub community
######################################@######################################

# merge bmntd from subsets and abundant communities
list_bmntd<-append(list_sub_bmntd, list_abun_bmntd)

# Loop PERMANOVA across subsets (only for samples that have env. conditions)
list_sub_perma<-list()
list_sub_perma_bray<-list()
for (i in names(list_bmntd) ) {

  # ßMNTD distance
  common<-colnames(list_bmntd[[i]])[colnames(list_bmntd[[i]]) %in% rownames(na.omit(env))]
  comm<-list_bmntd[[i]][common,common]
  
  # Bray curtis distance
  if(range_all[range_all$name == i,"SF"] == "Micro" ){SF="Micro";thres=40; dsb<-micro}
  if(range_all[range_all$name == i,"SF"] == "Nano" ){SF="Nano";thres=55;dsb<-nano}
  if(range_all[range_all$name == i,"SF"] == "Pico" ){SF="Pico";thres=55; dsb<-pico}
  cbray<- dsb[rownames(dsb) %in% list_range_all[[i]], common]
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
  print(paste("done with = ", match(i, names(list_bmntd)), "/", length(list_bmntd),sep = ""))
}

# save files
#saveRDS(object = list_sub_perma, file = "POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_PERMANOVA.RData")
#saveRDS(object = list_sub_perma_bray, file = "POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_PERMANOVA_BC.RData")
# read files
list_sub_perma<-readRDS("POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_PERMANOVA.RData")

# unlist all datasets
sub_perma<-bind_rows(list_sub_perma, .id = "id") # unlist every dataset into a merged object
sub_perma$env<-gsub("\\..*","",rownames(sub_perma))
sub_perma$SF<-gsub("_.*", "",sub_perma$id)
sub_perma$subset<-gsub(".*_", "",sub_perma$id)
sub_perma[sub_perma$subset == "AB","subset"]<-10
sub_perma$subset<-as.numeric(sub_perma$subset)


ass_df<-acast(id~`Dominant Ecological Process`, data = sub_ass_count, value.var = "perc")
env_df<-acast(id~env, data = sub_perma, value.var = "R2")
ass_env<-merge(ass_df, env_df, by = "row.names")
ass_env$subset<-sub(".*_","",ass_env$Row.names)
ass_env$subset<-ifelse(ass_env$subset == "AB",10, as.numeric(ass_env$subset) )
ass_env$SF<-gsub("_.*", "",ass_env$Row.names)

mx_sf<-aggregate(subset ~ SF,ass_env ,max );colnames(mx_sf)<-c("SF", "max")
ass_env<-merge(ass_env,mx_sf , by = "SF")
ass_env$per_rac<-(as.numeric(ass_env$subset)/ass_env$max)*100
ass_env$per_rac<-factor(ass_env$per_rac, levels=sort(as.numeric(unique(ass_env$per_rac))), ordered = TRUE )

ggplot(data= ass_env[grep("AB", ass_env$Row.names, invert = TRUE),], aes(x = Total-Residual, y = `Variable Selection`, color= per_rac))+
  scale_color_viridis(direction = -1, name ="Rank along RAC\n(in % of RAC per SF)" )+
  facet_grid(SF~.)+
  xlab("Explaining power of Env. var.\n(PERMANOVA, [0-1])")+
  ylab("Proportion of\nVariable Selection (%)")+
  geom_point(size = 3, alpha = 0.75)+
  theme(panel.spacing.y  = unit(0.25, "lines"))+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text.y =  element_text(size = 16, face = "bold", angle =0, hjust = 0))+
  theme(strip.text.x =  element_text(size = 16, face = "bold", angle =0, vjust =0, hjust = 0.5 ))+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0))

# Format for better plotting
last<-unique(sub_perma[,c("SF","subset")]) %>% group_by(SF) %>% arrange(desc(subset)) %>% slice(2) 
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
View(aggregate(abs(R2) ~ id, data = sub_perma, FUN = sum))
View(sub_perma)
gsperma<-ggplot(data = sub_perma[sub_perma$subset !=10,], aes(x = subset, y = abs(R2), fill = env_cond))+
  facet_grid(SF~.)+
  geom_bar(stat = "identity", width = 2000)+
  scale_fill_manual(name = "Environmental\nVariable",values = pal_p)+
  scale_y_continuous(expand =c(0,0), limits = c(0, 0.45))+
  scale_x_continuous(expand =c(0,0), limits =c(0,75000), breaks = seq(0,80000, 10000))+
  geom_text(data = ag_perma[ag_perma$subset != 10,], aes(x =subset, y = `abs(R2)`+0.025 , group =id , label = nb_env_samples), inherit.aes = FALSE, size = 3.00)+
  #geom_text(aes(label = nb_env_samples, group = id, y = 0.32), angle = 0, hjust = 0.5, vjust = 0 , color = "gray15", size= 3.00)+
  xlab("Subsets of Rank Abundance Curves (5000 OTUs)")+ylab("PERMANOVA's R2")+
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

gsperma_ab<-ggplot(data = sub_perma[sub_perma$subset==10,], aes(x = "Abundant OTUs", y = abs(R2), fill = env_cond))+
  facet_grid(SF~.)+
  geom_bar(stat = "identity", width = 2000)+
  scale_fill_manual(name = "Environmental\nVariable",values = pal_p, guide = FALSE)+
  scale_y_continuous(expand =c(0,0), limits = c(0, 0.45))+
  geom_text(data = ag_perma[ag_perma$subset == 10,], aes(x =subset, y = `abs(R2)`+0.025 , group =id , label = nb_env_samples), inherit.aes = FALSE, size = 3.00)+
  xlab("")+ylab("PERMANOVA's R2")+
  ggtitle("PERMANOVA of\nmost Abundant OTUs")+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 18,hjust = 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text = element_blank())+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gsperma_ab

library(patchwork)
#gstaxo + gsass + gsperma + plot_layout(ncol = 1) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold'))

p1<-gstaxo_ab+ggtitle("") + gstaxo+ggtitle("")+ylab("") +  plot_layout(design = c("ABBBBBBBBBBB"))+ plot_annotation(title = 'A: Protistan Taxonomy of Abundant OTUs\nand along Rank Abundance Curves') & theme(title = element_text(size = 16, face = "bold"))
p2<-gsass_ab+ggtitle("") + gsass+ggtitle("")+ylab("") +  plot_layout(design = c("ABBBBBBBBBBB"))+ plot_annotation(title = 'B: Assembly Processes of Abundant OTUs\nand along Rank Abundance Curves') & theme(title = element_text(size = 16, face = "bold"))
p3<-gsperma_ab+ggtitle("") + gsperma+ggtitle("") +ylab("")+  plot_layout(design = c("ABBBBBBBBBBB"))+ plot_annotation(title = 'C: PERMANOVA of Abundant OTUs\nand along Rank Abundance Curves') & theme(title = element_text(size = 16, face = "bold"))

ggpubr::ggarrange(p1,p2,p3, ncol = 1, align = "hv")

######################################@######################################
# Compute ßNTI and RC for taxo subsets
######################################@######################################

# bNTI.big and RC.pc {iCAMP}: Loop to compute ßNTI and RC on subdatasets as phyloseq objects with a phylogenetic tree
list_div_bnti<-list() # prepare list to receive ßNTI per subdataset
list_div_rc<-list() # prepare list to receive RC per subdataset
list_div_bmntd<-list() # prepare list to receive bMNTD per subdataset

# Remove abundant taxa
ab<-unique(unlist(list_range_ab))
list_taxo1<-list()
for (i in 1:length(list_taxo)){
  list_taxo1[[names(list_taxo)[i]]]<-setdiff(list_taxo[[names(list_taxo)[i] ]], ab)
}

# LOOP
for (i in 1:nrow(div) ) {
  
  # dataset selection according to size-fraction
  if(div[i,"SF"] == "Micro"){ps_dsb<-ps_micro}
  if(div[i,"SF"] == "Nano"){ps_dsb<-ps_nano}
  if(div[i,"SF"] == "Pico"){ps_dsb<-ps_pico}
  
  # data subset
  sub_dsb<-phyloseq(otu_table(subset(otu_table(ps_dsb), rownames(otu_table(ps_dsb)) %in% list_taxo1[[div[i, "name"]]] )), phy_tree(phylo))
  sub_dsb<-prune_samples(colSums(sub_dsb@otu_table) > 0, sub_dsb)
  
  # make sure the names on the phylogeny are ordered the same as the names in otu table
  match.phylo.otu = match.phylo.data(sub_dsb@phy_tree, sub_dsb@otu_table)
  
  # New command to compute BNTI
  save.wd="/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/TEMP1/"
  pd.big=pdist.big(tree = match.phylo.otu$phy, wd=save.wd, nworker = 40 )
  bNTI=bNTI.big(comm=t(match.phylo.otu$data), pd.desc=pd.big$pd.file,pd.spname=pd.big$tip.label,pd.wd=pd.big$pd.wd,
                spname.check=TRUE, nworker=40, weighted=TRUE, exclude.consp=FALSE,rand=999,
                output.dtail=TRUE, RC=FALSE, trace=FALSE)
  
  # New command to compute Raup-Crick
  RC=RC.pc(comm=t(match.phylo.otu$data), rand = 999,
           nworker = 40, weighted = TRUE,
           sig.index="RC", silent = TRUE)
  
  # save files in list
  list_div_bnti[[i]]<-bNTI$bNTI
  names(list_div_bnti)[[i]]<-div[i,3]
  list_div_bmntd[[i]]<-bNTI$bMNTD
  names(list_div_bmntd)[[i]]<-div[i,3]
  list_div_rc[[i]]<-RC
  names(list_div_rc)[[i]]<-div[i,3]
  
  # remove by-products
  unlink("/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/TEMP1/", recursive = TRUE)
  
  # progress
  print(paste("done with = ", i, "/", nrow(div),sep = ""))
}

# save files
#saveRDS(object = list_div_rc, file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_RC.RData")
#saveRDS(object = list_div_bnti, file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_BNTI.RData")
#saveRDS(object = list_div_bmntd, file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_BMNTD.RData")

#saveRDS(object = list_div_rc, file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_RC_NAB.RData")
#saveRDS(object = list_div_bnti, file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_BNTI_NAB.RData")
#saveRDS(object = list_div_bmntd, file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_BMNTD_NAB.RData")

# Import files
list_div_rc<-readRDS( file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_RC.RData")
list_div_bnti<-readRDS( file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_BNTI.RData")
list_div_bmntd<-readRDS( file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_BMNTD.RData")

# Import files
list_div_rc<-readRDS( file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_RC_NAB.RData")
list_div_bnti<-readRDS( file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_BNTI_NAB.RData")
list_div_bmntd<-readRDS( file = "/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/FIN_DIVISION_BMNTD_NAB.RData")


# Loop Assembly process inference across subsets 
list_div_ass<-list() # prepare list to receive the proportions of assembly processes per subdataset
for (i in 1:nrow(div) ) {
  
  # extract and merge ßNTI and ßRC
  eco<-merge(melt(list_div_bnti[[i]]),melt(list_div_rc[[i]]$index), by = c("Var1", "Var2")); colnames(eco)<-c("sp1", "sp2", "bNTI", "bRCBC")
  
  # infer assembly process per community comparison pairs
  eco$`Dominant Ecological Process`<-ifelse(eco$bNTI >= 2, "Variable Selection", 
                                            ifelse(eco$bNTI <= -2, 'Homogeneous Selection',
                                                   ifelse(abs(eco$bNTI) <= 2 & eco$bRCBC >= 0.95, "Dispersal Limitation",
                                                          ifelse(abs(eco$bNTI) <= 2 & eco$bRCBC <= -0.95, "Homogenizing Dispersal",
                                                                 ifelse(abs(eco$bNTI) <= 2 & abs(eco$bRCBC) <= 0.95, "Undominated Scenario", NA)))))
  
  # save files in list
  list_div_ass[[i]]<-eco
  names(list_div_ass)[[i]]<-div[i,3]
  remove(eco)
  
  # progress
  print(paste("done with = ", i, "/", nrow(div),sep = ""))
}

# format data
div_ass<-bind_rows(list_div_ass, .id = "id") # unlist every dataset into a merged object
div_ass[is.na(div_ass$`Dominant Ecological Process`),]
div_ass<-div_ass[div_ass$sp1 != div_ass$sp2,]

# count assembly processes per subset
div_ass_count<-div_ass %>% count(id, `Dominant Ecological Process`)
div_ass_count$SF<-gsub("_.*", "",div_ass_count$id)
div_ass_count$Division<-factor(gsub(".*_", "",div_ass_count$id), levels = labtax$Division, ordered = TRUE)
div_ass_count<-na.omit(div_ass_count) 
div_ass_count$`Dominant Ecological Process`<-factor(div_ass_count$`Dominant Ecological Process`, levels = rev(c("Undominated Scenario","Homogenizing Dispersal", "Dispersal Limitation","Homogeneous Selection", "Variable Selection")), ordered = TRUE)
# NA occurs increasingly towards rare subsets 
# because less samples contain OTUs (a sample without OTUs cannot be compared to other samples)

# compute proportions %
div_ass_count$tot_sub<-NA
div_ass_count$perc<-NA
for (i in 1:nrow(div_ass_count) ){
  div_ass_count[i, "tot_sub"]<-sum(div_ass_count[div_ass_count$id %in% div_ass_count[i,"id"],"n"])
  div_ass_count[i, "perc"]<-(div_ass_count[i, "n"]/div_ass_count[i, "tot_sub"])*100
}
div_ass_count<-merge(div_ass_count, labtax, by = "Division")
head(div_ass_count)
# format assembly processes name and types of rarity
div_ass_count$`Assembly Processes and Types of Rarity`<-ifelse(div_ass_count$`Dominant Ecological Process` == "Variable Selection", "Variable Selection\nConditionally Rare",
                                                               ifelse(div_ass_count$`Dominant Ecological Process` == "Homogeneous Selection", "Homogeneous Selection\nPermanently Rare",
                                                                      ifelse(div_ass_count$`Dominant Ecological Process` == "Dispersal Limitation", "Dispersal Limitation\nTransiently Rare",
                                                                             ifelse(div_ass_count$`Dominant Ecological Process` == "Homogenizing Dispersal", "Homogenizing Dispersal\nPermanently Rare (with periodic distribution)",
                                                                                    ifelse(div_ass_count$`Dominant Ecological Process` == "Undominated Scenario", "Undominated Scenario\nTransiently Rare",NA)))))
div_ass_count$`Assembly Processes and Types of Rarity`<-factor(div_ass_count$`Assembly Processes and Types of Rarity`, levels = c("Variable Selection\nConditionally Rare", "Homogeneous Selection\nPermanently Rare", "Dispersal Limitation\nTransiently Rare", "Homogenizing Dispersal\nPermanently Rare (with periodic distribution)","Undominated Scenario\nTransiently Rare"), ordered = TRUE)
levels(div_ass_count$`Assembly Processes and Type of Rarity`)

# plot the assembly processes per subdataset and SF
str(div_ass_count)
table(labtax[,2])==1
View(div_ass_count[div_ass_count$Division == "Sagenista",])
View(aggregate(perc ~ SF+Division , data = div_ass_count[div_ass_count$`Dominant Ecological Process` %in% c("Homogenizing Dispersal", "Homogeneous Selection"),], FUN =sum))
div_ass_count[div_ass_count$`Dominant Ecological Process` %in% c("Variable Selection"),]

gdass<-ggplot(data = div_ass_count, aes(x = perc, y = Division, fill = `Assembly Processes and Types of Rarity`))+
  facet_grid(Supergroup~SF, scale = "free", space = "free")+
  geom_bar(stat = "identity")+
  #geom_text(data = ag_sub_ass[ag_sub_ass$subset != 10,], aes(x =subset, y = 105 , group =id , label = nb_samples), inherit.aes = FALSE, size = 3.00, vjust =0)+
  #geom_text(aes(label = nb_samples, group = id, y = 1), angle = 0, hjust = 0.5, vjust = 0 , color = "gray15", size= 3.00)+
  scale_x_continuous(expand =c(0,0), limits = c(0, 110), breaks = seq(0,100, 25))+
  scale_fill_manual(values = c("#063175","#7aadff", "#a60202", "#ff491c", "#cccccc"))+
  xlab("Frequency (%)")+ylab("")+
  ggtitle("Assembly Processes across\nProtistan Divisions")+
  theme(panel.spacing.y  = unit(0.25, "lines"))+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text.y =  element_text(size = 16, face = "bold", angle =0, hjust = 0))+
  theme(strip.text.x =  element_text(size = 16, face = "bold", angle =0, vjust =0, hjust = 0.5 ))+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gdass

ga<-gdass+ggtitle("A: Assembly Processes across\nProtistan Divisions\nWith Abundant OTUs")+theme(strip.text.y = element_blank())+theme(legend.position = "none")
gb<-gdass1+ggtitle("B: Assembly Processes across\nProtistan Divisions\nWithout Abundant OTUs")+theme(axis.text.y = element_blank(),axis.title.y = element_blank())
ga+gb

######################################@######################################
# position on RAC of otus for each Division
######################################@######################################

df_range_div<-merge(df_range_div, labtax, by = "Division")
mx_sf<-aggregate(range ~ SF, range,max );colnames(mx_sf)<-c("SF", "max")
df_range_div<-merge(df_range_div,mx_sf , by = "SF")
df_range_div$per_rac<-(as.numeric(df_range_div$range)/df_range_div$max)*100
df_range_div$per_rac<-factor(df_range_div$per_rac, levels=sort(as.numeric(unique(df_range_div$per_rac))), ordered = TRUE )
df_range_div<-df_range_div[order(df_range_div$per_rac),]

ag_sf_div<-aggregate(notus ~ Division+SF,df_range_div, sum)

ggdf<-ggplot(data = df_range_div, aes(x = notus, y = Division, fill = per_rac))+
  facet_grid(Supergroup~SF, scale = "free", space = "free")+
  geom_bar(stat = "identity", position = "fill")+
  scale_fill_viridis_c(direction = -1, name ="Rank along RAC\n(in % of RAC per SF)" )+
  xlab("# OTUs")+ylab("")+
  ggtitle("Assembly Processes across\nProtistan Divisions")+
  theme(panel.spacing.y  = unit(0.25, "lines"))+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text.y =  element_text(size = 16, face = "bold", angle =0, hjust = 0))+
  theme(strip.text.x =  element_text(size = 16, face = "bold", angle =0, vjust =0, hjust = 0.5 ))+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));ggdf

(ggdf+ggdf1) + plot_layout(guides = "collect")

######################################@######################################
# Compute ßNTI and RC for one set of taxa
######################################@######################################

for (i in c("Micro", "Nano", "Pico") ) {
  
  # dataset selection according to size-fraction
  if(i == "Micro"){ps_dsb<-ps_micro}
  if(i == "Nano"){ps_dsb<-ps_nano}
  if(i == "Pico"){ps_dsb<-ps_pico}
  
  # data subset
  sub_dsb<-phyloseq(otu_table(subset(otu_table(ps_dsb), rownames(otu_table(ps_dsb)) %in% taxo[taxo$Class == "Labyrinthulomycetes", "cid"] )), phy_tree(phylo))
  sub_dsb<-prune_samples(colSums(sub_dsb@otu_table) > 0, sub_dsb)
  
  # make sure the names on the phylogeny are ordered the same as the names in otu table
  match.phylo.otu = match.phylo.data(sub_dsb@phy_tree, sub_dsb@otu_table)
  
  # New command to compute BNTI
  save.wd="/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/TEMP1/"
  pd.big=pdist.big(tree = match.phylo.otu$phy, wd=save.wd, nworker = 40 )
  bNTI=bNTI.big(comm=t(match.phylo.otu$data), pd.desc=pd.big$pd.file,pd.spname=pd.big$tip.label,pd.wd=pd.big$pd.wd,
                spname.check=TRUE, nworker=40, weighted=TRUE, exclude.consp=FALSE,rand=999,
                output.dtail=TRUE, RC=FALSE, trace=FALSE)
  
  # New command to compute Raup-Crick
  RC=RC.pc(comm=t(match.phylo.otu$data), rand = 999,
           nworker = 40, weighted = TRUE,
           sig.index="RC", silent = TRUE)
  
  # save files in list
  list_taxo_bnti[[i]]<-bNTI$bNTI
  names(list_taxo_bnti)[[i]]<-i
  list_taxo_bmntd[[i]]<-bNTI$bMNTD
  names(list_taxo_bmntd)[[i]]<-i
  list_taxo_rc[[i]]<-RC
  names(list_taxo_rc)[[i]]<-i
  
  # remove by-products
  unlink("/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/TEMP1/", recursive = TRUE)
  
  # progress
  print(paste("done with = ", i, "/", nrow(div),sep = ""))
}



######################################@##################################################
# Compute ßNRI across size-fraction: delineate assembly processes across protistan phyla
######################################@############################################@#####

# loop iCAMP command across size-fraction datasets
list_bnri<-list
for (i in c("Micro", "Nano", "Pico") ) {
  
  if(i == "Micro"){SF="Micro"; ps_dsb<-ps_micro}
  if(i == "Nano"){SF="Nano";ps_dsb<-ps_nano}
  if(i == "Pico"){SF="Pico";ps_dsb<-ps_pico}
  
  # parameter for iCAMP command
  save.wd="/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/TEMP/"
  pd.wd=paste0(save.wd,"/pdbig")
  nworker=40 # parallel computing thread number
  rand.time=999 # usually use 1000 for real data.
  bin.size.limit=100 # min nb of OTUs in phylogenetic bins
  
  # make sure the names on the phylogeny are ordered the same as the names in otu table
  match.phylo.otu = match.phylo.data(ps_dsb@phy_tree, ps_dsb@otu_table)
  
  # iCAMP commands
  setwd(save.wd)
  icamp.out=icamp.big(comm=t(match.phylo.otu$data),tree=match.phylo.otu$phy,pd.wd=pd.wd,
                    rand=rand.time, nworker=nworker,
                    bin.size.limit=bin.size.limit)
  
  # save files in list
  list_bnri[[i]]<-icamp.out
  
  # erase temporary directory
  unlink("/export/lv6/user/pramond/POHEM/ANALYSIS/ASSEMBLY_PROCESS/TEMP/pdbig", recursive = TRUE)
}

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
# Compute functional richness across subsets of rank abundance curves
############################################################################################################

## empty table referencing OTUs and rank abundance subsets
sub_table<-matrix(ncol = length(unique(mm$cid)), nrow = length(list_range_all) )
colnames(sub_table)<-unique(mm$cid);rownames(sub_table)<-names(list_range_all)
# fill the table according to the occurence of each OTU in each susbset
for (i in rownames(sub_table) ){
  sub_table[i,]<-ifelse(colnames(sub_table) %in% list_range_all[[i]], 1, 0  )
  print(i)
}
sub_table_f<-sub_table[, colnames(sub_table) %in% rownames(traitf)];rowSums(sub_table_f)
sbsb<-merge(t(sub_table_f), trait[,c(1,2)], by.x = "row.names", by.y = "cid")
# aggregate per taxonomic reference
sbsb<-aggregate(x = sbsb[,-c(1, ncol(sbsb))], by = list(sbsb$lineage), FUN = sum ); rownames(sbsb)<-sbsb$Group.1; sbsb<-sbsb[,-1]; sbsb<-as.data.frame(t(sbsb))

## format trait table
# turn bivariable column into numeric
for (i in colnames(traitf_f[,apply(traitf_f, 2, function(x) length(unique(x)))==2])){
  traitf_f[,i]<-as.numeric(traitf_f[,i])
}
trtr<-merge(traitf_f, trait[,c(1,2)], by.x = "row.names", by.y = "cid")
trtr<-unique(trtr[,-c(1)]); rownames(trtr)<-trtr$lineage; trtr<-trtr[,-ncol(trtr)]

## order dataset as preparation for functional diversity functions
sbsb_f<-as.data.frame(sbsb[,colnames(sbsb) %in% rownames(trtr)])
trtr_f<-trtr[match(colnames(sbsb_f), rownames(trtr)),]

## compute Functional diversity indexes per subset
#FD<-dbFD(trtr_f, sbsb_f, corr = "cailliez", stand.FRic = TRUE) # takes about 30min
#saveRDS(object = FD, file = "/export/lv6/user/pramond/POHEM/ANALYSIS/TRAIT_RARE/FD_SUBSET.RData")
#FD<-readRDS("/export/lv6/user/pramond/POHEM/ANALYSIS/TRAIT_RARE/FD_SUBSET.RData")
# format FD data for plotting
# FD[["sub"]]<-names(FD$nbsp)
# FD<-as.data.frame(bind_rows(FD))
# FD$sub

## output of function from mFD package computed locally (could not install package on servers)
FD<-readRDS("/export/lv6/user/pramond/POHEM/ANALYSIS/TRAIT_RARE/mFD_pa.RData")
FD<-FD$functional_diversity_indices[grep("fide", colnames(FD$functional_diversity_indices), invert = TRUE)]
FD$sub<-rownames(FD)

# format and retrieve information
FD$SF<-gsub("_.*", "",FD$sub)
FD$subset<-gsub(".*_", "",FD$sub)
FD[FD$subset == "AB","subset"]<-"10"
FD$subset<-as.numeric(FD$subset)

# change value of last subset
last<-unique(FD[,c("SF","subset")]) %>% group_by(SF) %>% arrange(desc(subset)) %>% slice(2) 
last$last<-last$subset+2500
last<-merge(last, aggregate(subset ~ SF,data = FD, FUN = max), by = "SF")
for (i in 1:nrow(last) ){
  FD[FD$subset == last[i, 4], "subset"]<-last[i, 3]
}
FD$data<-ifelse(FD$subset == 10, "Abundant OTUs", "Rank Abundance Curve")
FD$data<-factor(FD$data, levels = c("Rank Abundance Curve", "Abundant OTUs"), ordered = TRUE)
FD<-merge(FD, list_all, by.x = "sub", by.y = "row.names")
sub_notus<-as.data.frame(sapply(list_range_all, FUN = length));colnames(sub_notus)<-"notus"
FD<-merge(FD, sub_notus, by.x = "sub", by.y = "row.names")
FD<-FD[FD$notus == 5000,]

summary(FD$fric)
summary(FD$fdis)
aggregate(fric~SF, data = FD, FUN = mean)
aggregate(fdis~SF, data = FD, FUN = mean)

## plot FD
# top 3 values per index
pairs(FD[,c(2:10)])
top3_fr<-FD %>% group_by(SF) %>% top_n(3, fric) # top3 values of FRic across SF
top3_fdis<-FD %>% group_by(SF) %>% top_n(3, fdis) # top3 values of FDiv across SF

# format data
mFD<-reshape2::melt(FD[,-c(1:2)] , id.vars = c("SF","subset","data","nb_samples","nb_env_samples", "notus") )
mFD$variable<-factor(mFD$variable, levels = c("fric",levels(mFD$variable)[-5]) , ordered = TRUE)
m3fr<-reshape2::melt(top3_fr[,-c(1:2)] , id.vars = c("SF","subset","data","nb_samples","nb_env_samples", "notus") )
m3fdis<-reshape2::melt(top3_fdis[,-c(1:2)] , id.vars = c("SF","subset","data","nb_samples","nb_env_samples","notus") )

# plot across Rank abundance Curves
scaleFUN <- function(x) sprintf("%.2f", x)
unique(mFD$variable)
gfd<-ggplot(data = mFD[mFD$variable %in% c("fric","fdis"),], aes(x = subset, y = value, color = variable))+
  facet_grid(SF~., scale = "free_x", space = "free")+
  geom_line(size= 1)+
  geom_point(data = m3fr[m3fr$variable == "fric",], shape = 21, size = 3, stroke = 1)+
  geom_point(data = m3fdis[m3fdis$variable == "fdis",], shape = 21, size = 3, stroke = 1)+
  scale_color_manual(values = c("#7aadff","#063175", "#ff491c", "#cccccc"), labels = c("Richness","Dispersion", "Divergence"))+
  scale_x_continuous(expand =c(0,0), limits =c(0,75000), breaks = seq(0,80000, 10000))+
  ylab("Value [0,1]")+
  xlab("Subsets of Rank Abundance Curves (5000 OTUs)")+
  ggtitle("Functional Diversity across Rank Abundance Curves")+
  labs(color = "Functional Index: ")+
  theme(legend.position = "bottom")+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5, face = "bold"))+
  theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text =  element_text(size = 18, face = "bold", angle =0, hjust = 0))+
  theme(strip.text.y =  element_text(size = 18, face = "bold", angle =0, hjust = 0))+
  theme(strip.background = element_blank())+
  theme(legend.key =  element_rect(fill = "#00000000"))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.grid.major.x = element_blank())+
  theme(panel.grid.minor.x = element_blank())+
  theme(panel.grid.minor.y = element_line(color = "#ffffff", size = 0.5))+
  theme(panel.grid.major.y = element_line(color = "#ffffff", size = 1))+
  theme(panel.background = element_rect(fill = "#ebebfa",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0));gfd

# plot FDis vs nb of Taxonomic reference
head(FD)
scaleFUN <- function(x) sprintf("%.2f", x)
gfdis<-ggplot(data = FD, aes(x = sp_richn, y = fdis))+
  geom_point( shape = 19, size = 4, stroke = 1, alpha = 0.5, color = "#063175")+
  geom_text(data =top3_fdis, aes(label = substr(sub,1,1)), size = 2.5, fontface = "bold", color = "gray95")+
  scale_y_continuous(expand =c(0,0), limits=c(0.50,0.65),breaks = seq(0.55, 0.60, 0.05), labels = scaleFUN)+
  scale_x_continuous(expand =c(0,0), limits =c(340,600), breaks = seq(350,550, 100))+
  ylab("Functional\nDispersion")+
  xlab("# Taxonomic Reference\nin OTU Subset")+
  #theme(aspect.ratio = 1)+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5, face = "bold"))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text =  element_text(size = 18, face = "bold", angle =0, hjust = 0))+
  theme(strip.text.y =  element_text(size = 18, face = "bold", angle =0, hjust = 0))+
  theme(strip.background = element_blank())+
  theme(legend.key =  element_rect(fill = "#00000000"))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.grid.minor.x = element_line(color = "#ffffff", size = 0.5))+
  theme(panel.grid.major.x = element_line(color = "#ffffff", size = 1))+
  theme(panel.grid.minor.y = element_line(color = "#ffffff", size = 0.5))+
  theme(panel.grid.major.y = element_line(color = "#ffffff", size = 1))+
  theme(panel.background = element_rect(fill = "#ebebfa",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gfdis

# plot FRic vs nb of Taxonomic reference
head(FD)
scaleFUN <- function(x) sprintf("%.2f", x)
gfric<-ggplot(data = FD, aes(x = sp_richn, y = fric))+
  geom_point( shape = 19, size = 4, stroke = 1, alpha = 0.5, color = "#7aadff")+
  geom_text(data =top3_fr, aes(label = substr(sub,1,1)), size = 2.5, fontface = "bold", color = "gray5")+
  scale_y_continuous(expand =c(0,0), limits =c(0.55, 0.85), breaks = seq(0.60, 0.80, 0.1), labels = scaleFUN)+
  scale_x_continuous(expand =c(0,0), limits =c(340,600), breaks = seq(350,550, 100))+
  ylab("Functional\nRichness")+
  xlab("")+
  #theme(aspect.ratio = 1)+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5, face = "bold"))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text =  element_text(size = 18, face = "bold", angle =0, hjust = 0))+
  theme(strip.text.y =  element_text(size = 18, face = "bold", angle =0, hjust = 0))+
  theme(strip.background = element_blank())+
  theme(legend.key =  element_rect(fill = "#00000000"))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.grid.minor.x = element_line(color = "#ffffff", size = 0.5))+
  theme(panel.grid.major.x = element_line(color = "#ffffff", size = 1))+
  theme(panel.grid.minor.y = element_line(color = "#ffffff", size = 0.5))+
  theme(panel.grid.major.y = element_line(color = "#ffffff", size = 1))+
  theme(panel.background = element_rect(fill = "#ebebfa",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gfric

# 15 x 7.5
pp<-gfric/gfdis
gfd + pp +  
  plot_layout(design = "AAAB") +
  plot_annotation(tag_levels = 'A')  & 
  theme(plot.tag= element_text(size = 18, face = "bold"))

#####################
# END
#####################