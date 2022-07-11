#library(RMySQL)
library(readxl)
library(tidyverse)
library(Biostrings)
library(ggplot2)
library(zoo)
library(seqinr)
library(ggnewscale)
library(ggtree)
library("rjson")
library(ggmap)
library(rworldmap)
library(ape)
library(lubridate)
library(phytools)
library(maptools)
library(sf)
library(geojsonsf)
library(ggpubr)
library(treeio)
library(bdskytools)
library(circlize)
library("xlsx")
library(ggh4x)
library(phangorn)
library(seqinr)
library("ggalt")

setwd("/Users/ruopengxie/vjlab Dropbox/Ruopeng Xie/repos/fluduringCOVID-China/scripts")

#phylogeography

#3a1
tree <- read.beast("../results/MCC tree/China-flub-wholegenome.3a1.sub.early.afa.ca.mcc.tree")
(tip_location <- read_tsv("../results/MCC tree/3a1_sub_early_annotation.txt"))

node_location = c("North","East","Southwest","Northwest","Central","South")
location_color=c("#dde889","#81c77f","#d6ade5","#52a3cc","#f8da65","#fc8d62")

clade <- tree_subset(tree,node=436, levels_back=0)

(p <- ggtree(clade, mrsd = "2021-09-01",aes(color=region)) + 
    geom_range(range='height_0.95_HPD', size=1, color="#BABABA",alpha=0.3)+
    theme_tree2()+
    #geom_text(aes(label=node), hjust=-.3)+
    #geom_rootedge(rootedge = 0.01, color = "#cc5252")+            #root length, set root (Goldfields-Esperance) color
    scale_color_manual(values=location_color,breaks=node_location)+
    #geom_text2(aes(subset = !isTip, label=label),size = 1.5)+  #support value
    #geom_text(aes(label=round(as.numeric(location.prob), 2)))+
    scale_x_continuous(limits = c(2020.164,2021.7),breaks = c(2020.164,2020.331,2020.497,2020.667,2020.833,2021.003,2021.164,2021.331,2021.497,2021.667), labels = c("","","July","","","Jan","","","July",""))+
    scale_y_continuous(expand = c(0.02,0))
    #theme(axis.text.x = element_text(angle=25, hjust=1))
  #geom_hilight(node=1253, fill="steelblue", alpha=0.5) +
  #geom_hilight(node=1215, fill="#8852b8", alpha=0.5)
)

## tree
head(p$data$label)
p$data$tip_location <- sapply(p$data$label, function(x){
  tmp <- tip_location %>% filter(headers == x) %>% .$region
  if(length(tmp)==0){return(NA)}else{return(tmp)}
})

sort(table(p$data$tip_location))


mcc.tree.3a1 <- p + geom_tippoint(aes(fill = factor(tip_location)), size = 1.8, shape = 21, color="black",stroke = 0.4) + 
  scale_fill_manual(values=location_color,breaks=node_location)+
  #geom_nodelab(aes(label=round(as.numeric(region.prob), 2)),size=3,color="black",hjust=-0.1)+
  #scale_alpha_manual(values=c(1,1,1,1,0.7))+
  coord_cartesian(clip="off")+
  #geom_vline(xintercept = 2020.331,colour="#1b9e77")+
  #geom_vline(xintercept = 2021.331,colour="#1b9e77")+
  #geom_vline(xintercept = 2020.749,colour="#e7298a")+
  #geom_vline(xintercept = 2021.115,colour="red")+
  #theme(legend.position=c(0.85,0.75))+
  #guides(fill=guide_legend(title = "Location"),color="none")
  guides(fill=FALSE,color="none")
mcc.tree.3a1
ggsave("../results/3a1_mcc_tree.pdf",mcc.tree.3a1, width=9,height=20,units="cm")

#circlize
c.data.parent <- data.frame("node" = p$data$parent)
c.data <- data.frame("node" = p$data$node, "region_destination"= p$data$region,"parent_node" = p$data$parent, row.names = NULL)
c.data.parent <- merge(c.data.parent,c.data,by = "node", all.x = TRUE) %>% 
  subset(,c("node","region_destination")) %>% 
  rename(node="parent_node",region_destination="region_origin") %>% 
  unique()

circlize.data <- merge(c.data,c.data.parent,by = "parent_node", all.x = TRUE) %>% 
  subset(,c("region_origin","region_destination")) %>% 
  group_by(region_origin,region_destination) %>% 
  summarise(count = n()) %>%
  pivot_wider(id_cols="region_origin", names_from = "region_destination",values_from = count) %>%
  as.data.frame

circlize.data[is.na(circlize.data)] <- 0
rownames(circlize.data) <- circlize.data[,1]
circlize.data <- circlize.data[,-1]
circlize.data <- circlize.data[ , order(names(circlize.data))]

#write.table(circlize.data,"../analysis/whole_genome/Beast/3a1_transmission.tsv",sep="\t",row.names = T)

grid.col = c(North="#dde889",East="#81c77f",Southwest="#d6ade5",Northwest="#52a3cc",Central="#f8da65",South="#fc8d62")
circos.par(gap.after = c("North" = 5, "East" = 5, "Southwest" = 5, "Northwest" = 5, "Central" = 5,"South" = 5))

pdf(file="../results/circlize_3a1.pdf", width=3.8,height=3.8)
chordDiagram(data.matrix(circlize.data), grid.col = grid.col, annotationTrack = "grid",annotationTrackHeight = c(0.1),
             directional = 1, direction.type = c("arrows","diffHeight"),link.arr.type = "big.arrow",diffHeight = -mm_h(2))
circos.clear()
dev.off()

#3a2
tree <- read.beast("../results/MCC tree/China-flub-wholegenome.3a2.afa.ca.mcc.tree")
(tip_location <- read_tsv("../results/MCC tree/3a2_annotation.txt"))

node_location = c("North","East","Southwest","Northwest","Central","South")
location_color=c("#dde889","#81c77f","#d6ade5","#52a3cc","#f8da65","#fc8d62")

(p <- ggtree(tree, mrsd = "2021-09-02",aes(color=region)) + 
    geom_range(range='height_0.95_HPD', size=1, color="#BABABA",alpha=0.3)+
    theme_tree2()+
    #geom_rootedge(rootedge = 0.01, color = "#cc5252")+            #root length, set root (Goldfields-Esperance) color
    scale_color_manual(values=location_color,breaks=node_location)+
    #geom_text2(aes(subset = !isTip, label=label),size = 1.5)+  #support value
    #geom_text(aes(label=round(as.numeric(location.prob), 2)))+
    #scale_x_continuous(limits = c(2019.4,2021.7),breaks = c(2019.497,2019.667,2019.833,2020.003,2020.164,2020.331,2020.497,2020.667,2020.833,2021.003,2021.164,2021.331,2021.497,2021.667), labels = c("Jul","","","Jan","","","July","","","Jan","","","July",""))+
    scale_x_continuous(limits = c(2020.63,2021.7),breaks = c(2020.667,2020.915,2021.164,2021.415,2021.667), labels = c("Sep-2020","","Mar-2021","","Sep-2021"))+
    scale_y_continuous(expand = c(0.02,0))
  #geom_hilight(node=1253, fill="steelblue", alpha=0.5) +
  #geom_hilight(node=1215, fill="#8852b8", alpha=0.5)
)

## tree
head(p$data$label)
p$data$tip_location <- sapply(p$data$label, function(x){
  tmp <- tip_location %>% filter(headers == x) %>% .$region
  if(length(tmp)==0){return(NA)}else{return(tmp)}
})

sort(table(p$data$tip_location))


mcc.tree.3a2 <- p + geom_tippoint(aes(fill = factor(tip_location)), size = 2, shape = 21, color="black",stroke = 0.4) + 
  scale_fill_manual(values=location_color,breaks=node_location)+
  #geom_nodelab(aes(label=round(as.numeric(region.prob), 2)),size=3,color="black",hjust=-0.1)+
  #scale_alpha_manual(values=c(1,1,1,1,0.7))+
  coord_cartesian(clip="off")+
  #geom_vline(xintercept = 2021.331,colour="#1b9e77")+
  #geom_vline(xintercept = 2020.749,colour="#e7298a")+
  #geom_vline(xintercept = 2021.115,colour="red")+
  #theme(legend.position=c(0.85,0.75))+
  #guides(fill=guide_legend(title = "Location"),color="none")
  guides(fill=FALSE,color="none")

mcc.tree.3a2
ggsave("../results/3a2_mcc_tree.pdf",mcc.tree.3a2, width=9,height=14,units="cm")

#circlize
c.data.parent <- data.frame("node" = p$data$parent)
c.data <- data.frame("node" = p$data$node, "region_destination"= p$data$region,"parent_node" = p$data$parent, row.names = NULL)
c.data.parent <- merge(c.data.parent,c.data,by = "node", all.x = TRUE) %>% 
  subset(,c("node","region_destination")) %>% 
  rename(node="parent_node",region_destination="region_origin") %>% 
  unique()

circlize.data <- merge(c.data,c.data.parent,by = "parent_node", all.x = TRUE) %>% 
  subset(,c("region_origin","region_destination")) %>% 
  group_by(region_origin,region_destination) %>% 
  summarise(count = n()) %>%
  pivot_wider(id_cols="region_origin", names_from = "region_destination",values_from = count) %>%
  as.data.frame

circlize.data[is.na(circlize.data)] <- 0
rownames(circlize.data) <- circlize.data[,1]
circlize.data <- circlize.data[,-1]
circlize.data <- circlize.data[ , order(names(circlize.data))]

#write.table(circlize.data,"../analysis/whole_genome/Beast/3a2_transmission.tsv",sep="\t",row.names = T)

grid.col = c(North="#dde889",East="#81c77f",Southwest="#d6ade5",Northwest="#52a3cc",Central="#f8da65",South="#fc8d62")
circos.par(gap.after = c("North" = 5, "East" = 5, "Southwest" = 5, "Northwest" = 5, "Central" = 5,"South" = 5))

pdf(file="../results/circlize_3a2.pdf", width=3.8,height=3.8)
chordDiagram(data.matrix(circlize.data), grid.col = grid.col, annotationTrack = "grid",annotationTrackHeight = c(0.1),
             directional = 1, direction.type = c("arrows","diffHeight"),link.arr.type = "big.arrow",diffHeight = -mm_h(2))

circos.clear()
dev.off()

#V1A
tree <- read.beast("../results/MCC tree/China-flub-wholegenome.V1A.afa.ca.mcc.tree")
(tip_location <- read_tsv("../results/MCC tree/V1A_annotation.txt"))

node_location = c("North","East","Southwest","Northwest","Central","South")
location_color=c("#dde889","#81c77f","#d6ade5","#52a3cc","#f8da65","#fc8d62")

(p <- ggtree(tree, mrsd = "2021-03-15",aes(color=region)) + 
    geom_range(range='height_0.95_HPD', size=1, color="#BABABA",alpha=0.3)+
    theme_tree2()+
    #geom_rootedge(rootedge = 0.01, color = "#cc5252")+            #root length, set root (Goldfields-Esperance) color
    scale_color_manual(values=location_color,breaks=node_location)+
    #geom_text2(aes(subset = !isTip, label=label),size = 1.5)+  #support value
    #geom_text(aes(label=round(as.numeric(location.prob), 2)))+
    #scale_x_continuous(limits = c(2020.6,2021.8),breaks = c(2020.6,2020.9,2021.2,2021.5,2021.8))+
    scale_y_continuous(expand = c(0.02,0))
  #geom_hilight(node=1253, fill="steelblue", alpha=0.5) +
  #geom_hilight(node=1215, fill="#8852b8", alpha=0.5)
)

## tree
head(p$data$label)
p$data$tip_location <- sapply(p$data$label, function(x){
  tmp <- tip_location %>% filter(headers == x) %>% .$region
  if(length(tmp)==0){return(NA)}else{return(tmp)}
})

sort(table(p$data$tip_location))


mcc.tree.V1A <- p + geom_tippoint(aes(fill = factor(tip_location)), size = 1.8, shape = 21, color="black",stroke = 0.4) + 
  scale_fill_manual(values=location_color,breaks=node_location)+
  #geom_nodelab(aes(label=round(as.numeric(region.prob), 2)),size=3,color="black",hjust=-0.1)+
  #scale_alpha_manual(values=c(1,1,1,1,0.7))+
  coord_cartesian(clip="off")+
  #geom_vline(xintercept = 2018.128,colour="red", linetype = "longdash")+
  #geom_vline(xintercept = 2018.331,colour="red", linetype = "longdash")+
  #geom_vline(xintercept = 2018.749,colour="red", linetype = "longdash")+
  #geom_vline(xintercept = 2019.098,colour="red", linetype = "longdash")+
  #geom_vline(xintercept = 2019.331,colour="red", linetype = "longdash")+
  #geom_vline(xintercept = 2019.749,colour="red", linetype = "longdash")+
  #geom_vline(xintercept = 2020.068,colour="red", linetype = "longdash")+
  #theme(legend.position=c(0.85,0.75))+
  #guides(fill=guide_legend(title = "Location"),color="none")
  guides(fill=FALSE,color="none")

mcc.tree.V1A
ggsave("../results/V1A_mcc_tree.pdf",mcc.tree.V1A, width=9,height=12,units="cm")

#circlize
c.data.parent <- data.frame("node" = p$data$parent)
c.data <- data.frame("node" = p$data$node, "region_destination"= p$data$region,"parent_node" = p$data$parent, row.names = NULL)
c.data.parent <- merge(c.data.parent,c.data,by = "node", all.x = TRUE) %>% 
  subset(,c("node","region_destination")) %>% 
  rename(node="parent_node",region_destination="region_origin") %>% 
  unique()

circlize.data <- merge(c.data,c.data.parent,by = "parent_node", all.x = TRUE) %>% 
  subset(,c("region_origin","region_destination")) %>% 
  group_by(region_origin,region_destination) %>% 
  summarise(count = n()) %>%
  pivot_wider(id_cols="region_origin", names_from = "region_destination",values_from = count) %>%
  as.data.frame

circlize.data[is.na(circlize.data)] <- 0
rownames(circlize.data) <- circlize.data[,1]
circlize.data <- circlize.data[,-1]
circlize.data <- circlize.data[ , order(names(circlize.data))]

#write.table(circlize.data,"../analysis/whole_genome/Beast/V1A_transmission.tsv",sep="\t",row.names = T)

grid.col = c(North="#dde889",East="#81c77f",Southwest="#d6ade5",Northwest="#52a3cc",Central="#f8da65",South="#fc8d62")
circos.par(gap.after = c("North" = 5, "East" = 5, "Southwest" = 5, "Northwest" = 5, "Central" = 5,"South" = 5))

pdf(file="../results/circlize_V1A.pdf", width=3.8,height=3.8)

chordDiagram(data.matrix(circlize.data), grid.col = grid.col, annotationTrack = "grid",annotationTrackHeight = c(0.1),
             directional = 1, direction.type = c("arrows","diffHeight"),link.arr.type = "big.arrow",diffHeight = -mm_h(2))

circos.clear()
dev.off()
