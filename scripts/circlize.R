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
ggsave("../figures/V1A_mcc_tree.pdf",mcc.tree.V1A, width=9,height=12,units="cm")

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

pdf(file="../figures/circlize_V1A.pdf", width=3.8,height=3.8)

chordDiagram(data.matrix(circlize.data), grid.col = grid.col, annotationTrack = "grid",annotationTrackHeight = c(0.1),
             directional = 1, direction.type = c("arrows","diffHeight"),link.arr.type = "big.arrow",diffHeight = -mm_h(2))

circos.clear()
dev.off()


#transmission pattern in the map
#ref: https://www.codenong.com/js19935c44cbf7/
China_map = read_sf("../analysis/whole_genome/Beast/China.geojson.json")

China_map$region <- gsub(tolower("shanxi2|Gansu|Qinghai|Ningxia|Xinjiang"),"Northwest",China_map$drilldown)
China_map$region <- gsub(tolower("Beijing|Tianjin|Hebei|Shanxi|Neimenggu"),"North",China_map$region)
China_map$region <- gsub(tolower("Liaoning|Jilin|Heilongjiang"),"North",China_map$region)
China_map$region <- gsub(tolower("Shanghai|Jiangsu|Zhejiang|Anhui|Fujian|Jiangxi|Shandong"),"East",China_map$region)
China_map$region <- gsub(tolower("Henan|Hubei|Hunan"),"Central",China_map$region)
China_map$region <- gsub(tolower("Guangdong|Guangxi|Hainan"),"South",China_map$region)
China_map$region <- gsub(tolower("Chongqing|Sichuan|Guizhou|Yunnan|Xizang"),"Southwest",China_map$region)
China_map$region <- gsub("xianggang|aomen","Other",China_map$region)
China_map$region[is.na(China_map$region)] <- "Other"

unique(China_map$region)

g.china.map <- ggplot(China_map)+
  geom_sf(aes(fill = region),size=0.3) +
  scale_fill_manual(values=c(location_color,"#BABABA"),breaks=c(node_location,"Other"))+
  #geom_curve(aes(x = 3, y = 3, xend = 10, yend = 10),arrow = arrow(length = unit(0.03, "npc"), type="closed"),colour = "#EC7014",size = 1.2,angle = 90)+
  #coord_sf(xlim = c(80, 130),ylim = c(20, 50))+
  labs(title="",x="",y="")+
  theme(panel.grid =element_blank(),panel.background = element_blank(),legend.position="none") +   ## 删去网格线
  theme(axis.text = element_blank()) +   ## 删去刻度标签
  theme(axis.ticks = element_blank()) +   ## 删去刻度线
  theme(panel.border = element_blank())+   ## 删去外层边框
  theme(plot.margin = margin(t = 0,  # 顶部边缘距离
                             r = 0,  # 右边边缘距离
                             b = 0,  # 底部边缘距离
                             l = 0,  # 左边边缘距离
                             unit = "cm"))
g.china.map

ggsave("../figures/china_map.pdf",g.china.map, width=9,height=12,units="cm")


#BDSKY
lf.3a1 <- readLogfile("../analysis/whole_genome/Beast2-fixed/China-flub-wholegenome-3a1-early.log", burnin=0)

Re_sky.3a1    <- getSkylineSubset(lf.3a1, "reproductiveNumber")
Re_hpd.3a1   <- getMatrixHPD(Re_sky.3a1)
delta_hpd.3a1 <- getHPD(lf.3a1$becomeUninfectiousRate)

plotSkyline(1:10, Re_hpd.3a1, type='step', ylab="R")

origin_max <- mean(lf.3a1$origin_BDSKY_Serial)

timegrid.3a1 <- seq(0,origin_max,length.out=101)
Re_gridded.3a1 <- gridSkyline(Re_sky.3a1, lf.3a1$origin, timegrid.3a1)
Re_gridded_hpd.3a1 <- getMatrixHPD(Re_gridded.3a1)

df.bdsky.3a1 <- as.data.frame(t(as.data.frame(Re_gridded_hpd.3a1))) 

times.3a1 <- 2021.667-timegrid.3a1

rownames(df.bdsky.3a1) <- NULL
colnames(df.bdsky.3a1) <- c("Lower","Value","Upper")

df.bdsky.3a1$Time <- times.3a1
df.bdsky.3a1$clade <- "3a1"

df.bdsky.3a1 <- df.bdsky.3a1 %>% filter(Time > 2020.273)


#3a2
lf.3a2 <- readLogfile("../analysis/whole_genome/Beast2-fixed/China-flub-wholegenome-3a2.log", burnin=0)

Re_sky.3a2    <- getSkylineSubset(lf.3a2, "reproductiveNumber")
Re_hpd.3a2   <- getMatrixHPD(Re_sky.3a2)
delta_hpd.3a2 <- getHPD(lf.3a2$becomeUninfectiousRate)

plotSkyline(1:10, Re_hpd.3a2, type='step', ylab="R")

origin_max <- mean(lf.3a2$origin_BDSKY_Serial)

timegrid.3a2 <- seq(0,origin_max,length.out=101)
Re_gridded.3a2 <- gridSkyline(Re_sky.3a2, lf.3a2$origin, timegrid.3a2)
Re_gridded_hpd.3a2 <- getMatrixHPD(Re_gridded.3a2)

df.bdsky.3a2 <- as.data.frame(t(as.data.frame(Re_gridded_hpd.3a2)))

times.3a2 <- 2021.669-timegrid.3a2

rownames(df.bdsky.3a2) <- NULL
colnames(df.bdsky.3a2) <- c("Lower","Value","Upper")

df.bdsky.3a2$Time <- times.3a2
df.bdsky.3a2$clade <- "3a2"

#V1A
lf.V1A <- readLogfile("../analysis/whole_genome/Beast2/China-flub-wholegenome-V1A.log", burnin=0)

Re_sky.V1A    <- getSkylineSubset(lf.V1A, "reproductiveNumber")
Re_hpd.V1A   <- getMatrixHPD(Re_sky.V1A)
delta_hpd.V1A <- getHPD(lf.V1A$becomeUninfectiousRate)

plotSkyline(1:10, Re_hpd.V1A, type='step', ylab="R")

TreeHeight<- 2.436
origin_max <- 9.31
info_span <- 2.202
timegrid.V1A <- seq(0,info_span,length.out=101)
Re_gridded.V1A <- gridSkyline(Re_sky.V1A, lf.V1A$origin, timegrid.V1A)
Re_gridded_hpd.V1A <- getMatrixHPD(Re_gridded.V1A)

df.bdsky.V1A <- as.data.frame(t(as.data.frame(Re_gridded_hpd.V1A)))

times.V1A <- 2021.202-timegrid.V1A

rownames(df.bdsky.V1A) <- NULL
colnames(df.bdsky.V1A) <- c("Lower","Value","Upper")

df.bdsky.V1A$Time <- times.V1A
df.bdsky.V1A$clade <- "V1A"

#plot
df.bdsky <- rbind(df.bdsky.3a1,df.bdsky.3a2,df.bdsky.V1A)

flu.bdsky <- ggplot(df.bdsky)+
  geom_line(aes(x=as.Date(date_decimal(Time,tz = "EST")),y=Value,color=clade),)+
  geom_ribbon(aes(x=as.Date(date_decimal(Time,tz = "EST")),ymin=Lower,ymax=Upper,fill=clade), alpha=0.3) +
  scale_fill_manual(values=c("#40ab5d","#e06c9a","#7f7dbb"))+
  scale_color_manual(values=c("#40ab5d","#e06c9a","#7f7dbb"))+
  geom_hline(aes(yintercept=1),linetype = "dashed",color="grey")+
  labs(x = "Date", y = "Re")+
  #scale_fill_manual(values=c("#9970ab","#d7191c","#008837","#313695"))+
  scale_y_continuous(expand=c(0,0),breaks=c(0,2,4,6,8,10,12), limits = c(0,12),guide = "axis_minor")+
  scale_x_date(breaks="4 months", limits = c(as.Date("2019-01-01"),as.Date("2021-11-01")),date_labels = "%b",guide = "axis_minor")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
        legend.position=c(0.85,0.6))

flu.bdsky

ggsave("./bdskys_skyline.pdf",flu.bdsky,width=10,height=6,units="cm")


#seq distribution in geographic level over time
gisaid.vic.china.meta$province <- gsub("-.*|B/","",gisaid.vic.china.meta$Isolate_Name)
gisaid.vic.china.meta$province <- gsub("BVR","Sichuan",gisaid.vic.china.meta$province)
gisaid.vic.china.meta$province <- gsub("/.*","",gisaid.vic.china.meta$province)
gisaid.vic.china.meta$province <- gsub("Nanchang*","Jiangxi",gisaid.vic.china.meta$province)


gisaid.vic.china.meta$region <- gsub("Shaanxi|Gansu|Qinghai|Ningxia|Xinjiang","Northwest",gisaid.vic.china.meta$province)
gisaid.vic.china.meta$region <- gsub("Beijing|Tianjin|Hebei|Shanxi|Neimenggu","North",gisaid.vic.china.meta$region)
gisaid.vic.china.meta$region <- gsub("Liaoning|Jilin|Heilongjiang","North",gisaid.vic.china.meta$region)
gisaid.vic.china.meta$region <- gsub("Shanghai|Jiangsu|Zhejiang|Anhui|Fujian|Jiangxi|Shandong","East",gisaid.vic.china.meta$region)
gisaid.vic.china.meta$region <- gsub("Henan|Hubei|Hunan","Central",gisaid.vic.china.meta$region)
gisaid.vic.china.meta$region <- gsub("Guangdong|Guangxi|Hainan|Shenzhen","South",gisaid.vic.china.meta$region)
gisaid.vic.china.meta$region <- gsub("Chongqing|Sichuan|Guizhou|Yunnan|Xizang","Southwest",gisaid.vic.china.meta$region)
#gisaid.vic.china.meta$region <- gsub("xianggang|aomen","Other",gisaid.vic.china.meta$region)
#gisaid.vic.china.meta$region[is.na(gisaid.vic.china.meta$region)] <- "Other"

table(gisaid.vic.china.meta$region)

gisaid.vic.china.region.2019 <- gisaid.vic.china.meta %>%
  filter(Collection_Date > "2018-12-31" & region %in% c("Northwest","North","East","Central","South","Southwest")) %>%
  subset(,c("Collection_Date","province","region")) %>%
  group_by(region, province, Collection_Year=floor_date(as.Date(Collection_Date),"month")) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = "Collection_Year",values_from = count)

gisaid.vic.china.region.2019[is.na(gisaid.vic.china.region.2019)] <- 0

write_tsv(gisaid.vic.china.region.2019,"../figures/China_seq_geography.tsv")


#phylogeograph based on balance dataset
china.3a1.clade <- read.csv("../analysis/whole_genome/Beast_balance/3a1.b.headers.txt")
china.3a2.clade <- read.csv("../analysis/whole_genome/Beast_balance/3a2.b.headers.txt")

#china.3a1.clade$date <- gsub(".*\\|","",china.3a1.clade$headers)
#china.3a2.clade$date <- gsub(".*\\|","",china.3a2.clade$headers)

china.3a1.clade$province <- gsub("-.*|.*\\|B/","",china.3a1.clade$headers)
china.3a1.clade$province <- gsub("EPI_ISL_4197189\\|BVR","Sichuan",china.3a1.clade$province)
china.3a1.clade$province <- gsub("Guizhou/13587/2020\\|2020","Guizhou",china.3a1.clade$province)

china.3a2.clade$province <- gsub("-.*|.*\\|B/","",china.3a2.clade$headers)

china.3a1.clade$region <- gsub("Beijing|Tianjin|Hebei|Shanxi|Neimenggu","North",china.3a1.clade$province)
china.3a1.clade$region <- gsub("Liaoning|Jilin|Heilongjiang","North",china.3a1.clade$region)
china.3a1.clade$region <- gsub("Shanghai|Jiangsu|Zhejiang|Anhui|Fujian|Jiangxi|Shandong","East",china.3a1.clade$region)
china.3a1.clade$region <- gsub("Henan|Hubei|Hunan","Central",china.3a1.clade$region)
china.3a1.clade$region <- gsub("Guangdong|Guangxi|Hainan","South",china.3a1.clade$region)
china.3a1.clade$region <- gsub("Chongqing|Sichuan|Guizhou|Yunnan|Xizang","Southwest",china.3a1.clade$region)
china.3a1.clade$region <- gsub("Shaanxi|Gansu|Qinghai|Ningxia|Xinjiang","Northwest",china.3a1.clade$region)

china.3a2.clade$region <- gsub("Beijing|Tianjin|Hebei|Shanxi|Neimenggu","North",china.3a2.clade$province)
china.3a2.clade$region <- gsub("Liaoning|Jilin|Heilongjiang","North",china.3a2.clade$region)
china.3a2.clade$region <- gsub("Shanghai|Jiangsu|Zhejiang|Anhui|Fujian|Jiangxi|Shandong","East",china.3a2.clade$region)
china.3a2.clade$region <- gsub("Henan|Hubei|Hunan","Central",china.3a2.clade$region)
china.3a2.clade$region <- gsub("Guangdong|Guangxi|Hainan","South",china.3a2.clade$region)
china.3a2.clade$region <- gsub("Chongqing|Sichuan|Guizhou|Yunnan|Xizang","Southwest",china.3a2.clade$region)
china.3a2.clade$region <- gsub("Shaanxi|Gansu|Qinghai|Ningxia|Xinjiang","Northwest",china.3a2.clade$region)

table(china.3a1.clade$region)
table(china.3a2.clade$region)

table(china.3a1.clade$province)
table(china.3a2.clade$province)


write_tsv(china.3a1.clade,"../analysis/whole_genome/Beast_balance/3a1_annotation.txt")
write_tsv(china.3a2.clade,"../analysis/whole_genome/Beast_balance/3a2_annotation.txt")


#seqs distribution over provinces
province.gps <- read_tsv("../analysis/iqtree/province_gps.txt")

china.3a1.clade <- read_tsv("../analysis/whole_genome/Beast/3a1_annotation.txt") %>% 
  subset(,c("province")) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename(.,c("."="province","Freq"="count")) %>% 
  merge(.,province.gps,by = "province", all.x = TRUE)

china.3a1.b.clade <- read_tsv("../analysis/whole_genome/Beast_balance/3a1_annotation.txt") %>% 
  subset(,c("province")) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename(.,c("."="province","Freq"="count")) %>% 
  merge(.,province.gps,by = "province", all.x = TRUE)

china.3a2.clade <- read_tsv("../analysis/whole_genome/Beast/3a2_annotation.txt") %>% 
  subset(,c("province")) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename(.,c("."="province","Freq"="count")) %>% 
  merge(.,province.gps,by = "province", all.x = TRUE)

china.3a2.b.clade <- read_tsv("../analysis/whole_genome/Beast_balance/3a2_annotation.txt") %>% 
  subset(,c("province")) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename(.,c("."="province","Freq"="count")) %>% 
  merge(.,province.gps,by = "province", all.x = TRUE)

china.V1A.clade <- read_tsv("../analysis/whole_genome/Beast/V1A_annotation.txt") %>% 
  subset(,c("province")) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename(.,c("."="province","Freq"="count")) %>% 
  merge(.,province.gps,by = "province", all.x = TRUE)


API_pre = "http://xzqh.mca.gov.cn/data/"
China_map = st_read(dsn = paste0(API_pre, "quanguo.json"),
                    stringsAsFactors=FALSE)
st_crs(China_map) = 4326

province.map.3a1 <- ggplot(China_map)+
  geom_sf()+
  geom_point(aes(x=longitude, y=latitude,size=count),data=china.3a1.clade,shape=21,color="black",fill="#7ea34f")+
  scale_size_area(breaks = c(10,20,50,100,150),limits=c(1,180),max_size = 10)+
  #geom_scatterpie(aes(x=longitude, y=latitude, group=province,r=all_seqs/5),data=gisaid.vic.china.meta.filter, cols=c("V1A","3a1","3a2"), size=0.1,color="black", alpha = 0.8)+ #add pie chart
  #scale_fill_manual(values=c("#997ba5","#7ea34f","#50b2b2"),breaks=c("V1A","3a1","3a2"))+
  theme(panel.grid =element_blank(),panel.background = element_blank(),legend.position="none") +   ## 删去网格线
  theme(axis.text = element_blank()) +   ## 删去刻度标签
  theme(axis.ticks = element_blank()) +   ## 删去刻度线
  theme(panel.border = element_blank())+   ## 删去外层边框
  theme(plot.margin = margin(t = 0,  # 顶部边缘距离
                             r = 0,  # 右边边缘距离
                             b = 0,  # 底部边缘距离
                             l = 0,  # 左边边缘距离
                             unit = "cm"))+
  labs(x = "", y = "", title ="3a1")

province.map.3a1

province.map.3a1.b <- ggplot(China_map)+
  geom_sf()+
  geom_point(aes(x=longitude, y=latitude,size=count),data=china.3a1.b.clade,shape=21,color="black",fill="#7ea34f")+
  #scale_size_continuous(breaks = c(10,20,50,100,150),range=c(1,6))+
  scale_size_area(breaks = c(10,20,50,100,150),limits=c(1,180),max_size = 10)+
  #geom_scatterpie(aes(x=longitude, y=latitude, group=province,r=all_seqs/5),data=gisaid.vic.china.meta.filter, cols=c("V1A","3a1.b","3a2"), size=0.1,color="black", alpha = 0.8)+ #add pie chart
  #scale_fill_manual(values=c("#997ba5","#7ea34f","#50b2b2"),breaks=c("V1A","3a1.b","3a2"))+
  theme(panel.grid =element_blank(),panel.background = element_blank(),legend.position="none") +   ## 删去网格线
  theme(axis.text = element_blank()) +   ## 删去刻度标签
  theme(axis.ticks = element_blank()) +   ## 删去刻度线
  theme(panel.border = element_blank())+   ## 删去外层边框
  theme(plot.margin = margin(t = 0,  # 顶部边缘距离
                             r = 0,  # 右边边缘距离
                             b = 0,  # 底部边缘距离
                             l = 0,  # 左边边缘距离
                             unit = "cm"))+
  labs(x = "", y = "", title ="3a1-balance")

province.map.3a1.b

province.map.3a2 <- ggplot(China_map)+
  geom_sf()+
  geom_point(aes(x=longitude, y=latitude,size=count),data=china.3a2.clade,shape=21,color="black",fill="#4abac0")+
  #scale_size_continuous(breaks = c(5,10,20,40,50),range=c(1,3))+
  scale_size_area(breaks = c(10,20,50,100,150),limits=c(1,180),max_size = 10)+
  #geom_scatterpie(aes(x=longitude, y=latitude, group=province,r=all_seqs/5),data=gisaid.vic.china.meta.filter, cols=c("V1A","3a2","3a2"), size=0.1,color="black", alpha = 0.8)+ #add pie chart
  #scale_fill_manual(values=c("#997ba5","#7ea34f","#50b2b2"),breaks=c("V1A","3a2","3a2"))+
  theme(panel.grid =element_blank(),panel.background = element_blank(),legend.position="none") +   ## 删去网格线
  theme(axis.text = element_blank()) +   ## 删去刻度标签
  theme(axis.ticks = element_blank()) +   ## 删去刻度线
  theme(panel.border = element_blank())+   ## 删去外层边框
  theme(plot.margin = margin(t = 0,  # 顶部边缘距离
                             r = 0,  # 右边边缘距离
                             b = 0,  # 底部边缘距离
                             l = 0,  # 左边边缘距离
                             unit = "cm"))+
  labs(x = "", y = "", title ="3a2")

province.map.3a2

province.map.3a2.b <- ggplot(China_map)+
  geom_sf()+
  geom_point(aes(x=longitude, y=latitude,size=count),data=china.3a2.b.clade,shape=21,color="black",fill="#4abac0")+
  #scale_size_continuous(breaks = c(10,20,50,100,150),range=c(1,6))+
  scale_size_area(breaks = c(10,20,50,100,150),limits=c(1,180),max_size = 10)+
  #geom_scatterpie(aes(x=longitude, y=latitude, group=province,r=all_seqs/5),data=gisaid.vic.china.meta.filter, cols=c("V1A","3a2.b","3a2"), size=0.1,color="black", alpha = 0.8)+ #add pie chart
  #scale_fill_manual(values=c("#997ba5","#7ea34f","#50b2b2"),breaks=c("V1A","3a2.b","3a2"))+
  theme(panel.grid =element_blank(),panel.background = element_blank(),legend.position=c(1.2,0.4)) +   ## 删去网格线
  theme(axis.text = element_blank()) +   ## 删去刻度标签
  theme(axis.ticks = element_blank()) +   ## 删去刻度线
  theme(panel.border = element_blank())+   ## 删去外层边框
  theme(plot.margin = margin(t = 0,  # 顶部边缘距离
                             r = 0,  # 右边边缘距离
                             b = 0,  # 底部边缘距离
                             l = 0,  # 左边边缘距离
                             unit = "cm"))+
  labs(x = "", y = "", title ="3a2-balance")

province.map.3a2.b

province.map.V1A <- ggplot(China_map)+
  geom_sf()+
  geom_point(aes(x=longitude, y=latitude,size=count),data=china.V1A.clade,shape=21,color="black",fill="#987bab")+
  #scale_size_continuous(breaks = c(5,10,15), range=c(1,2))+
  scale_size_area(breaks = c(10,20,50,100,150),limits=c(1,180),max_size = 10)+
  #geom_scatterpie(aes(x=longitude, y=latitude, group=province,r=all_seqs/5),data=gisaid.vic.china.meta.filter, cols=c("V1A","V1A","V1A"), size=0.1,color="black", alpha = 0.8)+ #add pie chart
  #scale_fill_manual(values=c("#997ba5","#7ea34f","#50b2b2"),breaks=c("V1A","V1A","V1A"))+
  theme(panel.grid =element_blank(),panel.background = element_blank(),legend.position="none") +   ## 删去网格线
  theme(axis.text = element_blank()) +   ## 删去刻度标签
  theme(axis.ticks = element_blank()) +   ## 删去刻度线
  theme(panel.border = element_blank())+   ## 删去外层边框
  theme(plot.margin = margin(t = 0,  # 顶部边缘距离
                             r = 0,  # 右边边缘距离
                             b = 0,  # 底部边缘距离
                             l = 0,  # 左边边缘距离
                             unit = "cm"))+
  labs(x = "", y = "", title ="V1A.3-China")

province.map.V1A

# province.map <- ggplot(China_map)+
#   geom_sf()+
#   geom_text(data=province.gps,aes(x=longitude, y=latitude, label=province),position = "identity",size=2)+
#   theme(panel.grid =element_blank(),panel.background = element_blank(),legend.position="none") +   ## 删去网格线
#   theme(axis.text = element_blank()) +   ## 删去刻度标签
#   theme(axis.ticks = element_blank()) +   ## 删去刻度线
#   theme(panel.border = element_blank())+   ## 删去外层边框
#   theme(plot.margin = margin(t = 0,  # 顶部边缘距离
#                              r = 0,  # 右边边缘距离
#                              b = 0,  # 底部边缘距离
#                              l = 0,  # 左边边缘距离
#                              unit = "cm"))+
#   labs(x = "", y = "", title ="")
# 
# province.map

ggsave("../figures/3a1_seqs_distribution_province.pdf",province.map.3a1,width=20,height=11,units="cm")
ggsave("../figures/3a1b_seqs_distribution_province.pdf",province.map.3a1.b,width=20,height=11,units="cm")
ggsave("../figures/3a2_seqs_distribution_province.pdf",province.map.3a2,width=20,height=11,units="cm")
ggsave("../figures/3a2b_seqs_distribution_province.pdf",province.map.3a2.b,width=20,height=11,units="cm")
ggsave("../figures/V1A_seqs_distribution_province.pdf",province.map.V1A,width=20,height=11,units="cm")
#ggsave("../figures/map_province.pdf",province.map,width=20,height=11,units="cm")
#ggarrange(province.map.3a1,province.map.3a2,province.map.V1A,province.map.3a1.b,province.map.3a2.b, ncol = 3,nrow=2)


#BSP
#3a1
bsp.3a1 <- read_tsv("../analysis/whole_genome/Beast/3a1_early_bayesian_skyline_data") %>% 
  filter(Time > 2020.4)

p.bsp.3a1 <- ggplot(bsp.3a1)+
  geom_xspline(aes(x=Time,y=Median),size=1,color="#1F6A9E")+
  geom_point(aes(x=Time,y=Median),size=0,color="#1F6A9E")+
  geom_xspline(aes(x=Time,y=Upper),size=0.5,color="#1F6A9E")+
  geom_point(aes(x=Time,y=Upper),size=0,color="#1F6A9E")+
  geom_xspline(aes(x=Time,y=Lower),size=0.5,color="#1F6A9E")+
  geom_point(aes(x=Time,y=Lower),size=0,color="#1F6A9E")+
  #geom_vline(xintercept = 2020.249,colour="#1b9e77")+
  theme_bw()+
  #ylim(-1,5.5)+
  theme(panel.grid =element_blank(),panel.background = element_blank(),legend.position="none") +   ## 删去网格线
  theme(plot.margin = margin(t = 0,  # 顶部边缘距离
                             r = 0,  # 右边边缘距离
                             b = 0,  # 底部边缘距离
                             l = 0,  # 左边边缘距离
                             unit = "cm"))
p.bsp.3a1
#export pdf 

#3a2
bsp.3a2 <- read_tsv("../analysis/whole_genome/Beast/3a2_bayesian_skyline_data")

p.bsp.3a2 <- ggplot(bsp.3a2)+
  geom_xspline(aes(x=Time,y=Median),size=1,color="#1F6A9E")+
  geom_point(aes(x=Time,y=Median),size=0,color="#1F6A9E")+
  geom_xspline(aes(x=Time,y=Upper),size=0.5,color="#1F6A9E")+
  geom_point(aes(x=Time,y=Upper),size=0,color="#1F6A9E")+
  geom_xspline(aes(x=Time,y=Lower),size=0.5,color="#1F6A9E")+
  geom_point(aes(x=Time,y=Lower),size=0,color="#1F6A9E")+
  geom_vline(xintercept = 2020.749,colour="#1b9e77")+
  ylim(-1,5.5)+
  theme_bw()+
  theme(panel.grid =element_blank(),panel.background = element_blank(),legend.position="none") +   ## 删去网格线
  theme(plot.margin = margin(t = 0,  # 顶部边缘距离
                             r = 0,  # 右边边缘距离
                             b = 0,  # 底部边缘距离
                             l = 0,  # 左边边缘距离
                             unit = "cm"))
p.bsp.3a2
#export pdf 

#V1A
bsp.V1A <- read_tsv("../analysis/whole_genome/Beast/V1A_bayesian_skyline_data")

p.bsp.V1A <- ggplot(bsp.V1A)+
  geom_xspline(aes(x=Time,y=Median),size=1,color="#1F6A9E")+
  geom_point(aes(x=Time,y=Median),size=0,color="#1F6A9E")+
  geom_xspline(aes(x=Time,y=Upper),size=0.5,color="#1F6A9E")+
  geom_point(aes(x=Time,y=Upper),size=0,color="#1F6A9E")+
  geom_xspline(aes(x=Time,y=Lower),size=0.5,color="#1F6A9E")+
  geom_point(aes(x=Time,y=Lower),size=0,color="#1F6A9E")+
  #geom_vline(xintercept = 2020.749,colour="#1b9e77")+
  #ylim(-1,8)+
  theme_bw()+
  theme(panel.grid =element_blank(),panel.background = element_blank(),legend.position="none") +   ## 删去网格线
  theme(plot.margin = margin(t = 0,  # 顶部边缘距离
                             r = 0,  # 右边边缘距离
                             b = 0,  # 底部边缘距离
                             l = 0,  # 左边边缘距离
                             unit = "cm"))
p.bsp.V1A
#export pdf 

# #reassortment SH test
# clades = c("3a1_sub","3a2","V1A")
# flu_segments = c("PB2","PB1","PA","HA","NP","NA","MP","NS")
# 
# pml.list = list()
# clade <- "3a1_sub"
# whole.ml.tree <- read.tree(sprintf("../analysis/all_segments/China-flub-wholegenome.%s.afa.treefile",gsub("_",".",clade)))
# whole.seg.align <- read.alignment(sprintf("../analysis/all_segments/China-flub-wholegenome.%s.afa",gsub("_",".",clade)),"fasta")
# whole.fit <- pml(whole.ml.tree,as.phyDat(whole.seg.align))
# whole.fit <- optim.pml(whole.fit)
# pml.list <- append(pml.list,list(whole.fit))
# 
# for (flu_segment in flu_segments){
#   ml.tree <- read.tree(sprintf("../analysis/all_segments/vic_China_%s_%s.afa.treefile",clade,flu_segment))
#   seg.align <- read.alignment(sprintf("../analysis/all_segments/vic_China_%s_%s.afa",clade,flu_segment),"fasta")
#   fit <- pml(ml.tree,as.phyDat(seg.align))
#   fit <- optim.pml(fit)
#   pml.list <- append(pml.list,list(fit))
# }
# 
# SH.test(pml.list, B=1000)


# #map after 2020-12
# gisaid.vic.china.meta$province <- gsub("-.*|B/","",gisaid.vic.china.meta$Isolate_Name)
# gisaid.vic.china.meta$province <- gsub("BVR","Sichuan",gisaid.vic.china.meta$province)
# 
# province.gps <- read_tsv("../analysis/iqtree/province_gps.txt")
# 
# month.intervals <- c("2020-12","2021-01","2021-02","2021-03","2021-04","2021-05","2021-06","2021-07","2021-08","2021-09")
# 
# API_pre = "http://xzqh.mca.gov.cn/data/"
# China_map = st_read(dsn = paste0(API_pre, "quanguo.json"),
#                     stringsAsFactors=FALSE)
# st_crs(China_map) = 4326
# 
# (data <- fromJSON(file = "../analysis/Nextstrain/results/clades.json"))
# data <- as.data.frame(unlist(data$nodes))
# data <- tibble::rownames_to_column(data)
# names(data) <- c("Isolate_Id","clade")
# data$Taxa <- gsub(".clade_membership","",data$Isolate_Id)
# 
# for(i in 1:(length(month.intervals)-1)){
#   #print(month.intervals[i])
#   gisaid.vic.china.meta.filter <- gisaid.vic.china.meta %>%
#     filter(Collection_Date > month.intervals[i] & Collection_Date < month.intervals[i+1]) %>%
#     merge(.,data,by = "Isolate_Id", all.x = TRUE) %>%
#     subset(,c("Collection_Date","clade","province")) %>%
#     group_by(Collection_Year=floor_date(as.Date(Collection_Date),"month"),clade,province) %>%
#     summarise(count = n()) %>%
#     pivot_wider(names_from = clade,values_from = count) %>%
#     merge(.,province.gps,by = "province", all.x = TRUE)
# 
#   gisaid.vic.china.meta.filter[is.na(gisaid.vic.china.meta.filter)] <- 0
# 
#   if(!"V1A" %in% names(gisaid.vic.china.meta.filter)){
#     gisaid.vic.china.meta.filter$`V1A` <- 0
#   }
#   if(!"3a1" %in% names(gisaid.vic.china.meta.filter)){
#     gisaid.vic.china.meta.filter$`3a1` <- 0
#   }
#   if(!"3a2" %in% names(gisaid.vic.china.meta.filter)){
#     gisaid.vic.china.meta.filter$`3a2` <- 0
#   }
# 
#   gisaid.vic.china.meta.filter$all_seqs <- gisaid.vic.china.meta.filter$`V1A`+gisaid.vic.china.meta.filter$`3a1`+gisaid.vic.china.meta.filter$`3a2`
# 
#   print(month.intervals[i])
#   print(gisaid.vic.china.meta.filter[c("province","all_seqs")])
#   #world map
#   worldMap <- fortify(map_data("world"), region = "subregion")%>%
#     filter(region == "China")
#   #table(worldMap$subregion)
#   #table(worldMap$region)
#   #colnames(worldMap)
# 
#   p <- ggplot(China_map)+
#     geom_sf()+
#     geom_scatterpie(aes(x=longitude, y=latitude, group=province,r=all_seqs/5),data=gisaid.vic.china.meta.filter, cols=c("V1A","3a1","3a2"), size=0.1,color="black", alpha = 0.8)+ #add pie chart
#     scale_fill_manual(values=c("#997ba5","#7ea34f","#50b2b2"),breaks=c("V1A","3a1","3a2"))+
#     theme(panel.grid =element_blank(),panel.background = element_blank(),legend.position="none") +   ## 删去网格线
#     theme(axis.text = element_blank()) +   ## 删去刻度标签
#     theme(axis.ticks = element_blank()) +   ## 删去刻度线
#     theme(panel.border = element_blank())+   ## 删去外层边框
#     theme(plot.margin = margin(t = 0,  # 顶部边缘距离
#                                r = 0,  # 右边边缘距离
#                                b = 0,  # 底部边缘距离
#                                l = 0,  # 左边边缘距离
#                                unit = "cm"))+
#     labs(x = "", y = "", title = month.intervals[i])
#   assign(paste("p", i, sep = ""), p)
# }
# 
# 
# ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, ncol = 3, nrow = 3)
