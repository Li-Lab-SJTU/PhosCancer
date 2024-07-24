Cancertype<-'BRCA'
PhosNAcutoff=0.80
library(readxl)
clinical<-as.data.frame(read_xlsx(paste0("/dssg/home/Projects/PhosCancer/Clinical/FromPaper/",Cancertype,".xlsx"),sheet = 2))
NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroupsClinical<-data.frame(clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),gsub('X','',clinical$'Sample.IDs')),'NMF.Cluster',drop=F],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)
if(any(is.na(NormPhosphoGroupsClinical[,1]))){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(is.na(NormPhosphoGroupsClinical[,1])),]}
NormPhosphoGroupsfiltered<-NormPhosphoGroupsClinical[,-c(which(apply(NormPhosphoGroupsClinical[,-1],2,function(x){length(which(is.na(x)))})>nrow(NormPhosphoGroupsClinical[,-1])*PhosNAcutoff)+1)]
BRCA<-NormPhosphoGroupsfiltered

Cancertype<-'HNSCC'
library(readxl)
clinical<-as.data.frame(read_xlsx(paste0("/dssg/home/Projects/PhosCancer/Clinical/FromPaper/",Cancertype,".xlsx"),sheet = 2))
NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroupsClinical<-data.frame(clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),gsub('X','',clinical$'case_id')),'integrated_subtype',drop=F],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)
if(any(NormPhosphoGroupsClinical[,1]=='NA')){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(NormPhosphoGroupsClinical[,1]=='NA'),]}
NormPhosphoGroupsfiltered<-NormPhosphoGroupsClinical[,-c(which(apply(NormPhosphoGroupsClinical[,-1],2,function(x){length(which(is.na(x)))})>nrow(NormPhosphoGroupsClinical[,-1])*PhosNAcutoff)+1)]
HNSCC<-NormPhosphoGroupsfiltered

Cancertype<-'LSCC'
library(readxl)
clinical<-as.data.frame(read_xlsx(paste0("/dssg/home/Projects/PhosCancer/Clinical/FromPaper/",Cancertype,".xlsx"),sheet = 2))
NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroupsClinical<-data.frame(clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),gsub('\\.','-',clinical$'Sample.ID')),'NMF.subtype',drop=F],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)
if(any(is.na(NormPhosphoGroupsClinical[,1]))){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(is.na(NormPhosphoGroupsClinical[,1])),]}
NormPhosphoGroupsfiltered<-NormPhosphoGroupsClinical[,-c(which(apply(NormPhosphoGroupsClinical[,-1],2,function(x){length(which(is.na(x)))})>nrow(NormPhosphoGroupsClinical[,-1])*PhosNAcutoff)+1)]
LSCC<-NormPhosphoGroupsfiltered

Cancertype<-'LUAD'
library(readxl)
clinical<-as.data.frame(read_xlsx(paste0("/dssg/home/Projects/PhosCancer/Clinical/FromPaper/",Cancertype,".xlsx"),sheet = 2))
NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroupsClinical<-data.frame(clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),gsub('\\.','-',clinical$'Sample.ID')),'NMF.consensus',drop=F],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)
if(any(is.na(NormPhosphoGroupsClinical[,1]))){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(is.na(NormPhosphoGroupsClinical[,1])),]}
NormPhosphoGroupsfiltered<-NormPhosphoGroupsClinical[,-c(which(apply(NormPhosphoGroupsClinical[,-1],2,function(x){length(which(is.na(x)))})>nrow(NormPhosphoGroupsClinical[,-1])*PhosNAcutoff)+1)]
LUAD<-NormPhosphoGroupsfiltered


library(ggstatsplot)
library(ggplot2)

files<-c('BRCA','HNSCC','LSCC','LUAD')

Annotation <- read.table("TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
sites<-Annotation[,1]

library(reshape)
library(RColorBrewer)
if(length(sites)>0){
for(site in 1:length(sites)){
whichsite=sites[site]
Total<-c()
for(i in 1:length(files)){
singles <-get(files[i])
if(length(which(colnames(singles)==whichsite))!=0){
singles<-na.omit(singles[,c(1,which(colnames(singles)==whichsite))])
names(singles)<-c('Subtype','Site')
singles<-cbind(singles,Cancertype=files[i])
Total<-rbind(Total,singles)
	}
}
Total[,1]<-factor(Total[,1],levels=c("LumA-I","LumB-I","Basal-I","HER2-I","CIN","Immune","Basal","Proliferative primitive","EMT enriched","Classical","Basal inclusive","Inflamed secretory","C1","C2","C3","C4"))
colorss<-c(brewer.pal(9, "Blues")[c(3,5,7,9)],brewer.pal(9, "Greens")[c(3,5,7)],brewer.pal(9, "YlOrRd")[c(1,3,5,7,9)],brewer.pal(9, "Purples")[c(3,5,7,9)])
names(colorss)<-c("LumA-I","LumB-I","Basal-I","HER2-I","CIN","Immune","Basal","Proliferative primitive","EMT enriched","Classical","Basal inclusive","Inflamed secretory","C1","C2","C3","C4")

library(ggpubr)
if(!is.null(Total)){
jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/MajorSubtypediffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_ALL.jpeg'),width = 1000, height = 600)
print(ggplot(Total,aes(Cancertype,Site,color = Subtype)) + stat_boxplot(mapping=aes(x=Cancertype,y=Site),geom='errorbar',size=0.8)+
  geom_boxplot(outlier.shape = 21,outlier.fill  = 'black',outlier.size = 0.25,size=1.6) + 
  labs(x="",y="Phosphorylation level")+
  theme_grey() +
  theme(legend.position="top",
        strip.text.x = element_text(size = 12, colour = "black",face='bold'),
        panel.grid.major.x = element_line(size = 1),
        panel.grid.major.y = element_line(size = 1),
        panel.grid.minor.x = element_line(size = 0.8),
        panel.grid.minor.y = element_line(size = 0.8),
        panel.background = element_blank(),  
        axis.line.x = element_line(colour = 'black', size = 0.8),
        axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = 'black', size = 1.5),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour = 'black', size = 18,face='bold'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(size = 0)
		) + 
  theme(text = element_text(size=16)) + 
  scale_color_manual(values=colorss)+ 
  stat_compare_means(size=6,label = "p.signif")+
  facet_grid(.~Cancertype,scales="free") )
dev.off()
	}
}
}


