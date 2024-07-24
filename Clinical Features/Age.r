Cancertype<-'BRCA'
PhosNAcutoff=0.80
Annotation <- read.table("TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
library(readxl)
clinical <-read.table(paste0("/dssg/home/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,".tsv"), header = TRUE,fill=T,sep='\t',check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroupsClinical<-data.frame(clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),gsub('X','',clinical$'Cases Submitter ID')),'Age at Diagnosis',drop=F],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)

NormPhosphoGroupsClinical$'Age at Diagnosis'<-as.numeric(NormPhosphoGroupsClinical$'Age at Diagnosis')

NormPhosphoGroupsClinical[,1]
if(any(is.na(NormPhosphoGroupsClinical[,1]))){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(is.na(NormPhosphoGroupsClinical[,1])),]}
if(any(NormPhosphoGroupsClinical[,1]==0)){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(NormPhosphoGroupsClinical[,1]==0),]}
NormPhosphoGroupsClinical[,1]
NormPhosphoGroupsfiltered<-NormPhosphoGroupsClinical[,-c(which(apply(NormPhosphoGroupsClinical[,-1],2,function(x){length(which(is.na(x)))})>nrow(NormPhosphoGroupsClinical[,-1])*PhosNAcutoff)+1)]
dim(NormPhosphoGroupsfiltered)

if(nrow(NormPhosphoGroupsfiltered)!=0){
Agediffer<-c()
for(i in 2:ncol(NormPhosphoGroupsfiltered)){ 

SiteData<-na.omit(NormPhosphoGroupsfiltered[,c(1,i)])
Control<-which(SiteData[,1]<median(SiteData[,1]))
Case<-which(SiteData[,1]>=median(SiteData[,1]))

TumorSampleSize<-nrow(SiteData)
AgeCutoff<-median(SiteData[,1])
YoungerSampleMean<-mean(2^SiteData[Control,2])
OlderSampleMean<-mean(2^SiteData[Case,2])
log2FC<-log2(mean(2^SiteData[Case,2])/mean(2^SiteData[Control,2]))
wilcoxP<-wilcox.test(SiteData[Case,2],SiteData[Control,2])$p.value
Agediffer<-rbind(Agediffer,data.frame(TumorSampleSize,AgeCutoff,OlderSampleMean,YoungerSampleMean,log2FC,wilcoxP))

library(ggdist)      
library(gghalves)    
library(ggplot2)
library("ggpubr")
whichsite<-colnames(SiteData)[2]
boxplotData<-SiteData
boxplotData[Control,1]<-'Control'
boxplotData[Case,1]<-'Case'
boxplotData<-data.frame(Exp=boxplotData[,2],Type=boxplotData[,1])
colnames(SiteData)[2]
jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/Agediffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_',Cancertype,'.jpeg'),width = 700, height = 500)
print(ggplot(boxplotData, aes(x = Type, y = Exp, color = Type, fill = Type)) +
  scale_color_manual(limits=c("Control","Case"), values=c("#5F80B4","#922927"))+
  scale_fill_manual(limits=c("Control","Case"), values=c("#5F80B4","#922927")) + 
  ggdist::stat_halfeye(  trim = F,
    adjust = 0.5,
    width = .5, alpha=0.5,
    position = position_nudge(x = .05)
  ) +
  gghalves::geom_half_point(
    side = "1", 
    range_scale = .4, 
    alpha = .5, size = 3
  ) +
    stat_boxplot(geom="errorbar",width=0.1,size=1.5,position=position_nudge(x=0.05),color="black")+
    geom_boxplot(
    width = .15, color = "black",position=position_nudge(x=0.05),
    size = 1.5, outlier.shape = NA
  ) +
  stat_compare_means(aes(group=Type),method='wilcox.test',vjust=-20.5,hjust=1,size=8)+
  theme_classic(  
    base_line_size = 1 
  )+scale_x_discrete(labels = c("Case" =paste0("Older",' (',table(boxplotData[,2])[names(table(boxplotData[,2]))=='Case'],') '),"Control" = paste0("Younger",' (',table(boxplotData[,2])[names(table(boxplotData[,2]))=='Control'],') ')))+
  labs(x="",y="Phosphorylation level")+
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5,
								  face = "bold"),
        axis.title.x = element_text(size = 19, 
                                    color = "black",
                                    face = "bold", 
									),
        legend.title = element_text(color="black", 
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", 
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 18,  
                                   color = "black", 
                                   face = "bold",
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0), 
        axis.text.y = element_text(size = 19,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )+theme(legend.position = 'none')+coord_flip()
)
dev.off()
}
rownames(Agediffer)<-colnames(NormPhosphoGroupsfiltered)[2:ncol(NormPhosphoGroupsfiltered)]
Agediffer$wilcoxFDR<-p.adjust(Agediffer$wilcoxP,'BH')
write.table(Agediffer,file=paste0("/dssg/home/Projects/PhosCancer/Analysis/Quantitative/Agediffer/",Cancertype,"_Agediffer.tsv"),sep="\t",quote=F,row.names=T,col.names=T)
}