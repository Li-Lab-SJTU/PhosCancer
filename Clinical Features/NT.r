Cancertype='PDAC'

Annotation <- read.table("/dssg/home/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)

clinical <- read.table(paste0("/dssg/home/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,".tsv"), header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE,sep='\t')

NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
PhosNAcutoff=0.80
ProteinNAcutoff=0.50
NormPhosphoGroupsfiltered<-NormPhosphoGroups[-which(apply(NormPhosphoGroups,1,function(x){length(which(is.na(x[grep('Primary Tumor$',colnames(NormPhosphoGroups))])))})>length(grep('Primary Tumor$',colnames(NormPhosphoGroups)))*PhosNAcutoff|apply(NormPhosphoGroups,1,function(x){length(which(is.na(x[grep('Solid Tissue Normal$',colnames(NormPhosphoGroups))])))})>length(grep('Solid Tissue Normal$',colnames(NormPhosphoGroups)))*PhosNAcutoff),]
dim(NormPhosphoGroupsfiltered)

library(ggdist)        ## halfeye plots
library(gghalves)      ## off-set jitter
library(ggplot2)
library("ggpubr")

sites<-rownames(NormPhosphoGroupsfiltered)
if(length(sites)>0){
for(i in 1:length(sites)){
whichsite=sites[i]
boxplotData<-na.omit(data.frame(Exp=as.numeric(c(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,grep('Primary Tumor$',colnames(NormPhosphoGroupsfiltered))],NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,grep('Solid Tissue Normal$',colnames(NormPhosphoGroupsfiltered))])),Type=factor(c(rep('Tumor',length(grep('Primary Tumor$',colnames(NormPhosphoGroupsfiltered)))),rep('Normal',length(grep('Solid Tissue Normal$',colnames(NormPhosphoGroupsfiltered))))),levels=c('Tumor','Normal'))))

gsub('/','',Annotation[Annotation[,1]==whichsite,10])
jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/NTdiffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_',Cancertype,'.jpeg'),width = 700, height = 500)
print(ggplot(boxplotData, aes(x = Type, y = Exp, color = Type, fill = Type)) +
  scale_color_manual(limits=c("Normal","Tumor"), values=c("#5F80B4","#922927"))+
  scale_fill_manual(limits=c("Normal","Tumor"), values=c("#5F80B4","#922927")) + 
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
  )+scale_x_discrete(labels = c("Tumor" =paste0("Tumor",' (',table(boxplotData[,2])[1],') '),"Normal" = paste0("Normal",' (',table(boxplotData[,2])[2],') ')))+
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
}