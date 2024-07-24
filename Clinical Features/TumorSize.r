Cancertype<-'CCRCC'
PhosNAcutoff=0.80
Annotation <- read.table("/dssg/home/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
library(readxl)
clinical<-as.data.frame(read_xlsx(paste0("/dssg/home/Projects/PhosCancer/Clinical/FromPaper/",Cancertype,".xlsx"),sheet = 2))
NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroupsClinical<-data.frame(clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),gsub('X','',clinical$'Case_ID')),'Tumor_Size_cm',drop=F],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)

NormPhosphoGroupsfiltered<-NormPhosphoGroupsClinical[,-c(which(apply(NormPhosphoGroupsClinical[,-1],2,function(x){length(which(is.na(x)))})>nrow(NormPhosphoGroupsClinical[,-1])*PhosNAcutoff)+1)]
dim(NormPhosphoGroupsfiltered)

Agediffer<-c()
for(i in 2:ncol(NormPhosphoGroupsfiltered)){ 

SiteData<-na.omit(NormPhosphoGroupsfiltered[,c(1,i)])

TumorSampleSize<-nrow(SiteData)
SpearmanP<-cor.test(SiteData[,1],SiteData[,2],method='spearman')$p.value
SpearmanE<-cor.test(SiteData[,1],SiteData[,2],method='spearman')$estimate
Agediffer<-rbind(Agediffer,data.frame(TumorSampleSize,SpearmanP,SpearmanE))

library(ggdist)        
library(gghalves)     
library(ggplot2)
library("ggpubr")
whichsite<-colnames(SiteData)[2]
colnames(SiteData)<-c('BMI','Site')
jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/TumorSizeCor/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_',Cancertype,'.jpeg'),width = 500, height = 500)
print(
ggscatter(SiteData, x = "BMI", y = "Site",xlab="Tumor size (cm)",ylab='Phosphorylation level',
          size = 2,  shape=8,                     
          add = "reg.line",                  
          conf.int = TRUE,                  
          rug = TRUE ,                    
		  color='#EEAD0E',
          add.params = list(color = "black",  
                          fill = "#EEAD0E"),
          ggtheme = theme_bw()
)+
  stat_cor(size=8,fontface = "bold",method = "spearman",r.accuracy=0.01,p.accuracy=0.01) 
+
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5,
								  face = "bold"),
        axis.title = element_text(size = 19, 
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
  ))
dev.off()
}
rownames(Agediffer)<-colnames(NormPhosphoGroupsfiltered)[2:ncol(NormPhosphoGroupsfiltered)]
Agediffer$SpearmanFDR<-p.adjust(Agediffer$SpearmanP,'BH')
write.table(Agediffer,file=paste0("/dssg/home/Projects/PhosCancer/Analysis/Quantitative/TumorSizeCor/",Cancertype,"_TumorSizeCor.tsv"),sep="\t",quote=F,row.names=T,col.names=T)
