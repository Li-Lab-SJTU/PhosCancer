Cancertype<-'GBM'

Annotation <- read.table("/dssg/home/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)

PhosNAcutoff=0.80
ProteinNAcutoff=0.50
NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
gsva_matrix <- read.table(paste0("/dssg/home/Projects/PhosCancer/Proteome/Proteome_gsva/",Cancertype,".tsv"),comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)

intersamples<-intersect(grep('Primary Tumor$',names(NormPhosphoGroups),value=T),grep('Primary Tumor$',names(gsva_matrix),value=T))

NormPhosphoGroupsfiltered<-NormPhosphoGroups[,match(intersamples,names(NormPhosphoGroups))]
NormPhosphoGroupsfiltered<-NormPhosphoGroupsfiltered[-which(apply(NormPhosphoGroupsfiltered,1,function(x){length(which(is.na(x)))})>ncol(NormPhosphoGroupsfiltered)*PhosNAcutoff),]

gsva_matrix<-gsva_matrix[,match(intersamples,names(gsva_matrix))]

sites<-rownames(NormPhosphoGroupsfiltered)

#Table 
HallmarkCorTotal<-c()
HallmarkCorTop<-c()
for(i in 1:length(sites)){
whichsite=sites[i]
corP<-apply(gsva_matrix,1,function(x){cor.test(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,]),x,method='spearman')$p.value})
corE<-apply(gsva_matrix,1,function(x){cor.test(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,]),x,method='spearman')$estimate})

Pcolorss<-as.character(ifelse(corP<=0.05,'brown2','black'))#
Pfontfaces<-as.numeric(ifelse(corP<=0.05,2,1))
Psizes<-as.numeric(ifelse(corP<=0.05,4,3))
bardata<-data.frame(Hallmark=gsub('HALLMARK_','',names(corE)),Cor=abs(as.numeric(corE)),id=seq(1,length(corE))) 
library(ggplot2)
num<-length(corP)/2
bardata$angle<-ifelse(bardata$id<=num,96-bardata$id*(180/num),96-bardata$id*(180/num)+175)
bardata$hjust<-ifelse(bardata$id<=num,0,1)
bardata$Correlation<-as.character(ifelse(corE<=0,'Neg','Pos'))
bardata$Pcolors<-as.character(ifelse(corP<=0.05,'p<0.05','p≥0.05'))

limitsstart<-ifelse(max(abs(corE))>0.6,-0.25,-0.15)

jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/HallmarkCor/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_',Cancertype,'.jpeg'),width = 800, height = 800)
print(ggplot(bardata,aes(x=as.factor(id),y=Cor,fill=Correlation))+
  geom_bar(stat="identity")+
  coord_polar()+  scale_y_continuous(limits = c(limitsstart, max(abs(as.numeric(corE)))+0.1), expand = c(-0.1, 0)) +
  geom_text(aes(x=as.factor(id),y=abs(as.numeric(corE))+0.01,label=Hallmark, angle=angle,hjust=hjust,colour = Pcolors),
            size=Psizes,fontface = Pfontfaces)+
  theme_minimal()+ylab("The absolute value of correlation")+xlab("")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size= 15),
        axis.title.y = element_text(size= 20, face = "bold"),
		legend.position= c(0.9,0.12),
		legend.text = element_text(size=12),
        legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(0.6, 'cm'),
		legend.margin = margin(t = 0, r = 10, b = 0, l = 10),#上、右、下、左
		legend.title=element_blank()
		) +theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "in"))+
  scale_color_manual(values = c('brown2','black'))+		
  scale_fill_manual(values = c('#9FB3D2','#BE7F7D')))
dev.off()

SpearmanP=corP
SpearmanE=corE
PearsonP=apply(gsva_matrix,1,function(x){cor.test(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,]),x,method='pearson')$p.value})
PearsonE=apply(gsva_matrix,1,function(x){cor.test(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,]),x,method='pearson')$estimate})

whichhallmark<-names(SpearmanE[which(SpearmanP<0.05)][order(abs(SpearmanE[which(SpearmanP<0.05)]),decreasing=T)][1])
tempTop<-data.frame(Sites=rownames(NormPhosphoGroupsfiltered)[i],TumorSampleSize=ncol(NormPhosphoGroupsfiltered)-length(which(is.na(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,])))),
MostHallmark=gsub('HALLMARK_','',whichhallmark),
SpearmanP=SpearmanP[names(SpearmanP)==whichhallmark],
SpearmanFDR=p.adjust(SpearmanP,'BH')[names(SpearmanP)==whichhallmark],
SpearmanE=SpearmanE[names(SpearmanE)==whichhallmark],
PearsonP=PearsonP[names(PearsonP)==whichhallmark],
PearsonFDR=p.adjust(PearsonP,'BH')[names(SpearmanP)==whichhallmark],
PearsonE=PearsonE[names(PearsonE)==whichhallmark])
HallmarkCorTop<-rbind(tempTop,HallmarkCorTop)
#}
tempTotal<-data.frame(Sites=rownames(NormPhosphoGroupsfiltered)[i],TumorSampleSize=ncol(NormPhosphoGroupsfiltered)-length(which(is.na(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,])))),
SpearmanP=SpearmanP,
SpearmanFDR=p.adjust(SpearmanP,'BH'),
SpearmanE=SpearmanE,
PearsonP=PearsonP,
PearsonFDR=p.adjust(PearsonP,'BH'),
PearsonE=PearsonE,
Pathway=names(SpearmanP))
HallmarkCorTotal<-rbind(tempTotal,HallmarkCorTotal)
}
write.table(HallmarkCorTop,file=paste0("/dssg/home/Projects/PhosCancer/Analysis/Quantitative/",Cancertype,"_HallmarkCorTop.tsv"),sep="\t",quote=F,row.names=F,col.names=T)
write.table(HallmarkCorTotal,file=paste0("/dssg/home/Projects/PhosCancer/Analysis/Quantitative/",Cancertype,"_HallmarkCorTotal.tsv"),sep="\t",quote=F,row.names=F,col.names=T)
