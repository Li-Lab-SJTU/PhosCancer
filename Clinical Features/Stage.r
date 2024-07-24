Cancertype='CCRCC'# BRCA CCRCC CRC GBM OV PDAC UCEC     EOGC HNSCC LSCC LUAD PAAD HCC 

Annotation <- read.table("TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)

clinical <- read.table(paste0("/dssg/home/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,".tsv"), header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE,sep='\t')

NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
PhosNAcutoff=0.80
ProteinNAcutoff=0.50
NormPhosphoGroupsStage<-data.frame(Stage=clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),clinical$'Cases Submitter ID'),'Tumor Stage'],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)
unique(NormPhosphoGroupsStage[,1])

if(length(which(is.na(NormPhosphoGroupsStage[,1])))>0){
NormPhosphoGroupsStage<-NormPhosphoGroupsStage[-which(is.na(NormPhosphoGroupsStage[,1])),] }
unique(NormPhosphoGroupsStage[,1])

NormPhosphoGroupsfiltered<-NormPhosphoGroupsStage[,-c(which(apply(NormPhosphoGroupsStage[,-1],2,function(x){length(which(is.na(x)))})>nrow(NormPhosphoGroupsStage[,-1])*PhosNAcutoff)+1)]
dim(NormPhosphoGroupsfiltered)
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
sites<-rownames(NormPhosphoGroupsfiltered)
if(length(sites)>0){
for(i in 2:ncol(NormPhosphoGroupsfiltered)){
whichsite=colnames(NormPhosphoGroupsfiltered)[i]
boxplotData<-na.omit(NormPhosphoGroupsfiltered[,c(1,which(colnames(NormPhosphoGroupsfiltered)==whichsite))])
names(boxplotData)[2]<-'Site'
colorss<-c("#F0EC6F","#55B77F","#326070","#3F0748")
names(colorss)<-c("I","II","III","IIII")
cols<-colorss[unique(boxplotData[,1])]
jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/Stagediffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_',Cancertype,'.jpeg'),width = 600, height = 600)
my_plot<-ggbetweenstats(data = boxplotData, 
               x = Stage, 
               y = Site,
			   xlab = "",
			   ylab = 'Phosphorylation level', 
			   messages = FALSE,
			   results.subtitle=F,
			   type='nonparametric',
			   pairwise.comparisons = TRUE, 
			   pairwise.display = "all", 
               pairwise.annotation = "p.value", 
               p.adjust.method = "none", 
			   point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha= 0.8, size = 6, stroke = 0),
			   ggsignif.args = list(textsize = 4, tip_length = 0.015),
			   centrality.point.args = list(size = 5, color = "darkred"),
			   centrality.label.args = list(size = 6, nudge_x = 0.4, segment.linetype = 4,min.segment.length = 0),
)+ ggplot2::theme(plot.title = element_text(hjust = 0.5, size=20,face = "bold"),axis.title=element_text(vjust=2, size=24,face = "bold"),axis.text=element_text(vjust=1,size=20,face = "bold"),axis.title.y.right = element_blank(), axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank())
ymax <- ggplot_build(my_plot)$layout$panel_params[[1]]$y.range[2]
print(my_plot+stat_compare_means(label.y = ymax,label.x =0.8,size = 6 ) )
dev.off()
}
}

