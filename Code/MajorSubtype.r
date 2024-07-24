Cancertype<-'BRCA'
PhosNAcutoff=0.80
Annotation <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
# clinical <- read.table(paste0("E:/PhD/LiLab/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,".tsv"), header = TRUE,fill=T,sep='\t',check.names=F,stringsAsFactors = FALSE)
library(readxl)
clinical<-as.data.frame(read_xlsx(paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Clinical/FromPaper/",Cancertype,".xlsx"),sheet = 2))
NormPhosphoGroups <- read.table(paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroupsClinical<-data.frame(clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),gsub('X','',clinical$'Sample.IDs')),'NMF.Cluster',drop=F],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)


####
#缺失值过滤
NormPhosphoGroupsClinical[,1]
if(any(is.na(NormPhosphoGroupsClinical[,1]))){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(is.na(NormPhosphoGroupsClinical[,1])),]}
NormPhosphoGroupsClinical[,1]

NormPhosphoGroupsfiltered<-NormPhosphoGroupsClinical[,-c(which(apply(NormPhosphoGroupsClinical[,-1],2,function(x){length(which(is.na(x)))})>nrow(NormPhosphoGroupsClinical[,-1])*PhosNAcutoff)+1)]

if(length(unique(NormPhosphoGroupsfiltered[,1]))>1){

Agediffer<-c()
for(i in 2:ncol(NormPhosphoGroupsfiltered)){ 

SiteData<-na.omit(NormPhosphoGroupsfiltered[,c(1,i)])
whichsite<-colnames(SiteData)[2]
names(SiteData)<-c('Race','Site')
if(length(unique(SiteData[,1]))>1){
#计算
TumorSampleSize<-nrow(SiteData)
TumorSampleMean<-mean(2^SiteData[,2])

if(length(unique(SiteData[,1]))>2){
kruskalP<-kruskal.test(Site~Race,data=SiteData)$p.value
	}
if(length(unique(SiteData[,1]))==2){
kruskalP<-wilcox.test(Site~Race,data=SiteData)$p.value
	}
temp<-data.frame(TumorSampleSize,TumorSampleMean,kruskalP)
rownames(temp)<-whichsite	
Agediffer<-rbind(Agediffer,temp)
		
#画图
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
boxplotData<-SiteData
#${UniProtID}_${Position}_${cancerType}.jpeg

jpeg(filename=paste0('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/PNGs/MajorSubtypediffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_',Cancertype,'.jpeg'),width = 600, height = 600)
my_plot<-ggbetweenstats(data = boxplotData, 
               x = Race, 
               y = Site,
			   xlab = "",
			   ylab = 'Phosphorylation level', 
			   messages = FALSE,
			   results.subtitle=F,
			   #point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6),alpha = 0.5, size = 1.5, stroke = 0),
			   type='nonparametric',
			   pairwise.comparisons = TRUE, # display results from pairwise comparisons
			   pairwise.display = "all", #significant display only significant pairwise comparisons
               pairwise.annotation = "p.value", # annotate the pairwise comparisons using p-values
               p.adjust.method = "none", # adjust p-values for multiple tests using this method
			   point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha= 0.8, size = 6, stroke = 0),
			   #title = "Phosphorylation level (Log2 intensity) by tumor stage", #TMT ratio 还是Relative abundance 
			   ggsignif.args = list(textsize = 4, tip_length = 0.015),
			   #boxplot.args = list(width = 0),
			   #violin.args = list(width = 0.6, alpha = 0.2),
			   centrality.point.args = list(size = 5, color = "darkred"),
			   centrality.label.args = list(size = 6, nudge_x = 0.4, segment.linetype = 4,min.segment.length = 0),
 			   #ggtheme = ggplot2::theme_bw(base_size = 15),
)+ ggplot2::theme(plot.title = element_text(hjust = 0.5, size=20,face = "bold"),axis.title=element_text(vjust=2, size=24,face = "bold"),axis.text=element_text(vjust=1,size=20,face = "bold"),axis.title.y.right = element_blank(), axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank())# 改变点颜色
ymax <- ggplot_build(my_plot)$layout$panel_params[[1]]$y.range[2]
print(my_plot+stat_compare_means(label.y = ymax,label.x =0.8,size = 6 ) )
dev.off()

}
}
Agediffer$kruskalFDR<-p.adjust(Agediffer$kruskalP,'BH')
write.table(Agediffer,file=paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/Quantitative/MajorSubtypediffer/",Cancertype,"_MajorSubtypediffer.tsv"),sep="\t",quote=F,row.names=T,col.names=T)
}
