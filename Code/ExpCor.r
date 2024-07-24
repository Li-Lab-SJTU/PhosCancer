Cancertype='BRCA'# BRCA CCRCC CRC GBM OV PDAC UCEC     EOGC HNSCC LSCC LUAD PAAD HCC 
clinical <- read.table(paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,".tsv"), header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE,sep='\t')
PhosNAcutoff=0.80
ProteinNAcutoff=0.50

NormPhosphoGroups <- read.table(paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormProtein <- read.table(paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Proteome/Proteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE) #标准化之后的蛋白表达谱，磷酸化的交集样本，之后过滤NA比例
intersamples<-intersect(grep('Primary Tumor$',names(NormPhosphoGroups),value=T),grep('Primary Tumor$',names(NormProtein),value=T))
NormProteinfiltered<-NormProtein[,match(intersamples,names(NormProtein))]
NormProteinfiltered<-NormProteinfiltered[-which(apply(NormProteinfiltered,1,function(x){length(which(is.na(x)))})>ncol(NormProteinfiltered)*ProteinNAcutoff),]
NormPhosphoGroupsfiltered<-NormPhosphoGroups[,match(intersamples,names(NormPhosphoGroups))]
NormPhosphoGroupsfiltered<-NormPhosphoGroupsfiltered[-which(apply(NormPhosphoGroupsfiltered,1,function(x){length(which(is.na(x)))})>ncol(NormPhosphoGroupsfiltered)*PhosNAcutoff),]

sites<-setdiff(rownames(NormPhosphoGroupsfiltered),gsub('\\+','\\|',gsub(paste0('=',Cancertype),'',gsub('.png','',list.files('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/PNGs/ExpCor/',pattern=paste0('=',Cancertype,'.png'))))))

library(ggstatsplot)
library(ggplot2)
if(length(sites)>0){
for(i in 1:length(sites)){
whichsite=sites[i]
if(any(rownames(NormPhosphoGroupsfiltered)==whichsite)&any(rownames(NormProteinfiltered)==gsub('\\|p.*','',whichsite))){

plotdata<- data.frame(Phosphorylation=as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,]),Expression=as.numeric(NormProteinfiltered[rownames(NormProteinfiltered)==gsub('\\|p.*','',whichsite),]))

png(filename=paste0('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/PNGs/ExpCor/',gsub('\\|','+',whichsite),'=',Cancertype,'.png'),width = 600, height = 600)
print(ggscatterstats(
  plotdata,
  x = Phosphorylation,
  y = Expression,
  xlab = "Phosphorylation level (Log2 intensity)", # label for x axis
  ylab = "Protein abundance (Log2 intensity)", # label for y axis
  line.size = 1.5,
  line.color = "#FFD700", 
  marginal = TRUE, 
  marginal.type = "density",
  marginal.size = 5, 
  margins = "both",#c("both", "x", "y"), 
  width.jitter = NULL,
  height.jitter = NULL, xfill = "#81B0D0", yfill = "#85C4B0",
  centrality.para = 'median', type = "pearson", #spearman
  results.subtitle = F,
  #title = 'Correlation between  protein abundance versus phosphoprotein', caption = NULL, nboot = 100, beta = 0.1, k = 3,
  axes.range.restrict = FALSE, ggtheme = ggplot2::theme_minimal(),
  messages = F)#+ ggplot2::theme(plot.title = element_text(hjust = 0.5, size=18,face = "bold"),axis.title=element_text(vjust=2, size=18,face = "bold"),axis.text=element_text(vjust=1,size=18,face = "bold"))
)
  dev.off()
}
}
}