PhosNAcutoff=0.80

readtable<-function(Cancertype){
Annotation <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
# clinical <- read.table(paste0("E:/PhD/LiLab/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,".tsv"), header = TRUE,fill=T,sep='\t',check.names=F,stringsAsFactors = FALSE)
library(readxl)
clinical <-read.table(paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,".tsv"), header = TRUE,fill=T,sep='\t',check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroups <- read.table(paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroupsClinical<-data.frame(clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),gsub('X','',clinical$'Cases Submitter ID')),'Race',drop=F],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)

if(Cancertype=='LUAD'){
library(stringr)
NormPhosphoGroupsClinical[,1]<-str_to_title(NormPhosphoGroupsClinical[,1])
NormPhosphoGroupsClinical[,1]<-gsub(' Or ',' or ',NormPhosphoGroupsClinical[,1])
}
####
#缺失值过滤
NormPhosphoGroupsClinical[,1]
if(any(is.na(NormPhosphoGroupsClinical[,1]))){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(is.na(NormPhosphoGroupsClinical[,1])),]}
if(any(NormPhosphoGroupsClinical[,1]=='Not Reported')){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(NormPhosphoGroupsClinical[,1]=='Not Reported'),]}
if(any(NormPhosphoGroupsClinical[,1]=='Unknown')){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(NormPhosphoGroupsClinical[,1]=='Unknown'),]}
if(any(NormPhosphoGroupsClinical[,1]=='Other')){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(NormPhosphoGroupsClinical[,1]=='Other'),]}
NormPhosphoGroupsClinical[,1]

NormPhosphoGroupsfiltered<-NormPhosphoGroupsClinical[,-c(which(apply(NormPhosphoGroupsClinical[,-1],2,function(x){length(which(is.na(x)))})>nrow(NormPhosphoGroupsClinical[,-1])*PhosNAcutoff)+1)]
if(nrow(NormPhosphoGroupsfiltered)!=0){
if(length(unique(NormPhosphoGroupsfiltered[,1]))>1){
return(NormPhosphoGroupsfiltered)
}
}
}
BRCA<-readtable('BRCA')
CCRCC<-readtable('CCRCC')
CRC<-readtable('CRC')
EOGC<-readtable('EOGC')
GBM<-readtable('GBM')
HCC<-readtable('HCC')
HCC_cnhpp<-readtable('HCC_cnhpp')
HNSCC<-readtable('HNSCC')
LSCC<-readtable('LSCC')
LUAD<-readtable('LUAD')
LUAD_cnhpp<-readtable('LUAD_cnhpp')
OV<-readtable('OV')
PDAC<-readtable('PDAC')
UCEC<-readtable('UCEC')

library(ggstatsplot)
library(ggplot2)

files<-c('BRCA','CCRCC','CRC','EOGC','GBM','HCC','HCC_cnhpp', 'HNSCC','LSCC','LUAD','LUAD_cnhpp','OV','PDAC','UCEC' )

Annotation <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
sites<-Annotation[,1]

library(reshape)
if(length(sites)>0){
for(site in 1:length(sites)){
whichsite=sites[site]
Total<-c()
for(i in 1:length(files)){
singles <-get(files[i])
if(length(which(colnames(singles)==whichsite))!=0){
singles<-na.omit(singles[,c(1,which(colnames(singles)==whichsite))])
names(singles)[2]<-'Site'
singles<-cbind(singles,Cancertype=files[i])
Total<-rbind(Total,singles)
	}
}
colorss<-c("#7570B3","#D95F02","#1B9E77","#E7298A")
names(colorss)<-c("Asian","White","Black or African American","American Indian or Alaska Native")

library(ggpubr)
if(!is.null(Total)){
jpeg(filename=paste0('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/PNGs/Racediffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_ALL.jpeg'),width = 1000, height = 600)
print(ggplot(Total,aes(Cancertype,Site,color = Race)) + stat_boxplot(mapping=aes(x=Cancertype,y=Site),geom='errorbar',size=0.8)+# stat_boxplot() 函数用来给箱子添加median ±1.5 IQR 的最大值和最小值的bar, 注意分组用的参数的设置position,根据图像大小设置bar间距width参数
  geom_boxplot(outlier.shape = 21,outlier.fill  = 'black',outlier.size = 0.25,size=1.6) + 
  labs(x="",y="Phosphorylation level")+# 添加标题，x轴，y轴内容title="Distribution of phosphorylation level across different cancertype",
  theme_grey() +
  theme(legend.position="top",
        strip.text.x = element_text(size = 12, colour = "black",face='bold'),
        panel.grid.major.x = element_line(size = 1),
        panel.grid.major.y = element_line(size = 1),
        panel.grid.minor.x = element_line(size = 0.8),
        panel.grid.minor.y = element_line(size = 0.8),
        panel.background = element_blank(),  # 去掉背景色
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


