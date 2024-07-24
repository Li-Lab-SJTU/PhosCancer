PhosNAcutoff=0.80

readtable<-function(Cancertype){
library(readxl)
clinical <-read.table(paste0("/dssg/home/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,".tsv"), header = TRUE,fill=T,sep='\t',check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroupsClinical<-data.frame(clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),gsub('X','',clinical$'Cases Submitter ID')),'Gender',drop=F],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)

if(Cancertype=="HCC_cnhpp"){
NormPhosphoGroupsClinical[NormPhosphoGroupsClinical[,1]=='M',1]<-'Male'
NormPhosphoGroupsClinical[NormPhosphoGroupsClinical[,1]=='F',1]<-'Female'
}
if(Cancertype=="LUAD_cnhpp"){
NormPhosphoGroupsClinical[NormPhosphoGroupsClinical[,1]==1,1]<-'Male'
NormPhosphoGroupsClinical[NormPhosphoGroupsClinical[,1]==2,1]<-'Female'
}
if(Cancertype=="LUAD"){
NormPhosphoGroupsClinical[NormPhosphoGroupsClinical[,1]=='male',1]<-'Male'
NormPhosphoGroupsClinical[NormPhosphoGroupsClinical[,1]=='female',1]<-'Female'
}

NormPhosphoGroupsClinical[,1]
if(any(NormPhosphoGroupsClinical[,1]=='Not Reported')){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(NormPhosphoGroupsClinical[,1]=='Not Reported'),]}
NormPhosphoGroupsClinical[,1]

if(length(unique(NormPhosphoGroupsClinical[,1]))>1){
NormPhosphoGroupsfiltered<-NormPhosphoGroupsClinical[,-which(apply(NormPhosphoGroupsClinical[NormPhosphoGroupsClinical[,1]=='Male',],2,function(x){length(which(is.na(x)))})>length(which(NormPhosphoGroupsClinical[,1]=='Male'))*PhosNAcutoff|apply(NormPhosphoGroupsClinical[NormPhosphoGroupsClinical[,1]=='Female',],2,function(x){length(which(is.na(x)))})>length(which(NormPhosphoGroupsClinical[,1]=='Female'))*PhosNAcutoff)]
return(NormPhosphoGroupsfiltered)
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

files<-c('BRCA','CCRCC','CRC','EOGC','GBM','HCC','HCC_cnhpp', 'HNSCC','LSCC','LUAD','LUAD_cnhpp','OV','PDAC','UCEC' )

Annotation <- read.table("/dssg/home/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
sites<-Annotation[,1]
library(reshape)
if(length(sites)>0){
for(site in 1:length(sites)){
whichsite=sites[site]
Total<-c()
for(i in 1:length(files)){
singles <-get(files[i])
if(length(which(colnames(singles)==whichsite))!=0){
SiteData<-na.omit(singles[,c(1,which(colnames(singles)==whichsite))])
SiteData<-data.frame(value=SiteData[,2],Sampletype=SiteData[,1],Cancertype=files[i])
Total<-rbind(Total,SiteData)
	}
}

library(ggpubr)
if(!is.null(Total)){
jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/Genderdiffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_ALL.jpeg'),width = 1000, height = 600)
print(ggplot(Total,aes(Cancertype,value,fill = Sampletype)) + stat_boxplot(mapping=aes(x=Cancertype,y=value),geom='errorbar',size=1.0)+
  geom_boxplot(outlier.shape = 21,outlier.fill  = 'black',outlier.size = 0.25,color = "black",size=1.0) + 
  labs(x="",y="Phosphorylation level")+
  theme_grey() +
  theme(legend.position="top",
        strip.text.x = element_text(size = 12, colour = "black",face='bold'),
        panel.grid.major.x = element_line(size = 1),
        panel.grid.major.y = element_line(size = 1),
        panel.grid.minor.x = element_line(size = 0.8),
        panel.grid.minor.y = element_line(size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = 'black', size = 1.5),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour = 'black', size = 18,face='bold'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(size = 0)
		) + 
  theme(text = element_text(size=16)) + 
  scale_fill_manual(limits=c("Female","Male"),values=c("#5F80B4","#922927"))+ 
  stat_compare_means(aes(group = Sampletype,label = ..p.signif..),method='wilcox.test',vjust=-0.5,size=6)+
  facet_grid(.~Cancertype,scales="free") )
dev.off()
}
 }
}