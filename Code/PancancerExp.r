files<-c('BRCA','CCRCC','CRC','EOGC','GBM','HCC','HCC_cnhpp', 'HNSCC','LSCC','LUAD','LUAD_cnhpp','OV','PDAC','UCEC' )
BRCA <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/BRCA.tsv',sep='\t',fill=T,header=T,check.names=F)
CCRCC <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/CCRCC.tsv',sep='\t',fill=T,header=T,check.names=F)
CRC <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/CRC.tsv',sep='\t',fill=T,header=T,check.names=F)
EOGC <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/EOGC.tsv',sep='\t',fill=T,header=T,check.names=F)
GBM <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/GBM.tsv',sep='\t',fill=T,header=T,check.names=F)
HCC <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/HCC.tsv',sep='\t',fill=T,header=T,check.names=F)
HCC_cnhpp <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/HCC_cnhpp.tsv',sep='\t',fill=T,header=T,check.names=F)
HNSCC <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/HNSCC.tsv',sep='\t',fill=T,header=T,check.names=F)
LSCC <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/LSCC.tsv',sep='\t',fill=T,header=T,check.names=F)
LUAD <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/LUAD.tsv',sep='\t',fill=T,header=T,check.names=F)
LUAD_cnhpp <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/LUAD_cnhpp.tsv',sep='\t',fill=T,header=T,check.names=F)
OV <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/OV.tsv',sep='\t',fill=T,header=T,check.names=F)
PDAC <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/PDAC.tsv',sep='\t',fill=T,header=T,check.names=F)
UCEC <-read.table('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/UCEC.tsv',sep='\t',fill=T,header=T,check.names=F)

library(reshape)
library("ggplot2")
library("ggpubr")
library(RColorBrewer)
load('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/PNGs/colorss.RData')

Annotation <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
sites<-Annotation[,1]

#whichsite=Annotation[Annotation[,4]=='P27361'&Annotation[,3]=='202',1]
if(length(sites)>0){
for(site in 1:length(sites)){
whichsite=sites[site]
Total<-c()
medians<-c()
for(i in 1:length(files)){
singles <-get(files[i])
if(any(rownames(singles)==whichsite)){
singles<-data.frame(melt(singles[which(rownames(singles)==whichsite),grep('Primary Tumor$',colnames(singles))]),Cancertype=files[i])
if(length(na.omit(singles[,2])!=0)){
singles[,3]<-paste0(files[i],' (',length(na.omit(singles[,2])),')')
Total<-rbind(Total,singles)
temp<-median(na.omit(singles[,2]))
names(temp)<-paste0(files[i],' (',length(na.omit(singles[,2])),')')
medians<-c(medians,temp)}
	}
}
if(!is.null(Total)){
Total[,3]<-factor(Total[,3],levels=names(medians)[order(medians,decreasing=F)])
selectedcolorss<-colorss[match(gsub(' \\(.*','',levels(Total[,3])),names(colorss))]
names(selectedcolorss)<-levels(Total[,3])

xtextangle = ifelse(length(levels(Total[,3]))>8,30,0)
jpeg(filename=paste0('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/PNGs/PancancerExp/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_Tumor.jpeg'),width = 1000, height = 600)
par(mar = c(5, 7, 4, 2) + 0.1) 
print(ggviolin(Total, x='Cancertype', y='value', fill = 'Cancertype', 
palette = selectedcolorss, xlab = '', ylab = 'Phosphorylation level',
font.y = c(22, "bold", "black"),font.tickslab = c(20, "bold", "black"), x.text.angle=xtextangle,
add = "boxplot", add.params = list(fill="white"),#,width = 0.1
legend='none')+ theme(plot.margin = unit(c(0.1, 0.35, 0.1, 0.35), "in")))#    t、r、b、l分别表示上、右、下、左四侧边距；
dev.off()
}
Total<-c()
medians<-c()
for(i in 1:length(files)){
singles <-get(files[i])
if(any(rownames(singles)==whichsite)){
singles<-data.frame(melt(singles[which(rownames(singles)==whichsite),grep('Solid Tissue Normal$',colnames(singles))]),Cancertype=files[i])
if(length(na.omit(singles[,2])!=0)){
singles[,3]<-paste0(files[i],' (',length(na.omit(singles[,2])),')')
Total<-rbind(Total,singles)
temp<-median(na.omit(singles[,2]))
names(temp)<-paste0(files[i],' (',length(na.omit(singles[,2])),')')
medians<-c(medians,temp)}
	}
}
if(!is.null(Total)){
Total[,3]<-factor(Total[,3],levels=names(medians)[order(medians,decreasing=F)])
selectedcolorss<-colorss[match(gsub(' \\(.*','',levels(Total[,3])),names(colorss))]
names(selectedcolorss)<-levels(Total[,3])

xtextangle = ifelse(length(levels(Total[,3]))>8,30,0)

jpeg(filename=paste0('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/PNGs/PancancerExp/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_Normal.jpeg'),width = 1000, height = 600)
print(ggviolin(Total, x='Cancertype', y='value', fill = 'Cancertype', 
palette = selectedcolorss, xlab = '', ylab = 'Phosphorylation level',
font.y = c(22, "bold", "black"),font.tickslab = c(20, "bold", "black"), x.text.angle=xtextangle,
add = "boxplot", add.params = list(fill="white"),#,width = 0.1
legend='none')+ theme(plot.margin = unit(c(0.1, 0.35, 0.1, 0.35), "in")))#    t、r、b、l分别表示上、右、下、左四侧边距；
dev.off()
		}
	}
}
