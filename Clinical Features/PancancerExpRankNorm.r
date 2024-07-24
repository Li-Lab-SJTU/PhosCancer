files<-c('BRCA','CCRCC','CRC','EOGC','GBM','HCC','HCC_cnhpp', 'HNSCC','LSCC','LUAD','LUAD_cnhpp','OV','PDAC','UCEC' )
BRCA <-read.table('BRCA.tsv',sep='\t',fill=T,header=T,check.names=F)
CCRCC <-read.table('CCRCC.tsv',sep='\t',fill=T,header=T,check.names=F)
CRC <-read.table('CRC.tsv',sep='\t',fill=T,header=T,check.names=F)
EOGC <-read.table('EOGC.tsv',sep='\t',fill=T,header=T,check.names=F)
GBM <-read.table('GBM.tsv',sep='\t',fill=T,header=T,check.names=F)
HCC <-read.table('HCC.tsv',sep='\t',fill=T,header=T,check.names=F)
HCC_cnhpp <-read.table('HCC_cnhpp.tsv',sep='\t',fill=T,header=T,check.names=F)
HNSCC <-read.table('HNSCC.tsv',sep='\t',fill=T,header=T,check.names=F)
LSCC <-read.table('LSCC.tsv',sep='\t',fill=T,header=T,check.names=F)
LUAD <-read.table('LUAD.tsv',sep='\t',fill=T,header=T,check.names=F)
LUAD_cnhpp <-read.table('LUAD_cnhpp.tsv',sep='\t',fill=T,header=T,check.names=F)
OV <-read.table('OV.tsv',sep='\t',fill=T,header=T,check.names=F)
PDAC <-read.table('PDAC.tsv',sep='\t',fill=T,header=T,check.names=F)
UCEC <-read.table('UCEC.tsv',sep='\t',fill=T,header=T,check.names=F)

library(reshape)
library("ggplot2")
library("ggpubr")
library(RColorBrewer)
load('colorss.RData')

Annotation <- read.table("TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
sites<-Annotation[,1]

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
jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/PancancerExpRankNorm/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_Tumor.jpeg'),width = 1000, height = 600)
par(mar = c(5, 7, 4, 2) + 0.1) 
print(ggviolin(Total, x='Cancertype', y='value', fill = 'Cancertype', 
palette = selectedcolorss, xlab = '', ylab = 'Phosphorylation level',
font.y = c(22, "bold", "black"),font.tickslab = c(20, "bold", "black"), x.text.angle=xtextangle,
add = "boxplot", add.params = list(fill="white"),
legend='none')+ theme(plot.margin = unit(c(0.1, 0.35, 0.1, 0.35), "in")))
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

jpeg(filename=paste0('/dssg/home/PhosCancer/Analysis/PNGs/PancancerExpRankNorm/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_Normal.jpeg'),width = 1000, height = 600)
print(ggviolin(Total, x='Cancertype', y='value', fill = 'Cancertype', 
palette = selectedcolorss, xlab = '', ylab = 'Phosphorylation level',
font.y = c(22, "bold", "black"),font.tickslab = c(20, "bold", "black"), x.text.angle=xtextangle,
add = "boxplot", add.params = list(fill="white"),
legend='none')+ theme(plot.margin = unit(c(0.1, 0.35, 0.1, 0.35), "in")))
dev.off()
		}
	}
}
