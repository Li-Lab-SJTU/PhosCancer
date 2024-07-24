PhosNAcutoff=0.80
ProteinNAcutoff=0.50

readtable<-function(Cancertype){
NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
gsva_matrix <- read.table(paste0("/dssg/home/Projects/PhosCancer/Proteome/Proteome_gsva/",Cancertype,".tsv"),comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
#
intersamples<-intersect(grep('Primary Tumor$',names(NormPhosphoGroups),value=T),grep('Primary Tumor$',names(gsva_matrix),value=T))

NormPhosphoGroupsfiltered<-NormPhosphoGroups[,match(intersamples,names(NormPhosphoGroups))]
NormPhosphoGroupsfiltered<-NormPhosphoGroupsfiltered[-which(apply(NormPhosphoGroupsfiltered,1,function(x){length(which(is.na(x)))})>ncol(NormPhosphoGroupsfiltered)*PhosNAcutoff),]

all(names(NormPhosphoGroupsfiltered)==names(gsva_matrix))
return(list(NormPhosphoGroupsfiltered,gsva_matrix))
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

Hallmark.rename <- read.table('/dssg/home/Projects/PhosCancer/Analysis/Hallmark.rename.tsv',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE, encoding = "latin1")

Annotation <- read.table("/dssg/home/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
sites<-Annotation[,1]

library(reshape)
if(length(sites)>0){
for(site in 1:length(sites)){
print(site)
whichsite=sites[site]
Total<-c()
for(i in 1:length(files)){
singles <-get(files[i])
Phos<-singles[[1]]
Pathway<-singles[[2]]
if(length(which(rownames(Phos)==whichsite))!=0){
corP<-apply(Pathway,1,function(x){cor.test(as.numeric(Phos[rownames(Phos)==whichsite,]),x,method='spearman')$p.value})
corE<-apply(Pathway,1,function(x){cor.test(as.numeric(Phos[rownames(Phos)==whichsite,]),x,method='spearman')$estimate})

temp<-cbind(corP,corE)
colnames(temp)<-c(paste(files[i],'pvalue'),files[i])
temp<-temp[match(Hallmark.rename[,2],rownames(temp)),]
rownames(temp)<-Hallmark.rename[,2]
Total<-cbind(Total,temp)
	}
}
if(!is.null(Total)){
if(any(apply(Total,1,function(x){length(which(is.na(x)))})==ncol(Total))){
 Total<-Total[-which(apply(Total,1,function(x){length(which(is.na(x)))})==ncol(Total)),]
}
 Total<-as.data.frame(Total)
 library(pheatmap)
 library(RColorBrewer) 
 Total[,grep('pvalue',colnames(Total))]<-ifelse(Total[,grep('pvalue',colnames(Total))]<0.05,'*','')
 paletteLength <- 100
 myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(paletteLength)
 myBreaks <- c(seq(min(na.omit(Total[,-grep('pvalue',colnames(Total))])), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(na.omit(Total[,-grep('pvalue',colnames(Total))]))/paletteLength, max(na.omit(Total[,-grep('pvalue',colnames(Total))])), length.out=floor(paletteLength/2)))
jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/HallmarkCor/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_ALL.jpeg'),width = 1160, height = nrow(Total)*8)

 if(length(grep('pvalue',colnames(Total)))>1){
  p2<-pheatmap(t(Total[,-grep('pvalue',colnames(Total)),drop=F]),gap_row=c(1,2),show_colnames = T,labels_col=sapply(1:nrow(Total), function(i) {as.expression(bquote(bold(.(rownames(Total)[i]))))}),labels_row=sapply(1:(ncol(Total)/2), function(i) {as.expression(bquote(bold(.(colnames(Total)[-grep('pvalue',colnames(Total))][i]))))}),scale = "none",cluster_row = T, cluster_col = T, border=NA,display_numbers = t(Total[,grep('pvalue',colnames(Total)),drop=F]),fontsize_number = 15, fontsize_row = 15, fontsize_col=13,number_color = "white",col=colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100), breaks=myBreaks,angle_col = 45,border_color=NA,main = "Spearman’s correlation coefficient in pan-cancer",fontsize = 16)}else{
  p2<-pheatmap(t(Total[,-grep('pvalue',colnames(Total)),drop=F]),gap_row=c(1,2),show_colnames = T,labels_col=sapply(1:nrow(Total), function(i) {as.expression(bquote(bold(.(rownames(Total)[i]))))}),labels_row=sapply(1:(ncol(Total)/2), function(i) {as.expression(bquote(bold(.(colnames(Total)[-grep('pvalue',colnames(Total))][i]))))}),scale = "none",cluster_row = F, cluster_col = F, border=NA,display_numbers = t(Total[,grep('pvalue',colnames(Total)),drop=F]),fontsize_number = 15, fontsize_row = 15,fontsize_col=13,number_color = "white",col=colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100), breaks=myBreaks,angle_col = 45,border_color=NA,main = "Spearman’s correlation coefficients in pan-cancer",fontsize = 16)
  }
  if(length(grep('pvalue',colnames(Total)))>1){
  p1<-pheatmap(matrix(1:(ncol(Total)/2),ncol(Total)/2,1),col='white', na_color='white',cluster_row = F,cluster_col = F,legend=F,border_color=NA)
   }else{p1<-pheatmap(matrix(1:2,2,1),col='white', na_color='white',cluster_row = F,cluster_col = F,legend=F,border_color=NA)
}
 print(cowplot::plot_grid(p1$gtable, p2$gtable, ncol= 2, rel_widths = c(1, 10)))
 dev.off()
 }
  }
}


