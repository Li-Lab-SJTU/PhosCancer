library(reshape2)
library(ggplot2)
setwd('E:/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed')

singles <-read.table('LUAD_cnhpp.tsv',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
names(singles)<-gsub('_N$','-Solid Tissue Normal',gsub('_T$','-Primary Tumor',gsub('LADC','LUAD',gsub('_phos','',names(singles)))))
write.table(singles,file="E:/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/LUAD_cnhpp.tsv",sep="\t",quote=F,row.names=F,col.names=T)

singles <-read.table('HCC_cnhpp.tsv',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
names(singles)<-gsub('__cancer','',gsub('_cancer_TiO2_','',gsub('_Human__cancer_TiO2','',gsub('_human__cancer','',gsub('Human__Cancer_','',gsub('_human_liver_cancer','',gsub('.*_L','L',gsub('Liver','',names(singles)))))))))
names(singles)<-gsub('_P$','-Solid Tissue Normal',gsub('_T$','-Primary Tumor',names(singles)))
write.table(singles,file="E:/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/HCC_cnhpp.tsv",sep="\t",quote=F,row.names=F,col.names=T)

singles <-read.table('EOGC.tsv',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
names(singles)<-gsub('-Tumor','-Primary Tumor',gsub('_N','-Solid Tissue Normal',names(singles)))
write.table(singles,file="E:/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/EOGC.tsv",sep="\t",quote=F,row.names=F,col.names=T)


files<-list.files(,pattern='tsv')
for(whichcancer in 1:length(files)){
Phosphoproteome <- read.table(files[whichcancer],sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
rownames(Phosphoproteome)<-paste(Phosphoproteome[,1],Phosphoproteome[,2],sep='|')
Phosphoproteome[Phosphoproteome==0]<-NA

if(any(table(colnames(Phosphoproteome))>1)){ 
unqiuesamples<-unique(colnames(Phosphoproteome)[-c(1:2)])
Totalsamples<-c()
for(i in 1:length(unqiuesamples)){
singlesample<-Phosphoproteome[,colnames(Phosphoproteome)%in%unqiuesamples[i],drop=F]
Totalsamples<-cbind(Totalsamples,rowMeans(singlesample))
	}
colnames(Totalsamples)<-unqiuesamples
}else{Totalsamples<-as.matrix(Phosphoproteome[-c(1:2)])}

Totalsamples<- Totalsamples[,c(grep('Tumor',colnames(Totalsamples)),grep('Normal',colnames(Totalsamples)))]
NormPhosphoGroups<-log2(sweep(Totalsamples,2,apply(Totalsamples,2,median,na.rm=T),FUN="/"))
BoxData<-melt(NormPhosphoGroups)
png(paste0(gsub('tsv','png',files[whichcancer])),width = 800, height = 600)
print(ggplot(BoxData,aes(x=Var2,y=value,fill=Var2))+ 
  geom_boxplot( outlier.size = 0.5)+
  theme_bw()+ theme(axis.text.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  xlab("")+ylab('Phosphorylation'))
dev.off()
write.table(NormPhosphoGroups,file=paste0("E:/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",files[whichcancer]),sep="\t",quote=F,row.names=T,col.names=T)
}


RankPercentile<-function(x){
rank_values <- rank(x, na.last = "keep", ties.method = "min")
total_genes <- sum(!is.na(x))
return(ifelse(!is.na(rank_values), rank_values / total_genes, NA))
}

setwd('E:/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed')
files<-list.files(,pattern='tsv')

for(i in 1:length(files)){
Phosphoproteome <- read.table(files[i],sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)#
NormPhosphoGroups<-apply(Phosphoproteome,2,RankPercentile)
write.table(NormPhosphoGroups,file=paste0("E:/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/PanCancer/",files[i]),sep="\t",quote=F,row.names=T,col.names=T)
}


library(reshape2)
library(ggplot2)
setwd('E:/Projects/PhosCancer/Proteome/Proteome_processed')
singles <-read.table('LUAD_cnhpp.tsv',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
names(singles)<-gsub('_N','-Solid Tissue Normal',gsub('_T','-Primary Tumor',gsub('LADC','LUAD',gsub('_phos','',names(singles)))))
write.table(singles,file="E:/Projects/PhosCancer/Proteome/Proteome_processed/LUAD_cnhpp.tsv",sep="\t",quote=F,row.names=F,col.names=T)
singles <-read.table('HCC_cnhpp.tsv',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
names(singles)<-gsub('__','_',gsub('T$','_T',gsub('P$','_P',gsub('CNHPP_HCC_LC_profiling_','',names(singles)))))
names(singles)<-gsub('_P$','-Solid Tissue Normal',gsub('_T$','-Primary Tumor',names(singles)))
write.table(singles,file="E:/Projects/PhosCancer/Proteome/Proteome_processed/HCC_cnhpp.tsv",sep="\t",quote=F,row.names=F,col.names=T)
singles <-read.table('EOGC.tsv',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
names(singles)<-gsub('-Tumor','-Primary Tumor',gsub('_N','-Solid Tissue Normal',names(singles)))
write.table(singles,file="E:/Projects/PhosCancer/Proteome/Proteome_processed/EOGC.tsv",sep="\t",quote=F,row.names=F,col.names=T)

files<-list.files(,pattern='tsv')
for(whichcancer in 1:length(files)){
Proteome <- read.table(files[whichcancer],sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE,row.names=1)
Proteome[Proteome==0]<-NA
if(any(table(colnames(Proteome))>1)){
unqiuesamples<-unique(colnames(Proteome)[-c(1:2)])
Totalsamples<-c()
for(i in 1:length(unqiuesamples)){
singlesample<-Proteome[,colnames(Proteome)%in%unqiuesamples[i],drop=F]
Totalsamples<-cbind(Totalsamples,rowMeans(singlesample))
	}
colnames(Totalsamples)<-unqiuesamples
}else{Totalsamples<-as.matrix(Proteome)}
Totalsamples<- Totalsamples[,c(grep('Tumor',colnames(Totalsamples)),grep('Normal',colnames(Totalsamples)))]

NormProtein<-log2(sweep(Totalsamples,2,apply(Totalsamples,2,median,na.rm=T),FUN="/"))

BoxData<-melt(NormProtein)
png(paste0(gsub('tsv','png',files[whichcancer])),width = 800, height = 600)
print(ggplot(BoxData,aes(x=Var2,y=value,fill=Var2))+ 
  geom_boxplot( outlier.size = 0.5)+
  theme_bw()+ theme(axis.text.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  xlab("")+ylab('Phosphorylation'))
dev.off()
write.table(NormProtein,file=paste0("E:/Projects/PhosCancer/Proteome/Proteome_processed/",files[whichcancer]),sep="\t",quote=F,row.names=T,col.names=T)
}

setwd('E:/Projects/PhosCancer/Clinical/')
Cancertypes<-gsub('.tsv','',list.files('E:/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed',pattern='tsv'))

for(i in 1:length(Cancertypes)){
Cancertype<-Cancertypes[i]
clinical <- read.csv(paste0("E:/Projects/PhosCancer/Clinical/",Cancertype,"/",grep('PDC_clinical_manifest_',list.files(paste0("E:/Projects/PhosCancer/Clinical/",Cancertype)),value=T)), header = TRUE,fill=T,,check.names=F,stringsAsFactors = FALSE)

clinical[clinical$'Tumor Stage'%in%c('Stage I','Stage1','Stage IA','Stage IB','Stage IC','IA','IB','IC','1A','1B','1C','Stage 1B','IA1','IA2','IA3'),'Tumor Stage']<-'I'
clinical[clinical$'Tumor Stage'%in%c('Stage II','Stage IIA','Stage IIB','Stage IIC','IIA','IIB','IIC'),'Tumor Stage']<-'II'
clinical[clinical$'Tumor Stage'%in%c('Stage III','Stage IIIA','Stage IIIB','Stage IIIC','IIIA','IIIB','IIIC'),'Tumor Stage']<-'III'
clinical[clinical$'Tumor Stage'%in%c('Stage IV','Stage IVA','Stage IVB','Stage IVC','IVA','IVB','IVC'),'Tumor Stage']<-'IV'
print(Cancertype)
print(table( clinical[,'Tumor Stage']))
}

for(i in 1:length(Cancertypes)){
Cancertype<-Cancertypes[i]

clinical <- read.csv(paste0("E:/Projects/PhosCancer/Clinical/",Cancertype,"/",grep('PDC_clinical_manifest_',list.files(paste0("E:/Projects/PhosCancer/Clinical/",Cancertype)),value=T)), header = TRUE,fill=T,,check.names=F,stringsAsFactors = FALSE)
print(Cancertype)
clinical[clinical$'Tumor Stage'%in%c('Stage I','Stage1','Stage IA','Stage IB','Stage IC','IA','IB','IC','1A','1B','1C','Stage 1B','IA1','IA2','IA3'),'Tumor Stage']<-'I'
clinical[clinical$'Tumor Stage'%in%c('Stage II','Stage IIA','Stage IIB','Stage IIC','IIA','IIB','IIC'),'Tumor Stage']<-'II'
clinical[clinical$'Tumor Stage'%in%c('Stage III','Stage IIIA','Stage IIIB','Stage IIIC','IIIA','IIIB','IIIC'),'Tumor Stage']<-'III'
clinical[clinical$'Tumor Stage'%in%c('Stage IV','Stage IVA','Stage IVB','Stage IVC','IVA','IVB','IVC'),'Tumor Stage']<-'IV'
clinical[!clinical$'Tumor Stage'%in%c('I','II','III','IV'),'Tumor Stage']<-NA
print(table( clinical[,'Tumor Stage']))
if(!Cancertype%in%c('HCC_cnhpp','LUAD_cnhpp')){
clinical[clinical[,'Vital Status']=='Alive'&clinical[,'Days to Last Follow Up']!='null','Days to Death']<-clinical[clinical[,'Vital Status']=='Alive'&clinical[,'Days to Last Follow Up']!='null','Days to Last Follow Up']
clinical[,'Days to Death']<-as.numeric(clinical[,'Days to Death'])
clinical[clinical[,'Vital Status']=='Alive','Vital Status']<-0
clinical[clinical[,'Vital Status']=='Dead','Vital Status']<-1
clinical[,'Vital Status']<-as.numeric(clinical[,'Vital Status'])

clinical[clinical[,'Progression or Recurrence']%in%c('No','no')&clinical[,'Vital Status']==1&!is.na(clinical[,'Days to Death']),'Progression or Recurrence']<-'Yes'
clinical[clinical[,'Progression or Recurrence']%in%c('No','no')&clinical[,'Vital Status']==1&!is.na(clinical[,'Days to Death']),'Days to Recurrence']<-clinical[clinical[,'Progression or Recurrence']%in%c('No','no')&clinical[,'Vital Status']==1&!is.na(clinical[,'Days to Death']),'Days to Death']

if(Cancertype=='GBM'){
clinical[which(clinical[,'Progression or Recurrence']%in%c('No','no')&(as.numeric(clinical[,'Days to Recurrence'])>as.numeric(clinical[,'Days to Last Follow Up']))),'Days to Last Follow Up']<-clinical[which(clinical[,'Progression or Recurrence']%in%c('No','no')&(as.numeric(clinical[,'Days to Recurrence'])>as.numeric(clinical[,'Days to Last Follow Up']))),'Days to Recurrence']
clinical[clinical[,'Progression or Recurrence']%in%c('No','no')&clinical[,'Days to Last Follow Up']!='null',c('Days to Recurrence')]<-clinical[clinical[,'Progression or Recurrence']%in%c('No','no')&clinical[,'Days to Last Follow Up']!='null',c('Days to Last Follow Up')]
	
}
else{
clinical[clinical[,'Progression or Recurrence']%in%c('No','no')&clinical[,'Days to Last Follow Up']!='null','Days to Recurrence']<-clinical[clinical[,'Progression or Recurrence']%in%c('No','no')&clinical[,'Days to Last Follow Up']!='null','Days to Last Follow Up']
}

clinical[,'Days to Recurrence']<-as.numeric(clinical[,'Days to Recurrence'])
clinical[clinical[,'Progression or Recurrence']%in%c('No','no'),'Progression or Recurrence']<-0
clinical[clinical[,'Progression or Recurrence']%in%c('Yes','yes'),'Progression or Recurrence']<-1
clinical[,'Progression or Recurrence']<-as.numeric(clinical[,'Progression or Recurrence'])
}
write.table(clinical,file=paste0("E:/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,'.tsv'),sep='\t',row.names=F,col.names=T)			
}


setwd('E:/Projects/PhosCancer/Phosphoproteome/CorrectedwithProtein')
Cancertypes<-gsub('.tsv','',list.files('E:/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed',pattern='tsv'))

for(i in 1:length(Cancertypes)){
Cancertype<-Cancertypes[i]
print(Cancertype)
Phospho<-read.table(paste0('E:/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/',Cancertype,'.tsv'),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)

Protein<-read.table(paste0('E:/Projects/PhosCancer/Proteome/Proteome_processed/',Cancertype,'.tsv'),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)

interS<-intersect(colnames(Phospho),colnames(Protein))
Phospho<-Phospho[,match(interS,colnames(Phospho))]
Protein<-Protein[,match(interS,colnames(Protein))]
rownames(Protein)<-unlist(strsplit(rownames(Protein),'\\|'))[grep('sp',unlist(strsplit(rownames(Protein),'\\|')))+1]

TotalNormlized<-c()
for(j in 1:nrow(Phospho)){
BaseProtein<-unlist(strsplit(rownames(Phospho)[j],'\\|'))[grep('sp',unlist(strsplit(rownames(Phospho)[j],'\\|')))+1]

if(any(rownames(Protein)==BaseProtein)){

Data<-data.frame(cbind(t(Phospho[j,]),Protein=t(Protein[rownames(Protein)==BaseProtein,])))
colnames(Data)<-c('site','protein')
Residuals<-residuals(lm(site ~ protein, data=Data))

Normlized<-matrix(Residuals[match(interS,names(Residuals))],nrow=1)
rownames(Normlized)<-rownames(Phospho)[j]
colnames(Normlized)<-interS
TotalNormlized<-rbind(TotalNormlized,Normlized)
		}
write.table(TotalNormlized,file=paste0(Cancertypes[i],'.tsv'),sep="\t",quote=F,row.names=T,col.names=T)
	}
}
