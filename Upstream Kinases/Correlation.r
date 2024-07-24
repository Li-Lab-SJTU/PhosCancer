Cancertype<-'UCEC'

KinaseGenes <-read.table("/dssg/home/Projects/PhosCancer/Analysis/Quantitative/KinaseCor/KinHub.tsv",sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE,check.names=F)#[,-c(1,2,6)])
KinaseGenes<-KinaseGenes[KinaseGenes[,3]!='',]
KinaseGenes<-KinaseGenes[-grep('Domain',KinaseGenes[,2]),]
dim(KinaseGenes)

PhosNAcutoff=0.80
NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormProteomeGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Proteome/Proteome_processed/",Cancertype,".tsv"),comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
rownames(NormProteomeGroups)<-unlist(strsplit(rownames(NormProteomeGroups),'\\|'))[grep('_HUMAN',unlist(strsplit(rownames(NormProteomeGroups),'\\|')))-1]
KinaseExp<-NormProteomeGroups[rownames(NormProteomeGroups)%in%KinaseGenes[,8],]

intersamples<-intersect(grep('Primary Tumor$',names(NormPhosphoGroups),value=T),grep('Primary Tumor$',names(KinaseExp),value=T))

NormPhosphoGroupsfiltered<-NormPhosphoGroups[,match(intersamples,names(NormPhosphoGroups))]
NormPhosphoGroupsfiltered<-NormPhosphoGroupsfiltered[-which(apply(NormPhosphoGroupsfiltered,1,function(x){length(which(is.na(x)))})>ncol(NormPhosphoGroupsfiltered)*PhosNAcutoff),]

KinaseExp<-KinaseExp[,match(intersamples,names(KinaseExp))]

sites<-rownames(NormPhosphoGroupsfiltered)

KinaseCorTop<-c()
for(i in 1:length(sites)){
print(i)
whichsite=sites[i]

KinaseExp1<-KinaseExp[apply(KinaseExp[,!is.na(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,]))],1,function(x){length(which(!is.na(x)))})>10,]

if(nrow(KinaseExp1)!=0){
corP<-apply(KinaseExp1,1,function(x){cor.test(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,]),x,method='spearman')$p.value})
corE<-apply(KinaseExp1,1,function(x){cor.test(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,]),x,method='spearman')$estimate})

SpearmanP=corP
SpearmanE=corE
PearsonP=apply(KinaseExp1,1,function(x){cor.test(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,]),x,method='pearson')$p.value})
PearsonE=apply(KinaseExp1,1,function(x){cor.test(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,]),x,method='pearson')$estimate})

whichKinase<-names(SpearmanE[which(SpearmanP<0.05)][order(abs(SpearmanE[which(SpearmanP<0.05)]),decreasing=T)][1:10])
tempTop<-data.frame(Sites=rownames(NormPhosphoGroupsfiltered)[i],TumorSampleSize=apply(KinaseExp1[match(whichKinase,rownames(KinaseExp1)),!is.na(as.numeric(NormPhosphoGroupsfiltered[rownames(NormPhosphoGroupsfiltered)==whichsite,]))],1,function(x){length(which(!is.na(x)))}),
MostKinase=whichKinase,
SpearmanP=SpearmanP[match(whichKinase,names(SpearmanP))],
SpearmanFDR=p.adjust(SpearmanP,'BH')[match(whichKinase,names(SpearmanP))],
SpearmanE=SpearmanE[match(whichKinase,names(SpearmanE))],
PearsonP=PearsonP[match(whichKinase,names(PearsonP))],
PearsonFDR=p.adjust(PearsonP,'BH')[match(whichKinase,names(PearsonP))],
PearsonE=PearsonE[match(whichKinase,names(PearsonE))])
KinaseCorTop<-rbind(tempTop,KinaseCorTop)
	}
}
write.table(KinaseCorTop,file=paste0("/dssg/home/Projects/PhosCancer/Analysis/Quantitative/KinaseCor/",Cancertype,"_KinaseCorTop.tsv"),sep="\t",quote=F,row.names=F,col.names=T)
