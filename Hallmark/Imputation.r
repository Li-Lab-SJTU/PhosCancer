files<-list.files(,pattern='tsv')
ProteinNAcutoff=0.50
SwissProt <- read.table("uniprotkb_reviewed_true_AND_model_organ_2024_01_18.fasta",comment.char="",quote='',sep='\t',header = F,fill=T,check.names=F,stringsAsFactors = FALSE)
SwissProt<-SwissProt[intersect(grep('>sp',SwissProt[,1]),grep("GN=",SwissProt[,1])),]

MatchID<-data.frame('UniProt ID'=unlist(strsplit(SwissProt,split='\\|'))[grep('HUMAN',unlist(strsplit(SwissProt,split='\\|')))-1],GN=gsub('GN=','',grep("GN=",unlist(strsplit(SwissProt,split=" ")),value = TRUE)))

for(i in 1:length(files)){
Cancertype<-gsub('\\.tsv','',files[i])
NormPhosphoGroups <- read.table(paste0(Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)

NormProtein <- read.table(paste0(Cancertype,".tsv"),sep='\t', header = TRUE,fill=T,,check.names=F,stringsAsFactors = FALSE)
intersamples<-intersect(grep('Primary Tumor$',names(NormPhosphoGroups),value=T),grep('Primary Tumor$',names(NormProtein),value=T))

NormProteinfiltered<-NormProtein[,match(intersamples,names(NormProtein))]
NormProteinfiltered<-NormProteinfiltered[-which(apply(NormProteinfiltered,1,function(x){length(which(is.na(x)))})>ncol(NormProteinfiltered)*ProteinNAcutoff),]

library(impute)
NormProteinfilteredImpute<-impute.knn(as.matrix(NormProteinfiltered),
  k=5,  
  rowmax=0.5, 
  colmax=0.5)[[1]]

rownames(NormProteinfilteredImpute)<-MatchID[match(unlist(strsplit(rownames(NormProteinfilteredImpute),'\\|'))[grep('_HUMAN',unlist(strsplit(rownames(NormProteinfilteredImpute),'\\|')))-1],MatchID[,1]),2]
 
library(GSVA)    
library(clusterProfiler)
gmt=read.gmt("h.all.v7.5.1.symbols.gmt") 
Hallmark.rename <- read.table('Hallmark.rename.tsv',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE, encoding = "latin1")
gmt$term<-Hallmark.rename[match(gmt$term,Hallmark.rename[,1]),2]
genelist = split(gmt$gene, gmt$term)
gsva_matrix <- gsva(as.matrix(NormProteinfilteredImpute), genelist, kcdf="Gaussian",method = "ssgsea",min.sz=10)
write.table(gsva_matrix,file=paste0("E:/Proteome_gsva/",Cancertype,'.tsv'),sep="\t",quote=F,row.names=T,col.names=T)
}
