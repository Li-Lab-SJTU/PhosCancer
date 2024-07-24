Cancertype='UCEC'# BRCA CCRCC CRC GBM OV PDAC UCEC     EOGC HNSCC LSCC LUAD PAAD HCC 
Annotation <- read.table("/dssg/home/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
clinical <- read.table(paste0("/dssg/home/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,".tsv"), header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE,sep='\t')
library("survival") 
library("survminer")
NormPhosphoGroups <- read.table(paste0("/dssg/home/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
PhosNAcutoff=0.80
ProteinNAcutoff=0.50
NormPhosphoGroupsSur<-data.frame(clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),clinical$'Cases Submitter ID'),c('Days to Recurrence','Progression or Recurrence','Days to Death','Vital Status','Days to Last Follow Up')],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)
NormPhosphoGroupsfilteredraw<-NormPhosphoGroupsSur[,-c(which(apply(NormPhosphoGroupsSur[,-c(1:5)],2,function(x){length(which(is.na(x)))})>nrow(NormPhosphoGroupsSur[,-c(1:5)])*PhosNAcutoff)+5)]
dim(NormPhosphoGroupsfilteredraw)
names(NormPhosphoGroupsfilteredraw)[1:5]<-gsub(' ','',names(NormPhosphoGroupsfilteredraw)[1:5])

sites<-colnames(NormPhosphoGroupsfilteredraw)[6:ncol(NormPhosphoGroupsfilteredraw)]
if(length(sites)>0){
for(nsite in 1:length(sites)){
print(nsite)
whichsite<-sites[nsite]
NormPhosphoGroupsfiltered<-NormPhosphoGroupsfilteredraw[,c(1:5,which(colnames(NormPhosphoGroupsfilteredraw)==whichsite))]
names(NormPhosphoGroupsfiltered)[6]<-'site'
covariates <-names(NormPhosphoGroupsfiltered)[6:ncol(NormPhosphoGroupsfiltered)]

NormPhosphoGroupsfilteredOS<-NormPhosphoGroupsfiltered
splots<-list()
if(length(which(na.omit(NormPhosphoGroupsfiltered[,c(3,4,6)])[,2]==1))>1){
OSres.cut <- surv_cutpoint(NormPhosphoGroupsfilteredOS, time = "DaystoDeath", event = "VitalStatus",variables = covariates)
OSres.cut <- surv_categorize(OSres.cut)

KMres<-c()
for(i in 3:ncol(OSres.cut )){
fit <- survfit(Surv(DaystoDeath, VitalStatus) ~get(names(OSres.cut)[i]) , data = OSres.cut )
KMres<-c(KMres,surv_pvalue(fit)[1,2])
}
KMrescutpoint<-data.frame(KM.cutpoint.P=KMres)
rownames(KMrescutpoint)<-names(OSres.cut )[3:ncol(OSres.cut)]
dim(KMrescutpoint)

fit <- survfit(Surv(DaystoDeath, VitalStatus) ~site , data = OSres.cut )

splots[[1]]=ggsurvplot(fit,data = OSres.cut,size = 2.0,
                palette = c('#ED0000','#183F7D'),
                conf.int = F,conf.int.style='ribbon', conf.int.alpha=0.1,
                pval = T,pval.method = T,censor.shape="+",pval.size=8,pval.method.size=8,
                risk.table = T,risk.table.pos='in', risk.table.fontsize=6,
				surv.scale='percent',
                legend.title="",
				legend = 'right',
				legend.labs=c('high','low'),
				tables.theme=theme_cleantable() ,
                title="Overall survival curve (cutpoint)", 
				xlab='Survival time (Day)',
                ggtheme = theme_classic2() ,
				censor.size=10,
				font.main = c(18, "bold", "black"),
                font.x = c(16, "bold", "black"),
                font.y = c(16, "bold", "black"),
                font.tickslab = c(14, "black"),
                font.title  = c(16, "bold", "black"),
                font.legend  = c(16, "bold", "black")
 )


for(i in 6:ncol(NormPhosphoGroupsfilteredOS)){
high<-which(NormPhosphoGroupsfilteredOS[,i]>median(na.omit(NormPhosphoGroupsfilteredOS[,i])))
low<-which(NormPhosphoGroupsfilteredOS[,i]<=median(na.omit(NormPhosphoGroupsfilteredOS[,i])))
NormPhosphoGroupsfilteredOS[high,i]<-'High'
NormPhosphoGroupsfilteredOS[low,i]<-'Low'
}
NormPhosphoGroupsfilteredOS<-NormPhosphoGroupsfilteredOS[,-c(1,2,5)]
KMres<-c()
for(i in 3:ncol(NormPhosphoGroupsfilteredOS)){
fit <- survfit(Surv(DaystoDeath, VitalStatus) ~get(names(NormPhosphoGroupsfilteredOS)[i]) , data = NormPhosphoGroupsfilteredOS)
KMres<-c(KMres,surv_pvalue(fit)[1,2])
}
KMresmedian<-data.frame(KM.Median.P=KMres)
rownames(KMresmedian)<-names(NormPhosphoGroupsfilteredOS)[3:ncol(NormPhosphoGroupsfilteredOS)]
dim(KMresmedian)

fit <- survfit(Surv(DaystoDeath, VitalStatus) ~site , data = NormPhosphoGroupsfilteredOS)

splots[[2]]=ggsurvplot(fit,data = NormPhosphoGroupsfilteredOS,size = 2.0,
                palette = c('#ED0000','#183F7D'),
                conf.int = F,conf.int.style='ribbon', conf.int.alpha=0.1,
                pval = T,pval.method = T,censor.shape="+",pval.size=8,pval.method.size=8,
                risk.table = T,risk.table.pos='in', risk.table.fontsize=6,
				surv.scale='percent',
                legend.title="",
				legend = 'right',
				legend.labs=c('high','low'),
				tables.theme=theme_cleantable() ,
                title="Overall survival curve (median)", 
				xlab='Survival time (Day)',
                ggtheme = theme_classic2() ,
				censor.size=10,
				font.main = c(18, "bold", "black"),
                font.x = c(16, "bold", "black"),
                font.y = c(16, "bold", "black"),
                font.tickslab = c(14, "black"),
                font.title  = c(16, "bold", "black"),
                font.legend  = c(16, "bold", "black")
 )
}
NormPhosphoGroupsfilteredDFS<-NormPhosphoGroupsfiltered
if(length(which(na.omit(NormPhosphoGroupsfiltered[,c(1,2,6)])[,2]==1))>1){
DFSres.cut <- surv_cutpoint(NormPhosphoGroupsfilteredDFS, time = "DaystoRecurrence", event = "ProgressionorRecurrence",variables = covariates)
DFSres.cat <- surv_categorize(DFSres.cut)

KMres<-c()
for(i in 3:ncol(DFSres.cat )){
fit <- survfit(Surv(DaystoRecurrence, ProgressionorRecurrence) ~get(names(DFSres.cat)[i]) , data = DFSres.cat )
KMres<-c(KMres,surv_pvalue(fit)[1,2])
}
KMrescutpoint<-data.frame(KM.cutpoint.P=KMres)
rownames(KMrescutpoint)<-names(DFSres.cat )[3:ncol(DFSres.cat )]
dim(KMrescutpoint)

fit <- survfit(Surv(DaystoRecurrence, ProgressionorRecurrence) ~site , data = DFSres.cat )

splots[[3]]=ggsurvplot(fit,data = DFSres.cat,size = 2.0,
                palette = c('#ED0000','#183F7D'),
                conf.int = F,conf.int.style='ribbon', conf.int.alpha=0.1,
                pval = T,pval.method = T,censor.shape="+",pval.size=8,pval.method.size=8,
                risk.table = T,risk.table.pos='in', risk.table.fontsize=6,
				surv.scale='percent',
                legend.title="",
				legend = 'right',
				legend.labs=c('high','low'),
				tables.theme=theme_cleantable() ,
                title="Disease Free survival curve (cutpoint)", 
				xlab='Survival time (Day)',
                ggtheme = theme_classic2() ,
				censor.size=10,
				font.main = c(18, "bold", "black"),
                font.x = c(16, "bold", "black"),
                font.y = c(16, "bold", "black"),
                font.tickslab = c(14, "black"),
                font.title  = c(16, "bold", "black"),
                font.legend  = c(16, "bold", "black")
 )

NormPhosphoGroupsfilteredDFS<-NormPhosphoGroupsfilteredDFS[-union(which(is.na(NormPhosphoGroupsfilteredDFS[,'DaystoRecurrence'])),which(is.na(NormPhosphoGroupsfilteredDFS[,'ProgressionorRecurrence']))),]

for(i in 6:ncol(NormPhosphoGroupsfilteredDFS)){
high<-which(NormPhosphoGroupsfilteredDFS[,i]>median(na.omit(NormPhosphoGroupsfilteredDFS[,i])))
low<-which(NormPhosphoGroupsfilteredDFS[,i]<=median(na.omit(NormPhosphoGroupsfilteredDFS[,i])))
NormPhosphoGroupsfilteredDFS[high,i]<-'High'
NormPhosphoGroupsfilteredDFS[low,i]<-'Low'
}
NormPhosphoGroupsfilteredDFS<-NormPhosphoGroupsfilteredDFS[,-c(3,4,5)]
KMres<-c()
for(i in 3:ncol(NormPhosphoGroupsfilteredDFS)){
fit <- survfit(Surv(DaystoRecurrence, ProgressionorRecurrence) ~get(names(NormPhosphoGroupsfilteredDFS)[i]) , data = NormPhosphoGroupsfilteredDFS)
KMres<-c(KMres,surv_pvalue(fit)[1,2])
}
KMresmedian<-data.frame(KM.Median.P=KMres)
rownames(KMresmedian)<-names(NormPhosphoGroupsfilteredDFS)[3:ncol(NormPhosphoGroupsfilteredDFS)]
dim(KMresmedian)

fit <- survfit(Surv(DaystoRecurrence, ProgressionorRecurrence) ~site , data = NormPhosphoGroupsfilteredDFS)

splots[[4]]=ggsurvplot(fit,data = NormPhosphoGroupsfilteredDFS,size = 2.0,
                palette = c('#ED0000','#183F7D'),
                conf.int = F,conf.int.style='ribbon', conf.int.alpha=0.1,
                pval = T,pval.method = T,censor.shape="+",pval.size=8,pval.method.size=8,
                risk.table = T,risk.table.pos='in', risk.table.fontsize=6,
				surv.scale='percent',
                legend.title="",
				legend = 'right',
				legend.labs=c('high','low'),
				tables.theme=theme_cleantable(risk.table.height = 0.5) ,
                tables.height = 0.5,risk.table.height=0.5,
				title="Disease Free survival curve (median)", 
				xlab='Survival time (Day)',
                ggtheme = theme_classic2(base_size=16) ,
				censor.size=10,
				font.main = c(18, "bold", "black"),
                font.x = c(16, "bold", "black"),
                font.y = c(16, "bold", "black"),
                font.tickslab = c(14, "black"),
                font.title  = c(16, "bold", "black"),
                font.caption  = c(16, "bold", "black"),
                font.legend  = c(16, "bold", "black")# Í¼Àý×ÖÌå
 )
}
if(length(which(na.omit(NormPhosphoGroupsfiltered[,c(3,4,6)])[,2]==1))>1&length(which(na.omit(NormPhosphoGroupsfiltered[,c(1,2,6)])[,2]==1))>1){
 jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/Surdiffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_',Cancertype,'.jpeg'),width = 1100, height = 800)
 print(arrange_ggsurvplots(splots, print = TRUE,ncol = 2, nrow = 2))
 dev.off()} 
if(length(which(na.omit(NormPhosphoGroupsfiltered[,c(3,4,6)])[,2]==1))>1&length(which(na.omit(NormPhosphoGroupsfiltered[,c(1,2,6)])[,2]==1))<=1){
 jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/Surdiffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_',Cancertype,'.jpeg'),width = 1100, height = 400)
 print(arrange_ggsurvplots(splots[1:2], print = TRUE,ncol = 2, nrow = 1))
 dev.off()} 
if(length(which(na.omit(NormPhosphoGroupsfiltered[,c(3,4,6)])[,2]==1))<=1&length(which(na.omit(NormPhosphoGroupsfiltered[,c(1,2,6)])[,2]==1))>1){
 jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/Surdiffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_',Cancertype,'.jpeg'),width = 1100, height = 400)
 print(arrange_ggsurvplots(splots[3:4], print = TRUE,ncol = 2, nrow = 1))
 dev.off()}
 }
}
