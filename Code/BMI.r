Cancertype<-'CCRCC'
PhosNAcutoff=0.80
Annotation <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
library(readxl)
clinical<-as.data.frame(read_xlsx(paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Clinical/FromPaper/",Cancertype,".xlsx"),sheet = 2))
NormPhosphoGroups <- read.table(paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroupsClinical<-data.frame(clinical[match(gsub('-Primary Tumor','',grep('Primary Tumor$',names(NormPhosphoGroups),value=T)),gsub('X','',clinical$'Case_ID')),'BMI',drop=F],t(NormPhosphoGroups[,grep('Primary Tumor$',names(NormPhosphoGroups))]),check.names=F)

#预处理

####Age
#缺失值过滤
NormPhosphoGroupsfiltered<-NormPhosphoGroupsClinical[,-c(which(apply(NormPhosphoGroupsClinical[,-1],2,function(x){length(which(is.na(x)))})>nrow(NormPhosphoGroupsClinical[,-1])*PhosNAcutoff)+1)]
dim(NormPhosphoGroupsfiltered)

Agediffer<-c()
for(i in 2:ncol(NormPhosphoGroupsfiltered)){ 

SiteData<-na.omit(NormPhosphoGroupsfiltered[,c(1,i)])

#计算
TumorSampleSize<-nrow(SiteData)
SpearmanP<-cor.test(SiteData[,1],SiteData[,2],method='spearman')$p.value
SpearmanE<-cor.test(SiteData[,1],SiteData[,2],method='spearman')$estimate
Agediffer<-rbind(Agediffer,data.frame(TumorSampleSize,SpearmanP,SpearmanE))

#画图
library(ggdist)        ## halfeye plots
library(gghalves)      ## off-set jitter
library(ggplot2)
library("ggpubr")
whichsite<-colnames(SiteData)[2]
colnames(SiteData)<-c('BMI','Site')
#${UniProtID}_${Position}_${cancerType}.jpeg
jpeg(filename=paste0('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/PNGs/BMICor/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_',Cancertype,'.jpeg'),width = 500, height = 500)
print(
ggscatter(SiteData, x = "BMI", y = "Site",xlab="BMI",ylab='Phosphorylation level',
          size = 2,  shape=8,                       # 点的大小为1.5
          add = "reg.line",                   # 添加回归直线
          conf.int = TRUE,                    # 添加置信区间
          rug = TRUE ,                      # 添加轴须线
		  color='#0054A6',
          add.params = list(color = "black",  # 设置回归直线、置信域的颜色
                          fill = "#0054A6"),
          ggtheme = theme_bw()
)+
  stat_cor(size=8,fontface = "bold",method = "spearman",r.accuracy=0.01,p.accuracy=0.01) 
+
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5,
								  face = "bold"),
        axis.title = element_text(size = 19, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    # vjust = 1.9, 
                                    # hjust = 0.5, 
                                    # angle = 90
									),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 18,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 19,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  ))
dev.off()
}
rownames(Agediffer)<-colnames(NormPhosphoGroupsfiltered)[2:ncol(NormPhosphoGroupsfiltered)]
Agediffer$SpearmanFDR<-p.adjust(Agediffer$SpearmanP,'BH')
write.table(Agediffer,file=paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/Quantitative/BMICor/",Cancertype,"_BMICor.tsv"),sep="\t",quote=F,row.names=T,col.names=T)
