Cancertype<-'BRCA'
PhosNAcutoff=0.80
Annotation <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
# clinical <- read.table(paste0("E:/PhD/LiLab/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,".tsv"), header = TRUE,fill=T,sep='\t',check.names=F,stringsAsFactors = FALSE)
library(readxl)
clinical <-read.table(paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Clinical/Clinical_processed/",Cancertype,".tsv"), header = TRUE,fill=T,sep='\t',check.names=F,stringsAsFactors = FALSE)
NormPhosphoGroups <- read.table(paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Phosphoproteome/Phosphoproteome_processed/",Cancertype,".tsv"),sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
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

####
#缺失值过滤
NormPhosphoGroupsClinical[,1]
if(any(NormPhosphoGroupsClinical[,1]=='Not Reported')){NormPhosphoGroupsClinical<-NormPhosphoGroupsClinical[-which(NormPhosphoGroupsClinical[,1]=='Not Reported'),]}
NormPhosphoGroupsClinical[,1]

if(length(unique(NormPhosphoGroupsClinical[,1]))>1){
NormPhosphoGroupsfiltered<-NormPhosphoGroupsClinical[,-which(apply(NormPhosphoGroupsClinical[NormPhosphoGroupsClinical[,1]=='Male',],2,function(x){length(which(is.na(x)))})>length(which(NormPhosphoGroupsClinical[,1]=='Male'))*PhosNAcutoff|apply(NormPhosphoGroupsClinical[NormPhosphoGroupsClinical[,1]=='Female',],2,function(x){length(which(is.na(x)))})>length(which(NormPhosphoGroupsClinical[,1]=='Female'))*PhosNAcutoff)]

Agediffer<-c()
for(i in 2:ncol(NormPhosphoGroupsfiltered)){ 

SiteData<-na.omit(NormPhosphoGroupsfiltered[,c(1,i)])
Control<-which(SiteData[,1]=='Female')
Case<-which(SiteData[,1]=='Male')

#计算
TumorSampleSize<-nrow(SiteData)
FemaleSampleSize<-length(Control)
MaleSampleSize<-length(Case)

FemaleSampleMean<-mean(2^SiteData[Control,2])#非NA均值
MaleSampleMean<-mean(2^SiteData[Case,2])#非NA均值
log2FC<-log2(mean(2^SiteData[Case,2])/mean(2^SiteData[Control,2]))
wilcoxP<-wilcox.test(SiteData[Case,2],SiteData[Control,2])$p.value
Agediffer<-rbind(Agediffer,data.frame(TumorSampleSize,FemaleSampleSize,MaleSampleSize,FemaleSampleMean,MaleSampleMean,log2FC,wilcoxP))

#画图
library(ggdist)        ## halfeye plots
library(gghalves)      ## off-set jitter
library(ggplot2)
library("ggpubr")
whichsite<-colnames(SiteData)[2]
boxplotData<-SiteData
boxplotData[Control,1]<-'Control'
boxplotData[Case,1]<-'Case'
boxplotData<-data.frame(Exp=boxplotData[,2],Type=boxplotData[,1])
#${UniProtID}_${Position}_${cancerType}.jpeg
jpeg(filename=paste0('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/PNGs/Genderdiffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_',Cancertype,'.jpeg'),width = 700, height = 500)
print(ggplot(boxplotData, aes(x = Type, y = Exp, color = Type, fill = Type)) +
  scale_color_manual(limits=c("Control","Case"), values=c("#5F80B4","#922927"))+
  scale_fill_manual(limits=c("Control","Case"), values=c("#5F80B4","#922927")) + 
  ggdist::stat_halfeye(  trim = F,
    adjust = 0.5, ## bandwidth 调整 adjust 参数来增加带宽，以获得更平滑的曲线较小的带宽会产生更尖锐的曲线，而较大的带宽会产生更平滑的曲线。此外，还可以尝试调整 width 参数来控制半眼图的宽度，以及调整 position 参数来微调曲线的位置。
    width = .5, alpha=0.5,
    #slab_color = 'black', 
    position = position_nudge(x = .05)
  ) +
  gghalves::geom_half_point(
    side = "1", 
    range_scale = .4, 
    alpha = .5, size = 3
  ) +
    stat_boxplot(geom="errorbar",width=0.1,size=1.5,position=position_nudge(x=0.05),color="black")+
    geom_boxplot(
    width = .15, color = "black",position=position_nudge(x=0.05),
    size = 1.5, outlier.shape = NA#,coef = 1.5  # 设置 coef 参数来调整箱线图的显示范围
  ) +
  stat_compare_means(aes(group=Type),method='wilcox.test',vjust=-20.5,hjust=1,size=8)+
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+scale_x_discrete(labels = c("Case" =paste0("Male",' (',table(boxplotData[,2])[names(table(boxplotData[,2]))=='Case'],') '),"Control" = paste0("Female",' (',table(boxplotData[,2])[names(table(boxplotData[,2]))=='Control'],') ')))+
  labs(x="",y="Phosphorylation level")+# title="Comparison of phosphorylation level between tumors with normals",
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5,
								  face = "bold"),
        axis.title.x = element_text(size = 19, 
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
  )+theme(legend.position = 'none')+coord_flip()
)
dev.off()
}
rownames(Agediffer)<-colnames(NormPhosphoGroupsfiltered)[2:ncol(NormPhosphoGroupsfiltered)]
Agediffer$wilcoxFDR<-p.adjust(Agediffer$wilcoxP,'BH')
write.table(Agediffer,file=paste0("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/Quantitative/Genderdiffer/",Cancertype,"_Genderdiffer.tsv"),sep="\t",quote=F,row.names=T,col.names=T)
}