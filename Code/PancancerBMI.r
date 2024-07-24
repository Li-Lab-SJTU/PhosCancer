CCRCC <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/Quantitative/BMICor/CCRCC_BMICor.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
GBM <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/Quantitative/BMICor/GBM_BMICor.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
LUAD <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/Quantitative/BMICor/LUAD_BMICor.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
PDAC <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/Quantitative/BMICor/PDAC_BMICor.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
UCEC <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/Quantitative/BMICor/UCEC_BMICor.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)

files<-c('CCRCC','GBM','LUAD','PDAC','UCEC' )

Annotation <- read.table("/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
sites<-Annotation[,1]

library(reshape)
if(length(sites)>0){
for(site in 1:length(sites)){
whichsite=sites[site]
Total<-c()
for(i in 1:length(files)){
singles <-get(files[i])
if(length(which(rownames(singles)==whichsite))!=0){
SiteData<-singles[which(rownames(singles)==whichsite),]
Total<-rbind(Total,cbind(SiteData,CancerType=gsub('.tsv','',files[i])))
	}
}
if(!is.null(Total)){
jpeg(filename=paste0('/dssg/home/acct-clslj/clslj-2/Dongqun/Projects/PhosCancer/Analysis/PNGs/BMICor/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_','_ALL.jpeg'),width = 700, height = 300)

Total$Ecolors<-as.character(ifelse(Total$SpearmanE<=0,'Neg','Pos'))
Total$hjust<-ifelse(Total$SpearmanE<=0,0,1)
library(ggthemes)
library(ggplot2)
ggplot(Total,aes(x = reorder(CancerType,SpearmanE),y=SpearmanE,
           fill = Ecolors,
           label = paste0('p= ',round(SpearmanP,3))
		   ))+
    geom_col(show.legend = T,width=0.8) +
    coord_flip() +
    geom_text(nudge_x = 0.1,vjust=1,hjust=Total$hjust,size=6,face = "bold") +
    xlab('') +ylab('Correlation') +
    #ggtitle('Expression level of significant changed genes')+
    theme_clean()+
  scale_color_manual(values = c('brown2','black'))+		
  scale_fill_manual(values = c('#9FB3D2','#BE7F7D')) +
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
        legend.title = element_blank(),
 		#legend.position= c(0.9,0.12),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 16,  # 修改X轴上字体大小，
                                   color ='black', # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 16,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  ) 
dev.off() 
 }
  }
}