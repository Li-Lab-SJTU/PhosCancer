CCRCC <- read.table("/dssg/home/Projects/PhosCancer/Analysis/Quantitative/TumorSizeCor/CCRCC_TumorSizeCor.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
HNSCC <- read.table("/dssg/home/Projects/PhosCancer/Analysis/Quantitative/TumorSizeCor/HNSCC_TumorSizeCor.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
PDAC <- read.table("/dssg/home/Projects/PhosCancer/Analysis/Quantitative/TumorSizeCor/PDAC_TumorSizeCor.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
UCEC <- read.table("/dssg/home/Projects/PhosCancer/Analysis/Quantitative/TumorSizeCor/UCEC_TumorSizeCor.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)

files<-c('CCRCC','HNSCC','PDAC','UCEC' )

Annotation <- read.table("/dssg/home/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv",comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
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
jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/TumorSizeCor/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_','_ALL.jpeg'),width = 700, height = 300)

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
    theme_clean()+
  scale_color_manual(values = c('brown2','black'))+		
  scale_fill_manual(values = c('#9FB3D2','#BE7F7D')) +
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5,
								  face = "bold"),
        axis.title = element_text(size = 19, 
                                    color = "black",
                                    face = "bold", 
									),
        legend.title = element_blank(),
 		#legend.position= c(0.9,0.12),
        legend.text = element_text(color="black", 
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 16,  
                                   color ='black', 
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0),
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