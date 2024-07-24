setwd('/dssg/home/Projects/PhosCancer/Analysis/Quantitative/Surdiffer/')
Annotation <- read.table('/dssg/home/Projects/PhosCancer/Analysis/TotalPhosphoprotein_annotation.tsv',comment.char="",quote='',sep='\t',header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)

files<-list.files(,pattern='_Surdiffer')
Total<-c()
for(i in 1:length(files)){
singlefile <- read.table(files[i],comment.char="",quote='',sep='\t',row.names=1,header = TRUE,fill=T,check.names=F,stringsAsFactors = FALSE)
singlefile$OSlower<-as.numeric(gsub('-.*','',singlefile[,2]))
singlefile$OSupper<-as.numeric(gsub('.*-','',singlefile[,2]))
singlefile$DFSlower<-as.numeric(gsub('-.*','',singlefile[,7]))
singlefile$DFSupper<-as.numeric(gsub('.*-','',singlefile[,7]))

singlefile<-cbind(Sites=rownames(singlefile),singlefile)
singlefile<-cbind(Cancer=gsub('_Surdiffer.tsv','',files[i]),singlefile)
Total<-rbind(Total,singlefile)
}
Total<-cbind(Total,Annotation[match(Total[,2],Annotation[,1]),c(4,3)])

sites<-unique(Total[,2])
for(site in  1:length(sites)){ 
whichsite=sites[site]

Data<-rbind(data.frame(Cancer=Total[Total[,2]==whichsite,'Cancer'],Mean=as.numeric(Total[Total[,2]==whichsite,'OSHazardRatio']),lower=Total[Total[,2]==whichsite,'OSlower'],upper=Total[Total[,2]==whichsite,'OSupper'],pvalue=paste('p =',as.numeric(Total[Total[,2]==whichsite,'OSCoxP'])),Type='OS'),data.frame(Cancer=Total[Total[,2]==whichsite,'Cancer'],Mean=as.numeric(Total[Total[,2]==whichsite,'DFSHazardRatio']),lower=Total[Total[,2]==whichsite,'DFSlower'],upper=Total[Total[,2]==whichsite,'DFSupper'],pvalue=paste('p =',as.numeric(Total[Total[,2]==whichsite,'DFSCoxP'])),Type='DFS'))
Data$Type<-factor(Data$Type,levels=c('OS','DFS'))
Data<-na.omit(Data)
library(ggplot2)

load('/dssg/home/Projects/PhosCancer/Analysis/PNGs/colorss.RData')
selectedcolorss<-colorss[match(gsub('\\(.*','',Data$Cancer),names(colorss))]
names(selectedcolorss)<-Data$Cancer

errorheight=1.2/length(unique(Data[,1]))
pointsize=20/length(unique(Data[,1]))
pngheight=150+30*length(unique(Data[,1]))
dotCOLS = c("#a6d8f0","#f9b282")
barCOLS = c("#008fd5","#de6b35")
jpeg(filename=paste0('/dssg/home/Projects/PhosCancer/Analysis/PNGs/Surdiffer/',Annotation[Annotation[,1]==whichsite,4],'_',Annotation[Annotation[,1]==whichsite,3],'_ALL.jpeg'),width = 1000, height = pngheight)
print(ggplot(Data, aes(x=Cancer, y=Mean, ymin=lower, ymax=upper,col=Type,fill=Type)) + 
  geom_hline(yintercept=1, lty=2) +
  geom_linerange(size=6,position=position_dodge(width = 0.5)) +
  ylab('Hazard ratios with 95% CI') + xlab("")+
  geom_point(size=3.6, shape=21, colour="white", stroke = 1,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  coord_flip() +  facet_grid(~ Type, scales = "free_x") +
  theme_minimal()+
  theme(strip.text.x = element_text(size = 14,face = "bold", colour = "black")) + 
  theme(legend.position='none',axis.title.x = element_text(size = 14, 
                                    color = "black",
                                    face = "bold"),
        axis.text.x = element_text(size = 12,  
                                   color = "black", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0), 
        axis.text.y = element_text(size = 14,  
                                   color = "black",
                                   face = "bold", 
                                   angle = 0)) 
								   )
dev.off()
}

