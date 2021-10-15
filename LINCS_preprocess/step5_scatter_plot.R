library(ggplot2)
library(reshape2)
lincs_score <- read.table('/home/arajo0707/TCGA_LINCS/LINCS_preprocess/LINCS_score/Enrichment_score_total.txt',sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)
tcga_score <- read.table('/home/arajo0707/TCGA_LINCS/TCGA_preprocess/TCGA_score/Enrichment_score_TCGA_total.txt',sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)
rowName <- lincs_score$Transcription_factor

lincs_score <- lincs_score[,2:ncol(lincs_score), drop = FALSE]
tcga_score <- tcga_score[,2:ncol(tcga_score), drop = FALSE]
drugname_set <- colnames(lincs_score)
cancername_set <- colnames(tcga_score)

for(i in cancername_set){
   # i <- "BRCA"
  result_matrix <- cbind.data.frame(tcga_score[,i,drop=FALSE],lincs_score)
  melted_result <- melt(result_matrix, id.vars=i)
  print(i)
  for(j in drugname_set) {
    print(j)
    # j <- "10-DEBC"
    p <-ggplot(melted_result)+geom_point(aes(x=value,y=eval(parse(text = i)),col=variable),col = "#999999",size=.5)+
        theme_bw()+xlab('LINCS enrichment score')+ylab('TCGA enrichment score')+
        theme(legend.position = "None",axis.title = element_text(size=15),axis.text = element_text(size=10))+
        geom_vline(xintercept = 0,linetype='dashed',col='brown')+geom_hline(yintercept = 0,linetype='dashed',col='brown')+
        ggtitle(paste0(i,"_",j))+geom_point(data = subset(melted_result,variable == j),aes(x=value,y=eval(parse(text = i)),col="#E69F00"),size=0.5)
                                            
    ggsave(filename = paste0('/home/arajo0707/TCGA_LINCS/LINCS_preprocess/LINCS_score/scatterplot/',i,'/',i,'_',j,'.tiff'),plot = p,width = 200,height = 100,units = 'mm')
  }
  
}

