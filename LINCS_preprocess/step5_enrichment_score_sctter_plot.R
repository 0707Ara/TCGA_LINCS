################################
# Making LINCS enrichment file #
################################
filenames <- list.files("/home/arajo0707/TCGA_LINCS/LINCS_preprocess/LINCS_score")

final_matrix <- NULL
rownimi <- NULL

for(LINCS_score in filenames){
  #LINCS_score <- "LINCS_score_leflunomide.txt"
  drug_name <- gsub("LINCS_score_", "", LINCS_score)
  drug_name <- gsub(".txt", "", drug_name)
  print(drug_name)
  lincs_score<-read.table(paste0('/home/arajo0707/TCGA_LINCS/LINCS_preprocess/LINCS_score/',LINCS_score),sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)
  # lincs_score<-read.table("/home/arajo0707/TCGA_LINCS/LINCS_preprocess/LINCS_score/LINCS_score_leflunomide.txt",sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)
  if(length(final_matrix) == 0){
    rownimi <- c(rownimi,drug_name)
    final_matrix <- lincs_score[,6, drop=FALSE]
    
  }
  else{
    final_matrix <- cbind.data.frame(final_matrix,lincs_score[,6, drop=FALSE])
    rownimi <- c(rownimi,drug_name)
    colnames(final_matrix) <- rownimi
  }
}
write.table(final_matrix,paste0(argv[2],"/" ,"TCGA_score_",cancer_name,".txt"),sep='\t',row.names=F,quote=F)
final_matrix <- cbind.data.frame(lincs_score[,1, drop=FALSE],final_matrix)
write.table(final_matrix,"/home/arajo0707/TCGA_LINCS/LINCS_preprocess/LINCS_score/Enrichment_score_LINCS_total.txt",sep='\t',row.names=F,quote=F)



################################
# Making TCGA enrichment file #
################################
filenames <- list.files("/home/arajo0707/TCGA_LINCS/TCGA_preprocess/TCGA_score")

final_matrix <- NULL
rownimi <- NULL

for(TCGA_score in filenames){
  #LINCS_score <- "LINCS_score_leflunomide.txt"
  if(TCGA_score == "TCGA_score_outlier.txt"){
    next
  }
  drug_name <- gsub("TCGA_score_", "", TCGA_score)
  drug_name <- gsub(".txt", "", drug_name)
  print(drug_name)
  tcga_score<-read.table(paste0('/home/arajo0707/TCGA_LINCS/TCGA_preprocess/TCGA_score/',TCGA_score),sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)
  # lincs_score<-read.table("/home/arajo0707/TCGA_LINCS/LINCS_preprocess/LINCS_score/LINCS_score_leflunomide.txt",sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)
  if(length(final_matrix) == 0){
    rownimi <- c(rownimi,drug_name)
    final_matrix <- tcga_score[,6, drop=FALSE]
    
  }
  else{
    final_matrix <- cbind.data.frame(final_matrix,tcga_score[,6, drop=FALSE])
    rownimi <- c(rownimi,drug_name)
    colnames(final_matrix) <- rownimi
  }
}


final_matrix <- cbind.data.frame(tcga_score[,1, drop=FALSE],final_matrix)
write.table(final_matrix,"/home/arajo0707/TCGA_LINCS/TCGA_preprocess/TCGA_score/Enrichment_score_TCGA_total.txt",sep='\t',row.names=F,quote=F)





###############################################
# Plotting scatter plot with enrichment files #
###############################################

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


###############################################
#             Calculating L-score             #
###############################################

lincs_score <- read.table('/home/arajo0707/TCGA_LINCS/LINCS_preprocess/LINCS_score/Enrichment_score_total.txt',sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)
tcga_score <- read.table('/home/arajo0707/TCGA_LINCS/TCGA_preprocess/TCGA_score/Enrichment_score_TCGA_total.txt',sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)
lincs_score <- lincs_score[,2:ncol(lincs_score), drop = FALSE]
tcga_score <- tcga_score[,2:ncol(tcga_score), drop = FALSE]

sum_total <- NULL
sum_total_cancer <- NULL
final_result <- NULL
for(i in 1:ncol(tcga_score)){
  print(i)
  for(j in 1:ncol(lincs_score)){
    print(j)
    mul <- lincs_score[,j] * tcga_score[,i]
    sum_mul <- sum(mul, na.rm = FALSE)
    # print(sum_mul)
    sum_total <- append(sum_total,sum_mul)
    
  }
  sum_total2 <- sum_total
  final_result<-rbind.data.frame(final_result,sum_total)
  sum_total <- NULL
}
rownames(final_result) <- colnames(tcga_score)
colnames(final_result) <- colnames(lincs_score)
write.table(final_result,"/home/arajo0707/TCGA_LINCS/LINCS_preprocess/L-score_final_result_file.txt",sep='\t',row.names=F,quote=F)

