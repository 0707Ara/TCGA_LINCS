library(stringr)

zscore <- read.table('/home/arajo0707/TCGA_LINCS/LINCS_preprocess/ssgsea/LINCS_ssGSEA_z-score.txt',sep='\t',header=T,stringsAsFactors=F,quote='',check.names=FALSE)
TF_name <- c(zscore$Transcription_factor)
zscore <- zscore[,2:ncol(zscore)]

print("zscore readed")
index <- read.table("/home/arajo0707/TCGA_LINCS/LINCS_preprocess/index_sorted_descending.txt",sep='\t',header=T,stringsAsFactors = F,na.strings = "", fill = TRUE)
selected_cell <- c('A375', 'A549', 'BT20', 'HCC515', 'HELA', 'HEPG2', 'HS578T', 'HT29', 'JURKAT', 'LNCAP', 'MCF7', 'MDAMB231',  'PC3', 'SKBR3', 'YAPC')
index_sorted <- index[!duplicated(index$new_drug_name),]
index_sorted <- index_sorted[index_sorted$cell_id %in% selected_cell & index_sorted$time == 24,]
index_sorted$drugname <- gsub('DMSO', '_DMSO', index_sorted$drugname)
row.names(index_sorted) <- NULL



# Extract DMSO(Ctl) samples
ctl <- index_sorted[index_sorted$drugname == "_DMSO" ,]
ctl_zscore <- zscore[,colnames(zscore) %in% ctl$new_drug_name]
trt <- index_sorted[!(index_sorted$drugname == "_DMSO") ,]

outlier <- NULL

trt_drugnames <- unique(trt$drugname)
for (drugname in trt_drugnames) {
  #drugname <- "10-DEBC"
  print(drugname)
  
  # Extract treatment group sampele for each drug
  trt_drug <- trt[trt$drugname == drugname,]
  trt_drug <- trt_drug[trt_drug$dose == max(as.numeric(trt_drug$dose)),]
  trt_zscore <- zscore[,colnames(zscore) %in% trt_drug$new_drug_name, drop=FALSE]
  
  pvalue <- c()
  mean <- NULL
  
  if(ncol(trt_zscore) < 2 ){    
    outlier <- append(outlier,drugname)
    next
  }
  
  # Do the student ttest 
  for(i in 1:nrow(ctl_zscore)){
    ttest <- t.test(trt_zscore[i,],ctl_zscore[i,],alternative = "two.sided", paired = FALSE, var.equal = TRUE)
    pvalue <- c(pvalue,ttest$p.value)
    mean <- rbind.data.frame(mean,ttest$estimate)
  }
  
  result <- cbind.data.frame(TF_name,mean,pvalue, fdr_p <- p.adjust(pvalue,'fdr'), enrichment_score <- -log10(fdr_p))
  colnames(result) <- c('Transcription_factor','Trt_mean','Ctl_mean','pvalue', 'FDR', 'Enrichment_score')
  
  for(i in 1:nrow(result)){
    if(result$Trt_mean[i] < result$Ctl_mean[i]){
      result$Enrichment_score[i] <- result$Enrichment_score[i] * (-1)
    }
  }
  
  write.table(result,paste0("./LINCS_score/","LINCS_score_",drugname,".txt"),sep='\t',row.names=F,quote=F)
}
write.table(outlier,paste0("./LINCS_score/" ,"LINCS_score_outlier.txt"),sep='\t',quote = F)
