library(stringr)
argv <- commandArgs(trailingOnly = T)
# Rscript ttest_Tumor_normal.R ./ssgsea_symbol.txt BRCA_Summary.txt

# zscore<-read.table('/home/arajo0707/TCGA_LINCS/TCGA_preprocess/ssgsea/TCGA_BRCA_ssGSEA_z-score.txt',sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)


# Normal 10,11,12,14
# Tumor 1,3,5,9
tumor_col <- c("-01","-03","-05","-09")
normal_col <- c("-10","-11","-12","-14")
cancer_name <- NULL
outlier <- NULL



filenames <- list.files(argv[1])



for(TCGA_zscore in filenames){
  print(TCGA_zscore)
  cancer_name <- strsplit(TCGA_zscore, "_")[[1]][2]
  
  if(endsWith(TCGA_zscore, "z-score.txt")){
    zscore<-read.table(paste0(argv[1],"/",TCGA_zscore),sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)
    normal <- zscore[, grep(paste(normal_col, collapse="|"), colnames(zscore)),drop=FALSE]
    tumor <- zscore[, grep(paste(tumor_col, collapse="|"), colnames(zscore)),drop=FALSE]
    
    TF_name <- c()
    pvalue <- c()
    mean <- NULL

    if(ncol(normal) < 2 || ncol(tumor) < 2 ){    
      outlier <- append(outlier,cancer_name)
      next
    }
                                          
    for(i in 1:nrow(zscore)){
      ttest <- t.test(tumor[i,],normal[i,],alternative = "two.sided", paired = FALSE, var.equal = TRUE)
      TF_name <- c(TF_name,zscore$Transcription_factor[i])
      pvalue <- c(pvalue,ttest$p.value)
      mean <- rbind.data.frame(mean,ttest$estimate)
    }
    
    result <- cbind.data.frame(TF_name,mean,pvalue, fdr_p <- p.adjust(pvalue,'fdr'), enrichment_score <- -log10(fdr_p))
    colnames(result) <- c('Transcription_factor','Tumor_mean','Normal_mean','pvalue', 'FDR', 'Enrichment_score')
    
    for(i in 1:nrow(result)){
      if(result$Tumor_mean[i] < result$Normal_mean[i]){
        result$Enrichment_score[i] <- result$Enrichment_score[i] * (-1)
      }
    }
    
    write.table(result,paste0(argv[2],"/" ,"TCGA_score_",cancer_name,".txt"),sep='\t',row.names=F,quote=F)
    
  }
  
}

write.table(outlier,paste0(argv[2],"/" ,"TCGA_score_outlier.txt"),sep='\t',quote = F)
