library(GSVA)
library(stringr)

argv <- commandArgs(trailingOnly = T)

# Define functions for the program
add_first_row <- function(df) {
  df1 <- rbind(c(colnames(df)), df)
  final_rowname <- c('Transcription_factor',rownames(df))
  rownames(df1) <- final_rowname
  return(df1)
}

clean_data <- function(df){
  normal_mat <- df[2:nrow(df),4:ncol(df),drop=FALSE]
  return(normal_mat)
}

# Loading 
geneset <- read.table(argv[1], header = TRUE, sep='\t',na.strings = "", fill = TRUE)
geneset_list <- list()
for(l in 1:nrow(geneset)){
  l1 <- geneset[l,4:ncol(geneset)]
  l1 <- l1[!is.na(l1)]
  geneset_list[[geneset$TF[l]]] <- as.character(l1)
}

filenames <- list.files(argv[2])



for(cancer_name in filenames){
  print(cancer_name) 
  cancer_list <- list.files(paste0(argv[2],"/",cancer_name))
  normal_mat <- NULL
  tumor_mat <- NULL
  final_matrix <- NULL
  rowname <- NULL
  
  for(cancer_sample in cancer_list){
    
    #
    if(endsWith(cancer_sample, "_N.txt")){
      normal_mat <- read.table(paste0(argv[2],"/",cancer_name,"/",cancer_sample),sep='\t',header = T,stringsAsFactors = F,check.names=FALSE)
      rowname <- normal_mat$Symbol[2:nrow(normal_mat)]
      normal_mat <- clean_data(normal_mat)
      
    }else if(endsWith(cancer_sample, "_T.txt")) {
      tumor_mat <- read.table(paste0(argv[2],"/",cancer_name,"/",cancer_sample),sep='\t',header = T,stringsAsFactors = F,check.names=FALSE)
      rowname <- tumor_mat$Symbol[2:nrow(tumor_mat)]
      tumor_mat <- clean_data(tumor_mat)
    }
  }
  
  # Formatting final matrix for GSVA function
  if(is.null(ncol(tumor_mat))){
    final_matrix <- normal_mat
    
  }else if(is.null(ncol(normal_mat))){
    final_matrix <- tumor_mat
    
  }else{
    final_matrix <- cbind.data.frame(normal_mat,tumor_mat)
    
  }
  final_matrix <- as.matrix(final_matrix)
  final_matrix <- apply(final_matrix,2,as.numeric)
  row.names(final_matrix) <- rowname
  
  # Save GSVA input matrix
  final_matrix1 <- add_first_row(final_matrix)
  write.table(final_matrix1,paste0(argv[3],'/','TCGA_',cancer_name,'_full.txt'),sep='\t',quote = F, col.names=F)
  
  print("GSVA function is running")
  suppressWarnings({ 
    ssgsea_r <- gsva(final_matrix, gset.idx.list = geneset_list, method='ssgsea',
                     kcdf="Gaussian", abs.ranking=FALSE, min.sz=2,
                     max.sz=Inf, parallel.sz=10, mx.diff=TRUE,
                     ssgsea.norm=FALSE,
                     verbose=F) 
    })

  # Save GSVA output matrix
  ssgsea_r_final <- add_first_row(ssgsea_r)
  write.table(ssgsea_r_final,paste0(argv[3],'/','TCGA_',cancer_name,'_ssGSEA.txt'),sep='\t',quote = F, col.names=F)
  
  # Reformatting the ssgsea matrix into z-score
  ssgsea_r_zscore <- apply(ssgsea_r,2,scale)
  ssgsea_r_zscore <- as.data.frame(ssgsea_r_zscore)
  colnames(ssgsea_r_zscore) <- colnames(ssgsea_r)
  rownames(ssgsea_r_zscore) <- rownames(ssgsea_r)
  ssgsea_z_final <- add_first_row(ssgsea_r_zscore)
  write.table(ssgsea_z_final,paste0(argv[3],'/','TCGA_',cancer_name,'_ssGSEA_z-score.txt'),sep='\t',quote = F,col.names=F)
  
  
}
