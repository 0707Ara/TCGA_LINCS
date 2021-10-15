library(GSVA)
library(stringr)
library(pheatmap)

argv<-commandArgs(trailingOnly = T)
# Function for row name
add_first_row <- function(df) {
  df1 <- rbind(c(colnames(df)), df)
  final_rowname <- c('Transcription_factor',rownames(df))
  rownames(df1) <- final_rowname
  return(df1)
}


# Load the gene set and restructure into list
geneset<-read.table(argv[1], header = TRUE, sep='\t',na.strings = "", fill = TRUE)
geneset<-read.table("/home/arajo0707/TCGA_LINCS/LINCS_preprocess/tft_benchmark_symbol_10TG.txt", header = TRUE, sep='\t',na.strings = "", fill = TRUE)
geneset_list<-list()
for(l in 1:nrow(geneset)){
  l1 <- geneset[l,4:ncol(geneset),drop=FALSE]
  l1 <- l1[!is.na(l1)]
  geneset_list[[geneset$TF[l]]] <- as.character(l1)
}
filenames <- list.files(argv[2])


# Load DMSO;the control group
DMSO<-read.table(paste0(argv[2],'/_DMSO_cancer.txt'),sep='\t',header=T,stringsAsFactors = F,check.names=FALSE)
rownames(DMSO)<-DMSO[,2]
DMSO<-DMSO[,3:ncol(DMSO)]

# Remove duplicated columns and leave values with 24 hour result
index <- read.table("/home/arajo0707/TCGA_LINCS/LINCS_preprocess/index_sorted_descending.txt",sep='\t',header=T,stringsAsFactors = F,na.strings = "", fill = TRUE)
selected_cell <- c('A375', 'A549', 'BT20', 'HCC515', 'HELA', 'HEPG2', 'HS578T', 'HT29', 'JURKAT', 'LNCAP', 'MCF7', 'MDAMB231',  'PC3', 'SKBR3', 'YAPC')
index_sorted <- index[!duplicated(index$new_drug_name),]
index_sorted <- index_sorted[index_sorted$cell_id %in% selected_cell & index_sorted$time == 24,]
index_sorted$drugname <- gsub('DMSO', '_DMSO', index_sorted$drugname)

DMSO_col <- index_sorted[index_sorted$drugname == "_DMSO" ,"new_drug_name"]
DMSO <- DMSO[,colnames(DMSO) %in% DMSO_col]
no_result <- NULL
final_matrix <- data.frame(matrix(NA, nrow = 12326, ncol = 0))

# Make matrix for GSVA function
for(drug_f in filenames){
  drug <- strsplit(drug_f,'_cancer.txt')[[1]][1]
  print(drug)
  if(drug!='_DMSO'){
    drug_mat <- read.table(paste0(argv[2],"/",drug_f),sep='\t',header = T,stringsAsFactors = F,check.names=FALSE)
    
    if(ncol(drug_mat)>2){
      drug_mat1 <-as.data.frame(drug_mat[,3:ncol(drug_mat)])
      rownames(drug_mat1) <- rownames(DMSO)
      new_col <- index_sorted[index_sorted$drugname == drug,"new_drug_name"]
      drug_mat1<-drug_mat1[,colnames(drug_mat1) %in% new_col]
      
      if( ncol(drug_mat1) > 1 ){
        final_matrix <- cbind.data.frame(final_matrix,drug_mat1)
      }else{
        no_result <- append(no_result,drug)
      }
    }
  }else{
    # For DMSO 
    final_matrix <- cbind.data.frame(DMSO,final_matrix)
    
  }
}

rownames(final_matrix) <- rownames(DMSO)
final_matrix2 <- add_first_row(final_matrix)
final_matrix <- as.matrix(final_matrix)
write.table(final_matrix2,paste0(argv[3],'/','LINCS_ssGSEA_input.txt'),sep='\t',quote = F)

# Apply gsva function
ssgsea_r <- gsva(final_matrix, gset.idx.list = geneset_list,method='ssgsea',
                 kcdf="Gaussian", abs.ranking=FALSE, min.sz=2,
                 max.sz=Inf, parallel.sz=10, mx.diff=TRUE,
                 ssgsea.norm=FALSE,
                 verbose=F) 

rownames(ssgsea_r) <- rownames(final_matrix)
ssgsea_r_final <- add_first_row(ssgsea_r)
write.table(ssgsea_r_final,paste0(argv[3],'/','LINCS_ssGSEA.txt'),sep='\t',quote = F, col.names=F)

# reformatting into z-score
ssgsea_r_zscore <- apply(ssgsea_r,2,scale)
colnames(ssgsea_r_zscore) <- colnames(ssgsea_r)
rownames(ssgsea_r_zscore) <- rownames(ssgsea_r)
rownames(ssgsea_r_zscore) <- rownames(final_matrix)
ssgsea_z_final <- add_first_row(ssgsea_r_zscore)

write.table(ssgsea_r_zscore,paste0(argv[3],'/','LINCS_ssGSEA_z-score.txt'),sep='\t',quote = F)

write.table(no_result,paste0(argv[3], "step3_GSVA_outlier.txt"),sep='\t',quote = F, col.names=F)

print("DONE")
