library(GSVA)
library(openxlsx)

geneset<-read.xlsx('../../LINCS/20210415_New_GSVA/c3.tft.gtrd.v7.4.entrez_symbol.xlsx')

geneset_list<-list()
for(l in 1:nrow(geneset)){
  l1=geneset[l,8:ncol(geneset)]
  l1=l1[!is.na(l1)]
  geneset_list[[geneset$TF[l]]]=as.character(l1)
}


Tumor<-read.table('/home/jojo9103/LINCS/TCGA_XENA/BRCA/TCGA_RNA-Seq_BRCA_T.txt',sep='\t',header=T,stringsAsFactors = F,check.names = F)
Normal<-read.table('/home/jojo9103/LINCS/TCGA_XENA/BRCA/TCGA_RNA-Seq_BRCA_N.txt',sep='\t',header=T,stringsAsFactors = F,check.names = F)

Tumor1<-as.matrix(Tumor)
Normal1<-as.matrix(Normal)

Tumor1<-Tumor1[2:nrow(Tumor1),]
Normal1<-Normal1[2:nrow(Normal1),]


rownames(Tumor1)=Tumor1[,3]
rownames(Normal1)=Normal1[,3]

Tumor1<-Tumor1[,4:ncol(Tumor1)]
Normal1<-Normal1[,4:ncol(Normal1)]

result<-cbind(Normal1,Tumor1)
result<-as.matrix(result)

test=apply(result,2,as.numeric)
rownames(test)=rownames(result)
rownames(test)=trimws(rownames(test))

ssgsea_r<-gsva(test,geneset_list,method='ssgsea',kcdf="Gaussian",abs.ranking=FALSE,
                     min.sz=2,max.sz=Inf,parallel.sz=10,mx.diff=TRUE,
                     ssgsea.norm=TRUE,
                     verbose=TRUE)

write.table(ssgsea_r,'/home/jojo9103/LINCS/TCGA_XENA/BRCA/GSVA_result/ssgsea_symbol.txt',sep='\t',quote = F)

