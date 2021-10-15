import pandas as pd
import numpy as np
import re

LINCS='/home/arajo0707/TCGA_LINCS/LINCS_preprocess'
# selected cell line info
#entrez=pd.read_csv(f'{LINCS}/LINCS_Gene_noseq_entrez_symbol.txt',sep='\t')

col_index=pd.read_csv(f'{LINCS}/index_sorted_descending.txt',sep='\t')

selected_gene = ['A375', 'A549', 'BT20', 'HCC515', 'HELA', 'HEPG2', 'HS578T', 'HT29', 'JURKAT', 'LNCAP', 'MCF7', 'MDAMB231',  'PC3', 'SKBR3', 'YAPC']
col_order=col_index.drop_duplicates(['new_drug_name'],keep='first').reset_index(drop=True)
col_order=col_order.loc[col_order['cell_id'].isin(selected_gene)]

for dg in set(col_index['drugname']):
    result_mat=pd.DataFrame()
    dg = re.sub('[ ()/\'&]', '-', dg)
    print(dg)
    subset=pd.read_csv(f'{LINCS}/drug_subset_merge/{dg}_subset.txt.gz',sep='\t',compression='gzip',dtype='str')
    #subset=subset.iloc[:,2:]
    subset=subset.loc[:,col_order.loc[col_order["drugname"] == dg,"new_drug_name"].tolist()]
    #subset.to_csv(f'{LINCS}/drug_cancer_gz/{dg}_cancer.txt.gz',sep='\t',compression='gzip' ,index=False) 
    if not subset.empty:
        subset.to_csv(f'{LINCS}/drug_cancer_txt/{dg}_cancer.txt',sep='\t',index=False) 
    print(dg)
