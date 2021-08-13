import pandas as pd
import numpy as np
import re
from cmapPy.pandasGEXpress.parse import parse

LINCS='/home/arajo0707/TCGA_LINCS/LINCS_preprocess'

# LINCS_G Gene 리스트 (duplicated, withdrawn등)
# LINCS_G Gene list (including duplicated, withdrawn)
# col_order 링스 데이터 콜론 정렬
# col_order for the LINCS_data column order

LINCS_G=pd.read_csv(f'{LINCS}/LINCS_Gene.txt',sep='\t')
col_order=pd.read_csv(f'{LINCS}/index_sorted_descending.txt',sep='\t')
subset_col_order=col_order.drop_duplicates(['new_drug_name'],keep='first')[['drugname','new_drug_name']].reset_index(drop=True)
LINCS_data=parse(f'{LINCS}LINCS_raw.gctx')
LINCS_data=LINCS_data.data_df


# dupelicated gene 리스트 생성
# creating duplicated gene list
dup_g=set(LINCS_G[LINCS_G['Symbol'].duplicated()]['Symbol']) 

# duplicated genes는 평균값을 계산한후 첫번째 값에 저장; 나머지값은 제거
# For duplicated genes, calculate the mean and replace the first duplicated gene; keep the first duplicated gene colunmn only.
for g in dup_g:
    dup_g=LINCS_G.index[LINCS_G['Symbol'] == g] #Find index of the gene
    LINCS_data.iloc[dup_g[0],] = LINCS_data.iloc[dup_g,].mean() #Replace duplicated genes the with mean 
    dup_g=dup_g[1:] # Select except first duplicated gene column 
    LINCS_data=LINCS_data[~LINCS_data.index.isin(dup_g)]
    LINCS_G=LINCS_G[~LINCS_G.index.isin(dup_g)]


# WithDrawn gene 제거
# Removing withDrawn
WithDr=LINCS_G.index[LINCS_G['Entrez']=='Withdrawn']

for g in WithDr:
    LINCS_data=LINCS_data.loc[~LINCS_data.index.isin([str(g)])]
    LINCS_G=LINCS_G[~LINCS_G.index.isin([str(g)])]

LINCS_G.index=LINCS_G['pr_gene_id']
LINCS_G=LINCS_G.loc[:,['Entrez','Symbol']]

# 같은 drugname셋에서 같은 new_drug_name끼리 평균내서 한개의 콜론으로 추출, 각 약물마다 반복
# In the same drugname set, average for same new_drug_name column; repeat for every drug in drugname set. 

for dg in set(col_order['drugname']):
    result_mat=pd.DataFrame()
    drug_set=col_order[col_order['drugname']==dg]
    dg = re.sub('[ ()/\'&]', '-', dg)

    for new_drug in set(drug_set['new_drug_name']):
        new_drug_subset=drug_set[drug_set['new_drug_name']==new_drug]
        Ldata_sub=LINCS_data.loc[:,new_drug_subset['Raw']]
        Ldata_sub=Ldata_sub.apply(np.mean,1)
        result_mat[new_drug['new_drug_name'].values[0]]=Ldata_sub

    LINCS_data_mat=LINCS_data_mat=pd.concat([LINCS_G.reset_index(drop=True),result_mat.reset_index(drop=True)],axis=1)
    column_names = subset_col_order.loc[subset_col_order["drugname"] == dg,"new_drug_name"].tolist()
    LINCS_data_mat = LINCS_data_mat.reindex(columns=column_names)

    # 각 약물에 대해서 파일로 저장하기
    # Save the drug subset file.
    LINCS_data_mat.reset_index().to_csv("".join(['/home/arajo0707/TCGA_LINCS/LINCS_preprocess/drug_subset/', dg, '.txt.gz']),compression='gzip',sep='\t',index=False)

