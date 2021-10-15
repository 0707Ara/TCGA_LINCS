import pandas as pd
import numpy as np

LINCS='/home/arajo0707/TCGA_LINCS/LINCS_preprocess'
# LINCS_Cell = Cell line 정보
# LINCS_G = Gene들의 최신화 (duplication, withdrawn등)
# GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz = LINCS id 정보

LINCS_data=pd.read_csv(f'{LINCS}/LINCS_newcol_mia2mean_withdrawn_head6164_1head.txt.gz', compression='gzip', sep='\t')
new_inst_id=pd.read_csv(f'{LINCS}/index_sorted_descending.txt',sep='\t')
for_colorder=pd.read_csv(f'{LINCS}/for_colorder.txt',sep='\t')

for dg in set(new_inst_id['drugname']):
    print(dg)
    s_inst=new_inst_id[new_inst_id['drugname']==dg]
    result_mat=pd.DataFrame()
    dg = re.sub('[ ()/\'&]', '-', dg)
    for cl in set(s_inst['cell_id']):
        cl_inst=s_inst[s_inst['cell_id']==cl]

        for new_drug in set(cl_inst['new_drug_name']):
            new_drug=cl_inst[cl_inst['new_drug_name']==new_drug]

            Ldata_sub=LINCS_data.loc[:,new_drug['Raw']]
            Ldata_sub=Ldata_sub.apply(np.mean,1)
            result_mat[new_drug['new_drug_name'].values[0]]=Ldata_sub
   

    LINCS_data_mat=result_mat.reset_index(drop=True)
    column_names = for_colorder.loc[for_colorder["drugname"] == dg,"new_drug_name"].tolist()
    LINCS_data_mat = LINCS_data_mat.reindex(columns=column_names) 
 # 각 약물에 대해서 파일로 저장하기
    LINCS_data_mat.reset_index().to_csv("".join(['/home/arajo0707/TCGA_LINCS/LINCS_preprocess/drug_subset/', dg, '_subset.txt.gz']),compression='gzip',sep='\t',index=False)
