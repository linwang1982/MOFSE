*drug_all_info_664.pkl*: 从Drugbank中获取的每个药物的相关信息
*drug_info_664.csv*: 按频率矩阵中顺序的药物的SMILES及ID
*se994_info.xlsx*: 从ADReCS下载的副作用信息，按频率矩阵中的顺序
*drug_all_info_664.pkl*: 从DrugBank上获取的药物信息，包含其靶点信息

`get_SE_features.py`: -> *se_feature.mat*
    计算每个副作用的各个MedDRA术语级别的multi-hot编码 

`get_tani_sim.py`: 
    计算药物的 Tanimoto 相似度（用于可解释性分析部分）

`get_drugs_targets_info.py`: -> *domain_GO_664.mat*
    根据靶点信息，得到任一药物的所有Domain和GO(component/function/process)信息并编码为multi-hot

