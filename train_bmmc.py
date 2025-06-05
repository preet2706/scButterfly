import anndata
import scanpy as sc 
import pandas as pd 
import numpy as np
import scipy


from scButterfly.butterfly import Butterfly
butterfly = Butterfly()

RNA_data = anndata.read_h5ad("/workspace/scButterfly/data/bmmc/ds_RNA_bmmc.h5ad")
ATAC_data = anndata.read_h5ad("/workspace/scButterfly/data/bmmc/ds_ATAC_bmmc.h5ad")

from scButterfly.split_datasets import *
id_list = five_fold_split_dataset(RNA_data, ATAC_data)
train_id, validation_id, test_id = id_list[3]

butterfly.load_data(RNA_data, ATAC_data, train_id, test_id, validation_id)

butterfly.data_preprocessing()

chrom_list = []
last_one = ''
for i in range(len(butterfly.ATAC_data_p.var_names)):
    temp = butterfly.ATAC_data_p.var_names[i][0:5]
    if temp[0 : 3] == 'chr':
        if not temp == last_one:
            chrom_list.append(1)
            last_one = temp
        else:
            chrom_list[-1] += 1
    else:
        chrom_list[-1] += 1

butterfly.augmentation(aug_type=None)

butterfly.construct_model(chrom_list=chrom_list)

butterfly.train_model()

A2R_predict, R2A_predict = butterfly.test_model()
A2R_predict.write('/workspace/scButterfly/data/bmmc/pred_RNA_bmmc.h5ad')
R2A_predict.write('/workspace/scButterfly/data/bmmc/pred_ATAC_bmmc.h5ad')