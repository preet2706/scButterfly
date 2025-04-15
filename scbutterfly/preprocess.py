#!/usr/bin/env python
# coding: utf-8
  
import anndata
import scanpy as sc 
import pandas as pd 
import numpy as np
import scipy

from scButterfly.butterfly import Butterfly
butterfly = Butterfly()

RNA_data = anndata.read_h5ad("data/lymphoma_RNA.h5ad")
ATAC_data = anndata.read_h5ad("data/lymphoma_ATAC.h5ad")

from scButterfly.split_datasets import *
id_list = five_fold_split_dataset(RNA_data, ATAC_data)
train_id, validation_id, test_id = id_list[4]

butterfly.load_data(RNA_data, ATAC_data, train_id, test_id, validation_id)

butterfly.data_preprocessing(file_path='scbutterfly', save_data=True)