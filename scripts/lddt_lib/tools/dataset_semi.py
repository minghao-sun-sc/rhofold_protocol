"""Dataset for proteins."""
import shutil
import warnings
warnings.filterwarnings("ignore")

import os
import logging

import numpy as np
import random

from tqdm import tqdm
import torch
from torch.utils.data import Dataset
from rnafold.tools.dataset_config import DatasetConfig
from rnafold.tools.dataset import RNADataset, tns_keys
from rnafold.utils.comm_utils import exists

from rnafold.models.config import default_config

class RNADatasetSemi(Dataset):
    """The self-distillation dataset."""

    def __init__(self, config, ss_config, ratio = 0.75, is_train = True):
        """Constructor function."""

        super().__init__()

        self.ratio = ratio

        if isinstance(config, list):
            all_dataset = [RNADataset(config_, is_train=is_train) for config_ in config]
            self.dataset = torch.utils.data.ConcatDataset(all_dataset)
            crop_size = all_dataset[0].config.crop_size
        else:
            self.dataset = RNADataset(config=config, is_train=is_train)
            crop_size = self.dataset.config.crop_size

        self.semi_dataset = None
        if self.ratio > 0 and ss_config is not None:
            self.semi_dataset = RNADataset(config=ss_config, is_train=is_train)
            assert crop_size == self.semi_dataset.config.crop_size

        self.num = int(len(self.dataset) / (1 - ratio))
        # TODO: for fairseq training
        self.sizes = np.asarray([crop_size for i in range(self.num)])

    def __len__(self):
        """Get the number of elements in the dataset."""
        return self.num

    def __getitem__(self, idx):

        if self.ratio <= 0 or self.semi_dataset is None:
            return self.dataset.__getitem__(random.randint(0, len(self.dataset)-1))

        if random.uniform(0,1) < self.ratio:
            return self.semi_dataset.__getitem__(random.randint(0, len(self.semi_dataset)-1))
        else:
            return self.dataset.__getitem__(random.randint(0, len(self.dataset)-1))

def unitest_semi():
    msa_depth = 32                           # depth for MSA
    crop_size = 128
    pdb_version = 'v0.0'
    frame_version = 'v1.1'
    default_config.dataset.dist.bins = 40
    default_config.dataset.dist.max_dist = 39
    default_config.dataset.dist.min_dist = 0


    data_dir = '/public/home/taoshen/data/rna/RNA3D/Ash_Total_RNA_test/'
    id_fpath = f'{data_dir}/exp/valid_rnapuzzle_len20_drfoldv1.txt'

    dr_dpath = f'{data_dir}/rna_seq_rMSA_debug'
    dr_fea_dpath = f'{data_dir}/drfold_npzs'
    msa_dpath = f'{data_dir}/msa/'
    fas_dpath = f'{data_dir}/rna_seq/'          # directory path to SEQ files
    pdb_dpath = f'{data_dir}/rna_pdb_align/'    # directory path to KEYPDB files

    ds_train_config = DatasetConfig(id_fpath=id_fpath,
                                    fas_dpath=fas_dpath,
                                    msa_dpath=msa_dpath,
                                    pdb_dpath=pdb_dpath,
                                    crop_size=crop_size,
                                    msa_depth=msa_depth,
                                    pdb_version=pdb_version,
                                    frame_version=frame_version,
                                    dr_dpath=dr_dpath,
                                    dr_fea_dpath=dr_fea_dpath,
                                    is_dr_fea=True,
                                    is_dr_pred=True,
                                    dist_config=default_config.dataset.dist,
                                    )

    data_dir = '/public/home/taoshen/data/rna/bpRNA/'
    id_fpath = f'{data_dir}/v0.1.txt'
    dr_dpath = None
    dr_fea_dpath = f'{data_dir}/npzs'
    msa_dpath = f'{data_dir}/rMSA/'
    fas_dpath = f'{data_dir}/seq/'          # directory path to SEQ files
    pdb_dpath = f'{data_dir}/standard_model_1_unrefined_align/'    # directory path to KEYPDB files

    ds_ss_config = DatasetConfig(id_fpath=id_fpath,
                                    fas_dpath=fas_dpath,
                                    msa_dpath=msa_dpath,
                                    pdb_dpath=pdb_dpath,
                                    crop_size=crop_size,
                                    msa_depth=msa_depth,
                                    pdb_version=pdb_version,
                                    frame_version=frame_version,
                                    dr_dpath=dr_dpath,
                                    dr_fea_dpath=dr_fea_dpath,
                                    is_dr_fea=True,
                                    is_dr_pred=True,
                                    dist_config=default_config.dataset.dist,
                                    )

    ds_train_dataset = RNADatasetSemi(config = ds_train_config, ss_config = ds_ss_config, ratio=0.5)

    print(len(ds_train_dataset))
    for i in tqdm(range(len(ds_train_dataset))):
        data_dict = ds_train_dataset[i]
        for key in data_dict.keys():
            if key in set(tns_keys) and exists(data_dict[key]):
                print(i, key, data_dict[key].shape)
#
if __name__ == '__main__':
    unitest_semi()



