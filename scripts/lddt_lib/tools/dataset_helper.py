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
from scipy.spatial.distance import cdist
from rnafold.tools.rna_utils import parse_fas_file
from rnafold.tools.dataset_config import DatasetConfig
from rnafold.tools.pdb_parser import PdbParser
from rnafold.utils.comm_utils import exists
from rnafold.components.rna_esm.data import get_rna_fm_token
from rnafold.tools.alphabet import get_msa_feature
from rnafold.tools.gen_tpl_fea import parse_tpl_file
from rnafold.tools.cal_da_lbl import *


tns_keys = ['token',
            'rna_fm_token',
            't1ds',
            't2ds',
            'cord_gt',
            'cord_gt_mask',
            'dist_bins',
            'cmap',
            'cmap_msk',
            'dr_msa',
            'dr_seq',
            'dr_hmm',
            'dr_ss',
            'dr_seq_fea',
            'dr_pair_fea',
            'ssc_i',
            'ssc_j',
            'ssc_dist',
            ]

def get_split_ids(complex_id):
    if '_' not in complex_id:
        return [complex_id]
    item = complex_id.split('_')
    return [f'{item[0]}_{i}' for i in item[1:]]

from rnafold.tools.alphabet import Alphabet
from rnafold.components.rna_esm.data import Alphabet as RNAFMAlphabet

alphabet = Alphabet.from_architecture('RNA MSA Transformer')
rnafm_alphabet = RNAFMAlphabet.from_architecture('ESM-1b', theme="rna")

def sample_tokens(msa_tokens_full, msa_depth, is_train, mask_prob=0.15, mode = 'random',
                  rna_fm_tokens = None):
    """Sample the MSA w/ optional random masks applied."""

    # initialization
    device = msa_tokens_full.device
    # sample a fixed number of sequences from the original MSA data
    n_seqs = msa_tokens_full.shape[1]

    if n_seqs <= msa_depth:
        msa_tokens_true = msa_tokens_full

    elif mode == 'topN':
        # use top-K' sequences during evaluation
        msa_tokens_true = msa_tokens_full[:, :msa_depth]

    elif mode == 'random':
        # randomly select K' sequences during training
        idxs_seq = [0]
        idxs_seq = (idxs_seq + list(random.sample(range(1, n_seqs), msa_depth - 1))) \
            if n_seqs > 1 else idxs_seq
        msa_tokens_true = torch.stack([msa_tokens_full[:, x] for x in idxs_seq], dim=1)
    else:
        raise NotImplementedError

    # apply random masks on MSA tokens
    msa_tokens_pert = msa_tokens_true.detach().clone()
    msa_masks = (torch.rand(msa_tokens_pert.shape, device=device) < mask_prob).to(torch.int8)

    if is_train:
        msa_tokens_pert[msa_masks == 1] = alphabet.mask_idx
        if exists(rna_fm_tokens):
            rna_fm_tokens[msa_masks[:,0,:] == 1] = rnafm_alphabet.mask_idx

    return msa_tokens_true, msa_tokens_pert, msa_masks, rna_fm_tokens

