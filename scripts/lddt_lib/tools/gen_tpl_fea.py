import os
import random
import pickle
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from rnafold.tools.alphabet import _get_msa_feature

align_socre_dict = {
    ' ':0,
    '.':1,
    ':':2,
}

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

########################################################################################################################
def parse_tpl_file(pkl_path, fas_path,  tpl_topk, is_train):
    """Parse the TPL file (structural templates)."""

    idx_mask = -1

    f = open(fas_path, "rb")
    aa_seq = f.readlines()[1].strip()
    f.close()

    return get_tpl_feature(
            pkl_path, tpl_topk=tpl_topk, is_train=is_train, mask_index=idx_mask, aa_seq = aa_seq)

########################################################################################################################
def get_tpl_feature(tpl_path,
                    tpl_topk,
                    is_train,
                    mask_index,
                    aa_seq):

    d1d, d2d = 12, 5

    seq_len = len(aa_seq)
    t1ds = torch.ones([tpl_topk, seq_len, d1d]) * mask_index
    t2ds = torch.ones([tpl_topk, seq_len, seq_len, d2d]) * mask_index
    align_seq_tokens = None

    if tpl_path is None:
        return t1ds, t2ds, align_seq_tokens

    try:
        with open(tpl_path, 'rb') as f:
            feature_dict = pickle.load(f)

        return get_tpl_fea(
                       aa_seq,
                       feature_dict,
                       tpl_num = tpl_topk,
                       mask_index = mask_index,
                       is_train = is_train,
                       )
    except:
        print(f'sth is wrong in {tpl_path}')
        return t1ds, t2ds, align_seq_tokens

def _get_align_score_token(align_score):
    token = np.asarray([align_socre_dict[s] for s in align_score])
    return torch.LongTensor(token).view(1, -1)

# template feature v1
def get_tpl_fea(aa_seq, feature_dict, tpl_num = 4, mask_index = -1, is_train = False):
    #  Return:
    #   - t1d: 1D template info (B, T, L, 2)
    #   - t2d: 2D template info (B, T, L, L, 10)

    tpl_num_pkl = len(feature_dict)
    seq_len = len(aa_seq)

    id_list = list(feature_dict.keys())
    id_list.sort()

    t1ds = np.ones([tpl_num, seq_len, 12]) * mask_index
    t2ds = np.ones([tpl_num, seq_len, seq_len, 5]) * mask_index

    if tpl_num_pkl == 0:
        t1ds = torch.FloatTensor(t1ds)
        t2ds = torch.FloatTensor(t2ds)
        return t1ds, t2ds, None

    #for data augmentation
    if is_train and random.randint(0, 1) == 0:
        t1ds = torch.FloatTensor(t1ds)
        t2ds = torch.FloatTensor(t2ds)
        return t1ds, t2ds, None

    align_seq_tokens = []
    for tpl in range(tpl_num):
        tpl_index = tpl

        # for data augmentation
        if is_train and random.randint(0, 1) == 0:
            tpl = random.randint(0, tpl_num_pkl + tpl_num - 1)

        if tpl + 1 > tpl_num_pkl:
            break

        id = id_list[tpl]

        align_seq = feature_dict[id]['align_seqs'][1]
        align_seq_token = _get_msa_feature([[('xxx', align_seq)]], is_rm=False)
        align_seq_tokens.append(align_seq_token)

        align_score = feature_dict[id]['align_seqs'][2]
        align_score_token = _get_align_score_token(align_score)

        dist = feature_dict[id]['dist']
        align_seq_fea = F.one_hot(align_seq_token, 9)
        align_score_fea = F.one_hot(align_score_token, 3)
        t1d = np.concatenate([align_seq_fea, align_score_fea], axis = -1)

        t1ds[tpl_index, :, :] = t1d
        t2ds[tpl_index, :, :, :] = dist.transpose([1, 2, 0])

    align_seq_tokens = torch.cat(align_seq_tokens, dim=0) if len(align_seq_tokens) > 0 else None

    t1ds = torch.FloatTensor(t1ds)
    t2ds = torch.FloatTensor(t2ds)

    return t1ds, t2ds, align_seq_tokens

if __name__ == '__main__':
    root = '/public/home/taoshen/data/rna/RNA3D/Ash_Total_RNA_test'

    tpl_dir = f'{root}/rna_tpl'
    seq_dir = f'{root}/rna_seq'

    for tgt in os.listdir(tpl_dir):
        tgt = tgt.replace('.pkl','')
        fas_path = f'{seq_dir}/{tgt}.seq'
        pkl_path = f'{tpl_dir}/{tgt}.pkl'
        n_resds = 16
        is_train = False
        tpl_topk = 4

        try:
            parse_tpl_file(pkl_path, fas_path, tpl_topk, is_train)
        except:
            print('xxx')
        exit()





