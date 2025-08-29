"""Build additional da lbl files."""

try:
    import matplotlib.pyplot as plt
except:
    print('matplotlib is not installed')

import os
import logging
import random

import torch
from scipy.spatial.distance import cdist
from rnafold.tools.pdb_parser import PdbParser
import numpy as np

from p_tqdm import p_map
from rnafold.utils import zfold_init
from tqdm import tqdm

def read_tgts(file):
    ''' '''
    f = open(file)
    tgts = f.readlines()
    f.close()
    tgts = [t.strip() for t in tgts]
    return tgts

def write_tgts(file, tgts):
    ''' '''
    f = open(file, 'w')
    for t in tgts:
        f.write(f'{t}\n')
    f.close()

def build_dist_bin_labels(dist, mask, bins = 37, max_dist = 20, min_dist = 2):
    """Get distanace ground truth labels."""

    # trrosetta 37-20-2
    # deeprnafold 40-39-0

    label_dis = np.zeros(dist.shape)
    bins = bins - 1
    assert bins % (max_dist - min_dist) == 0
    dist_inter = (max_dist - min_dist) / bins
    for i in range(bins):
        s = min_dist + dist_inter * i
        e = min_dist + dist_inter * (i + 1)
        label_dis[(dist > s) & (dist < e)] = i + 1
    label_dis[mask == 0] = -1

    return label_dis

def calc_dihedral_angle(cord_1, cord_2, cord_3, cord_4):
    """Calculate the dihedral angle defined by 4 points' 3-D coordinates.

    Args:
    * cord_1: 3-D coordinate of the 1st point
    * cord_2: 3-D coordinate of the 2nd point
    * cord_3: 3-D coordinate of the 3rd point
    * cord_4: 3-D coordinate of the 4th point

    Returns:
    * rad: dihedral angle (in radian)
    """

    eps = 1e-6
    a1 = cord_2 - cord_1
    a2 = cord_3 - cord_2
    a3 = cord_4 - cord_3
    v1 = np.cross(a1, a2)
    v1 = v1 / np.sqrt((v1 * v1).sum(-1) + eps)
    v2 = np.cross(a2, a3)
    v2 = v2 / np.sqrt((v2 * v2).sum(-1) + eps)
    sign = np.sign((v1 * a3).sum(-1))
    rad = np.arccos(np.clip(
        (v1 * v2).sum(-1) / np.sqrt((v1 ** 2).sum(-1) * (v2 ** 2).sum(-1) + eps), -1.0, 1.0))
    if sign != 0:
        rad *= sign

    return rad

def calc_plane_angle(cord_1, cord_2, cord_3):
    """Calculate the plane angle defined by 3 points' 3-D coordinates.

    Args:
    * cord_1: 3-D coordinate of the 1st point
    * cord_2: 3-D coordinate of the 2nd point
    * cord_3: 3-D coordinate of the 3rd point

    Returns:
    * rad: planar angle (in radian)
    """

    eps = 1e-6
    a1 = cord_1 - cord_2
    a2 = cord_3 - cord_2
    rad = np.arccos(np.clip(
        np.dot(a1, a2) / (np.linalg.norm(a1) * np.linalg.norm(a2) + eps), -1.0, 1.0))

    return rad


def parse_lbl_file(path, nctc_pos='first', da_bins = [40,25], dist_min=0.0, dist_max=40.0):
    """Parse the NPZ file containing GT-labels for inter-residue distance & orientation."""

    # configurations
    n_bins_dist = da_bins[0]
    n_bins_angle = da_bins[1]

    # functions for building classification labels
    def _build_dist_idxs(dist_mat, n_bins, dist_min=dist_min, dist_max=dist_max):
        bin_wid = (dist_max - dist_min) / (n_bins - 1)
        idxs = np.clip(np.floor((dist_mat - dist_min) / bin_wid).astype(np.int64), 0, n_bins - 1)
        nctc_mat = (idxs == n_bins - 1).astype(np.int8)
        return idxs, nctc_mat

    def _build_dihd_idxs(angl_mat, nctc_mat, n_bins, angl_min=-np.pi, angl_max=np.pi):
        bin_wid = (angl_max - angl_min) / (n_bins - 1)
        idxs = np.clip(np.floor((angl_mat - angl_min) / bin_wid).astype(np.int64), 0, n_bins - 2)
        idxs[nctc_mat == 1] = n_bins - 1
        return idxs

    def _build_plan_idxs(angl_mat, nctc_mat, n_bins, angl_min=-np.pi, angl_max=np.pi):
        bin_wid = (angl_max - angl_min) / (n_bins - 1)
        idxs = np.clip(np.floor((angl_mat - angl_min) / bin_wid).astype(np.int64), 0, n_bins - 2)
        idxs[nctc_mat == 1] = n_bins - 1
        return idxs

    # parse the NPZ file
    with np.load(path) as npz_data:
        # build classification labels (non-contact bin is the last one)
        idxs_p, nctc_mat = _build_dist_idxs(npz_data['p-val'], n_bins_dist)
        idxs_c4_, _ = _build_dist_idxs(npz_data['c4_-val'], n_bins_dist)
        idxs_n, _ = _build_dist_idxs(npz_data['n-val'], n_bins_dist)

        idxs_om = _build_dihd_idxs(npz_data['om-val'], nctc_mat, n_bins_angle)
        idxs_th = _build_dihd_idxs(npz_data['th-val'], nctc_mat, n_bins_angle)
        idxs_ph = _build_dihd_idxs(npz_data['ph-val'], nctc_mat, n_bins_angle)

        # move the non-contact bin to the first one if needed
        assert nctc_pos in ['first', 'last'], 'unrecognized <nctc_pos>: ' + nctc_pos
        if nctc_pos == 'first':
            idxs_p = (idxs_p + 1) % n_bins_dist  # move the non-contact bin to the first one
            idxs_c4_ = (idxs_c4_ + 1) % n_bins_dist  # move the non-contact bin to the first one
            idxs_n = (idxs_n + 1) % n_bins_dist  # move the non-contact bin to the first one

            idxs_om = (idxs_om + 1) % n_bins_angle
            idxs_th = (idxs_th + 1) % n_bins_angle
            idxs_ph = (idxs_ph + 1) % n_bins_angle

        # pack into a dict
        data_dict = {
            'p-idx': torch.tensor(idxs_p, dtype=torch.int64),  # L x L
            'p-msk': torch.tensor(npz_data['p-msk'], dtype=torch.int8),  # L x L
            'c4_-idx': torch.tensor(idxs_c4_, dtype=torch.int64),  # L x L
            'c4_-msk': torch.tensor(npz_data['c4_-msk'], dtype=torch.int8),  # L x L
            'n-idx': torch.tensor(idxs_n, dtype=torch.int64),  # L x L
            'n-msk': torch.tensor(npz_data['n-msk'], dtype=torch.int8),  # L x L

            'om-idx': torch.tensor(idxs_om, dtype=torch.int64),  # L x L
            'om-msk': torch.tensor(npz_data['om-msk'], dtype=torch.int8),  # L x L
            'th-idx': torch.tensor(idxs_th, dtype=torch.int64),  # L x L
            'th-msk': torch.tensor(npz_data['th-msk'], dtype=torch.int8),  # L x L
            'ph-idx': torch.tensor(idxs_ph, dtype=torch.int64),  # L x L
            'ph-msk': torch.tensor(npz_data['ph-msk'], dtype=torch.int8),  # L x L
        }

    return data_dict

def build_npz_file_deeprnafold(fasta_fpath, pdb_fpath, npz_fpath):
    """Build a NPZ file w/ ground-truth labels for inter-residue distance & orientation."""

    if os.path.exists(npz_fpath):
        # logging.info('NPZ file exists: %s', npz_fpath)
        return

    # initialization
    # full atom pdb parser version 'v0.0', deeprnafold frame version 'v3.0' ["P", "C4'", 'N9',]
    pdb_version = 'v0.0'
    frame_version = 'v3.0'

    parser = PdbParser(pdb_version=pdb_version, frame_version=frame_version, check_mode='lenient')

    aa_seq, atom_cords, atom_masks, structure, error_msg = \
        parser.run(pdb_fpath, fas_fpath=fasta_fpath)

    n_resds = len(aa_seq)

    x_p, x_c4_, x_n = [x.squeeze_().numpy() for x in torch.split(torch.FloatTensor(atom_cords), 1, dim=1)]
    m_p, m_c4_, m_n = [x.squeeze_().numpy() for x in torch.split(torch.FloatTensor(atom_masks), 1, dim=1)]

    # build ground-truth labels
    labl_data = {}
    labl_data['p-val'] = cdist(x_p, x_p).astype(np.float16)
    labl_data['p-msk'] = np.outer(m_p, m_p).astype(np.int8)
    labl_data['c4_-val'] = cdist(x_c4_, x_c4_).astype(np.float16)
    labl_data['c4_-msk'] = np.outer(m_c4_, m_c4_).astype(np.int8)
    labl_data['n-val'] = cdist(x_n, x_n).astype(np.float16)
    labl_data['n-msk'] = np.outer(m_n, m_n).astype(np.int8)

    labl_data['om-val'] = np.zeros((n_resds, n_resds), dtype=np.float16)
    labl_data['om-msk'] = np.zeros((n_resds, n_resds), dtype=np.int8)
    labl_data['th-val'] = np.zeros((n_resds, n_resds), dtype=np.float16)
    labl_data['th-msk'] = np.zeros((n_resds, n_resds), dtype=np.int8)
    labl_data['ph-val'] = np.zeros((n_resds, n_resds), dtype=np.float16)
    labl_data['ph-msk'] = np.zeros((n_resds, n_resds), dtype=np.int8)

    # build ground-truth labels
    for ir in tqdm(range(n_resds)):
        for ic in range(n_resds):
            labl_data['om-val'][ir, ic] = calc_dihedral_angle(x_c4_[ir], x_n[ir], x_n[ic], x_c4_[ic])
            labl_data['om-msk'][ir, ic] = m_c4_[ir] * m_n[ir] * m_n[ic] * m_c4_[ic]
            labl_data['th-val'][ir, ic] = calc_dihedral_angle(x_p[ir], x_c4_[ir], x_n[ir], x_n[ic])
            labl_data['th-msk'][ir, ic] = m_p[ir] * m_c4_[ir] * m_n[ir] * m_n[ic]
            labl_data['ph-val'][ir, ic] = calc_dihedral_angle(x_p[ic], x_c4_[ic], x_n[ic], x_n[ir])
            labl_data['ph-msk'][ir, ic] = m_p[ic] * m_c4_[ic] * m_n[ic] * m_n[ir]

    is_show = False
    if is_show:
        plt.figure(figsize=(12, 9))
        plt.subplot(2, 3, 1)
        plt.imshow(labl_data['p-val'].astype(np.float32))
        plt.title('Predicted Distance')
        plt.subplot(2, 3, 2)
        plt.imshow(labl_data['c4_-val'].astype(np.float32))
        plt.title('Predicted Distance')
        plt.subplot(2, 3, 3)
        plt.imshow(labl_data['n-val'].astype(np.float32))
        plt.title('Predicted Distance')

        plt.subplot(2, 3, 4)
        plt.imshow(labl_data['om-val'].astype(np.float32))
        plt.title('Predicted Distance')
        plt.subplot(2, 3, 5)
        plt.imshow(labl_data['th-val'].astype(np.float32))
        plt.title('Predicted Distance')
        plt.subplot(2, 3, 6)
        plt.imshow(labl_data['ph-val'].astype(np.float32))
        plt.title('Predicted Distance')
        plt.show()

    # build a NPZ file
    os.makedirs(os.path.dirname(os.path.realpath(npz_fpath)), exist_ok=True)
    np.savez(npz_fpath, **labl_data)
    logging.info('lbl file built: %s', npz_fpath)

def build_npz_file_deeprnafold_cplx(fasta_fpath, pdb_fpath, npz_fpath):
    """Build a NPZ file w/ ground-truth labels for inter-residue distance & orientation."""

    assert isinstance(fasta_fpath, list)
    assert isinstance(pdb_fpath, list)

    if os.path.exists(npz_fpath):
        # logging.info('NPZ file exists: %s', npz_fpath)
        return

    # logging.info('NPZ file exists: %s', npz_fpath)

    # initialization
    # full atom pdb parser version 'v0.0', deeprnafold frame version 'v3.0' ["P", "C4'", 'N9',]
    pdb_version = 'v0.0'
    frame_version = 'v3.0'

    parser = PdbParser(pdb_version=pdb_version, frame_version=frame_version, check_mode='lenient')

    aa_seqs = []
    atom_cords_list = []
    atom_masks_list = []

    for fas, pdb in zip(fasta_fpath, pdb_fpath):
        aa_seq, atom_cords, atom_masks, structure, error_msg = parser.run(pdb, fas_fpath=fas)
        aa_seqs.append(aa_seq)
        atom_cords_list.append(atom_cords)
        atom_masks_list.append(atom_masks)

    aa_seq = '_'.join(aa_seqs)
    seq_len = len(aa_seq)

    cords = np.zeros([seq_len, parser.frame_info.ATOM_NUM_MAX, 3])
    cords_mask = np.zeros([seq_len, parser.frame_info.ATOM_NUM_MAX])
    i_s = 0
    for i, (s, c, cmsk) in enumerate(zip(aa_seqs, atom_cords_list, atom_masks_list)):
        c = c.reshape([len(s), parser.frame_info.ATOM_NUM_MAX, 3])
        cmsk = cmsk.reshape([len(s), parser.frame_info.ATOM_NUM_MAX])
        i_e = i_s + len(s)
        cords[i_s:i_e, :, :] = c
        cords_mask[i_s:i_e, :] = cmsk
        i_s = i_e + 1

    atom_cords = cords
    atom_masks = cords_mask

    n_resds = len(aa_seq)
    x_p, x_c4_, x_n = [x.squeeze_().numpy() for x in torch.split(torch.FloatTensor(atom_cords), 1, dim=1)]
    m_p, m_c4_, m_n = [x.squeeze_().numpy() for x in torch.split(torch.FloatTensor(atom_masks), 1, dim=1)]

    # build ground-truth labels
    labl_data = {}
    labl_data['p-val'] = cdist(x_p, x_p).astype(np.float16)
    labl_data['p-msk'] = np.outer(m_p, m_p).astype(np.int8)
    labl_data['c4_-val'] = cdist(x_c4_, x_c4_).astype(np.float16)
    labl_data['c4_-msk'] = np.outer(m_c4_, m_c4_).astype(np.int8)
    labl_data['n-val'] = cdist(x_n, x_n).astype(np.float16)
    labl_data['n-msk'] = np.outer(m_n, m_n).astype(np.int8)

    labl_data['om-val'] = np.zeros((n_resds, n_resds), dtype=np.float16)
    labl_data['om-msk'] = np.zeros((n_resds, n_resds), dtype=np.int8)
    labl_data['th-val'] = np.zeros((n_resds, n_resds), dtype=np.float16)
    labl_data['th-msk'] = np.zeros((n_resds, n_resds), dtype=np.int8)
    labl_data['ph-val'] = np.zeros((n_resds, n_resds), dtype=np.float16)
    labl_data['ph-msk'] = np.zeros((n_resds, n_resds), dtype=np.int8)

    # build ground-truth labels
    for ir in tqdm(range(n_resds)):
        for ic in range(n_resds):
            labl_data['om-val'][ir, ic] = calc_dihedral_angle(x_c4_[ir], x_n[ir], x_n[ic], x_c4_[ic])
            labl_data['om-msk'][ir, ic] = m_c4_[ir] * m_n[ir] * m_n[ic] * m_c4_[ic]
            labl_data['th-val'][ir, ic] = calc_dihedral_angle(x_p[ir], x_c4_[ir], x_n[ir], x_n[ic])
            labl_data['th-msk'][ir, ic] = m_p[ir] * m_c4_[ir] * m_n[ir] * m_n[ic]
            labl_data['ph-val'][ir, ic] = calc_dihedral_angle(x_p[ic], x_c4_[ic], x_n[ic], x_n[ir])
            labl_data['ph-msk'][ir, ic] = m_p[ic] * m_c4_[ic] * m_n[ic] * m_n[ir]

    is_show = False
    if is_show:
        plt.figure(figsize=(12, 9))
        plt.subplot(2, 3, 1)
        plt.imshow(labl_data['p-val'].astype(np.float32))
        plt.title('Predicted Distance')
        plt.subplot(2, 3, 2)
        plt.imshow(labl_data['c4_-val'].astype(np.float32))
        plt.title('Predicted Distance')
        plt.subplot(2, 3, 3)
        plt.imshow(labl_data['n-val'].astype(np.float32))
        plt.title('Predicted Distance')

        plt.subplot(2, 3, 4)
        plt.imshow(labl_data['om-val'].astype(np.float32))
        plt.title('Predicted Distance')
        plt.subplot(2, 3, 5)
        plt.imshow(labl_data['th-val'].astype(np.float32))
        plt.title('Predicted Distance')
        plt.subplot(2, 3, 6)
        plt.imshow(labl_data['ph-val'].astype(np.float32))
        plt.title('Predicted Distance')
        plt.show()

    # build a NPZ file
    os.makedirs(os.path.dirname(os.path.realpath(npz_fpath)), exist_ok=True)
    np.savez(npz_fpath, **labl_data)
    logging.info('lbl file built: %s', npz_fpath)


def gen_one(param):
    fasta_fpath, pdb_fpath, npz_fpath = param
    try:
        build_npz_file_deeprnafold(fasta_fpath, pdb_fpath, npz_fpath)
    except:
        print('sth is wrong')


def gen_one_cplx(param):
    fasta_fpath, pdb_fpath, npz_fpath = param
    try:
        build_npz_file_deeprnafold_cplx(fasta_fpath, pdb_fpath, npz_fpath)
    except:
        print('sth is wrong')

def gen_pdb_da_lbls():

    root = '/public/home/taoshen/data/rna/RNA3D/Ash_Total_RNA_test/'

    zfold_init()

    params = []

    tag = 'pad'
    tgts = os.listdir(f'{root}/rna_pdb_align_{tag}/')

    random.shuffle(tgts)
    for tgt in tgts:
        if tgt.endswith('.pdb_out'):
            tgt = tgt.replace('.pdb_out', '')
            fasta_fpath = f'{root}/rna_seq/{tgt}.seq'
            pdb_fpath = f'{root}/rna_pdb_align_{tag}/{tgt}.pdb_out'
            npz_fpath = f'{root}/da_npzs_{tag}/{tgt}.npz'
            params.append((fasta_fpath, pdb_fpath, npz_fpath))

    print(len(params))
    p_map(gen_one, params)

def gen_semi_da_lbls():

    for mode in ['standard_model_1_unrefined_align', 'standard_model_1_refined_align']:
        root = '/public/home/taoshen/data/rna/bpRNA/'

        zfold_init()
        params = []
        for tgt in os.listdir(f'{root}/{mode}/'):
            if tgt.endswith('.pdb_out'):
                tgt = tgt.replace('.pdb_out', '')
                fasta_fpath = f'{root}/seq/{tgt}.seq'
                pdb_fpath = f'{root}/{mode}/{tgt}.pdb_out'
                npz_fpath = f'{root}/{mode}_da_npzs/{tgt}.npz'
                params.append((fasta_fpath, pdb_fpath, npz_fpath))
        print(len(params))
        p_map(gen_one, params)

def gen_pdb_cplx_da_lbls():

    cplx_ids = r'/public/home/taoshen/data/rna/RNA3D/Ash_Total_RNA_test/exp/paper1.0/valid_cplx_ids'
    cplx_ids = read_tgts(cplx_ids)

    root = '/public/home/taoshen/data/rna/RNA3D/Ash_Total_RNA_test/'

    zfold_init()

    params = []

    tag = 'pad'

    def get_split(t):
        t = t.split('_')
        pdb_id = t[0]
        return [f'{pdb_id}{t_}' for t_ in t[1:]]

    def merge_split(ts):
        chains = []
        for t in ts:
            pdb = t[:4]
            chain = t[4:]
            chains.append(chain)
        _all = [pdb] + chains
        return '_'.join(_all)

    for tgt in cplx_ids:
        split_tgt = get_split(tgt)
        fasta_fpath = []
        pdb_fpath = []
        # random.shuffle(split_tgt)
        # split_tgt.reverse()
        for st in split_tgt:
            fasta_fpath.append(f'{root}/rna_seq/{st}.seq')
            pdb_fpath.append(f'{root}/rna_pdb_align_{tag}/{st}.pdb_out')

        new_tgt = merge_split(split_tgt)
        npz_fpath = f'{root}/da_npzs_{tag}/{new_tgt}.npz'
        params.append((fasta_fpath, pdb_fpath, npz_fpath))

    print(len(params))
    p_map(gen_one_cplx, params)


def _convert_pseudo_lbl(param):

    '''
    convert output to label
    '''

    npz_in, npz_out = param

    if os.path.exists(npz_out):
        return
    if not os.path.exists(npz_in):
        return

    try:
        data = np.load(npz_in)
        data = dict(data)

        # build ground-truth labels
        labl_data = {}
        labl_data['p-val'] = np.argmax(data['dist_p'], axis=-1).astype(np.float16)
        labl_data['p-msk'] = np.ones_like(labl_data['p-val'])
        labl_data['c4_-val'] = np.argmax(data['dist_c4'], axis=-1).astype(np.float16)
        labl_data['c4_-msk'] = np.ones_like(labl_data['c4_-val'])
        labl_data['n-val'] = np.argmax(data['dist_n'], axis=-1).astype(np.float16)
        labl_data['n-msk'] = np.ones_like(labl_data['n-val'])

        labl_data['om-val'] = np.argmax(data['omega'], axis=-1).astype(np.float16)
        labl_data['om-msk'] = np.ones_like(labl_data['om-val'])
        labl_data['th-val'] = np.argmax(data['theta'], axis=-1).astype(np.float16)
        labl_data['th-msk'] = np.ones_like(labl_data['th-val'])
        labl_data['ph-val'] = np.argmax(data['phi'], axis=-1).astype(np.float16)
        labl_data['ph-msk'] = np.ones_like(labl_data['ph-val'])

        # build a NPZ file
        os.makedirs(os.path.dirname(os.path.realpath(npz_out)), exist_ok=True)
        np.savez(npz_out, **labl_data)
        logging.info('lbl file built: %s', npz_out)

    except:
        return


def convert_pseudo_lbls():
    # _root = '/public/home/taoshen/data/rMSA_gen'
    # root = f'{_root}/drfold_npzs'
    # save_root = f'{_root}/pseudo_danpzs'

    _root = '/public/home/taoshen/data/rna/bpRNA'
    root = f'{_root}/npzs'
    save_root = f'{_root}/pseudo_danpzs'

    os.makedirs(save_root, exist_ok=True)

    params = []
    for tgt in tqdm(os.listdir(root)):
        npz_in = f'{root}/{tgt}/drfold.npz'
        npz_out = f'{save_root}/{tgt}.npz'
        params.append([npz_in, npz_out])

    p_map(_convert_pseudo_lbl, params)


def gen_10cv():

    root = '/public/home/taoshen/data/rna/RNA3D/mmcif2pdb_rna/'
    npz_dir = f'{root}/da_npzs'
    os.makedirs(npz_dir, exist_ok=True)

    zfold_init()
    params = []
    tgts = os.listdir(f'{root}/pdb/')
    random.shuffle(tgts)

    for tgt in tgts:
        if tgt.endswith('.pdb'):
            tgt = tgt.replace('.pdb', '')
            fasta_fpath = f'{root}/seq/{tgt}.seq'
            pdb_fpath = f'{root}/pdb/{tgt}.pdb'
            npz_fpath = f'{root}/da_npzs/{tgt}.npz'

            params.append((fasta_fpath, pdb_fpath, npz_fpath))

    print(len(params))
    p_map(gen_one, params)

def gen_pdb_10cv_cplx_da_lbls():

    root = '/public/home/taoshen/data/rna/RNA3D/mmcif2pdb_rna'
    cplx_ids = r'/public/home/taoshen/code/Projects/E2EFold3D_Dev/data_pipeline/stats/nrlist_3.226_4.0A/' \
               r'nrlist_3.226_4.0A_cplx'
    cplx_ids = read_tgts(cplx_ids)

    zfold_init()
    params = []

    def get_split(t):
        t = t.split('_')
        pdb_id = t[0]
        return [f'{pdb_id}_{t_}' for t_ in t[1:]]

    def merge_split(ts):
        chains = []
        for t in ts:
            pdb = t[:4]
            chain = t[5:]
            chains.append(chain)
        _all = [pdb] + chains
        return '_'.join(_all)

    for tgt in cplx_ids:
        split_tgt = get_split(tgt)
        fasta_fpath = []
        pdb_fpath = []

        split_tgt.reverse()

        for st in split_tgt:
            fasta_fpath.append(f'{root}/seq/{st}.seq')
            pdb_fpath.append(f'{root}/pdb/{st}.pdb')

        new_tgt = merge_split(split_tgt)
        print(new_tgt)

        npz_fpath = f'{root}/da_npzs/{new_tgt}.npz'
        params.append((fasta_fpath, pdb_fpath, npz_fpath))

    print(len(params))
    p_map(gen_one_cplx, params)

if __name__ == '__main__':
    # gen_pdb_da_lbls()
    # gen_pdb_da_lbls()
    # gen_semi_da_lbls()
    # gen_pdb_10cv_cplx_da_lbls()
    gen_10cv()

    # convert_pseudo_lbls()
    # gen_10cv()
    # plt.figure(figsize=(12, 9))
    # plt.subplot(2, 3, 1)
    # plt.imshow(np.argmax(data['dist_n'], axis=-1).astype(np.float32))
    # plt.title('Predicted Distance')
    # plt.subplot(2, 3, 2)
    # plt.imshow(np.argmax(data['dist_c4'], axis=-1).astype(np.float32))
    # plt.title('Predicted Distance')
    # plt.subplot(2, 3, 3)
    # plt.imshow(np.argmax(data['dist_p'], axis=-1).astype(np.float32))
    # plt.title('Predicted Distance')
    #
    # plt.subplot(2, 3, 4)
    # plt.imshow(np.argmax(data['omega'], axis=-1).astype(np.float32))
    # plt.title('Predicted Distance')
    # plt.subplot(2, 3, 5)
    # plt.imshow(np.argmax(data['theta'], axis=-1).astype(np.float32))
    # plt.title('Predicted Distance')
    # plt.subplot(2, 3, 6)
    # plt.imshow(np.argmax(data['phi'], axis=-1).astype(np.float32))
    # plt.title('Predicted Distance')
    #
    # plt.show()




