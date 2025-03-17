"""Dataset for proteins."""
import warnings
warnings.filterwarnings("ignore")

from rnafold.models.config import default_config
from rnafold.tools.dataset_helper import *
from drfold.tools.features import collect_features
from rnafold.tools.rna_constants import BSPAIR_CONSTRAINT

import logging

class DatasetConfig():
    """Configurations for the <Dataset> class."""

    def __init__(
            self,
            id_fpath,             # IDs file
            fas_dpath,            # directory path to FASTA files
            msa_dpath,            # directory path to MSA files
            is_complex=False,
            clstr_id_fpath=None,  # IDs file
            clstr_fpath=None,     # IDs file
            pdb_dpath=None,       # directory path to PDB files
            danpz_dpath=None,     # directory path to PDB files
            ss_dpath=None,        # directory path to ss npy files
            ff_dpath=None,
            tpl_dpath=None,
            dr_dpath=None,
            dr_fea_dpath=None,
            dr_ss_dpath=None,     # directory path to ss npy files
            dr_cm_dpath=None,     # directory path to ss npy files
            is_dr_fea=False,
            is_dr_pred=False,
            pdb_version='v0.0',
            frame_version='v4.3',
            crop_size=128,        # crop size
            msa_depth=64,         # msa_depth
            dist_config=None,
            num=None,
        ):
        """Constructor function."""

        # setup configurations
        self.id_fpath = id_fpath
        self.fas_dpath = fas_dpath
        self.msa_dpath = msa_dpath
        self.pdb_dpath = pdb_dpath
        self.ss_dpath = ss_dpath
        self.ff_dpath = ff_dpath
        self.tpl_dpath = tpl_dpath
        self.msa_depth = msa_depth

        self.danpz_dpath = danpz_dpath
        self.dr_ss_dpath = dr_ss_dpath
        self.dr_cm_dpath = dr_cm_dpath
        self.dr_dpath = dr_dpath
        self.dr_fea_dpath = dr_fea_dpath
        self.is_complex = is_complex
        self.dist_config = dist_config
        self.is_dr_fea = is_dr_fea
        self.is_dr_pred = is_dr_pred
        self.crop_size = crop_size
        self.pdb_version = pdb_version
        self.frame_version = frame_version

        self.rna_ids = None
        if self.id_fpath is not None:
            with open(self.id_fpath, 'r') as f:
                self.rna_ids = [tgt.strip() for tgt in f.readlines()]

        self.clstr_id_fpath = clstr_id_fpath
        self.clstr_fpath = clstr_fpath
        self.clstrs = None
        self.clstr_ids = None
        self.num = num

    def show(self):
        """Show detailed configurations."""
        logging.info('=== DatasetConfig - Start ===')
        logging.info('fas_dpath: %s', self.fas_dpath)
        logging.info('msa_dpath: %s', self.msa_dpath)
        logging.info('pdb_dpath: %s', self.pdb_dpath)
        logging.info('is_complex: %s', self.is_complex)
        logging.info('dr_dpath: %s', self.dr_dpath)
        logging.info('dr_ss_dpath: %s',  self.dr_ss_dpath)
        logging.info('dr_cm_dpath: %s',  self.dr_cm_dpath)
        logging.info('dr_fea_dpath: %s', self.dr_fea_dpath)
        logging.info('   is_fea: %s', self.is_dr_fea)
        logging.info('   is_pred: %s', self.is_dr_pred)

        logging.info('id_fpath: %s', self.id_fpath)
        logging.info('rna num: %d', len(self.rna_ids))
        logging.info('danpz_dpath : %s', self.danpz_dpath)
        logging.info('dist bins: %d', self.dist_config.bins)
        logging.info('     max_dist: %d', self.dist_config.max_dist)
        logging.info('     min dist: %d', self.dist_config.min_dist)
        logging.info('ss_dpath: %s', self.ss_dpath)
        logging.info('ff_dpath: %s', self.ff_dpath)
        logging.info('tpl_dpath: %s', self.tpl_dpath)
        logging.info('crop_size: %s', self.crop_size)
        logging.info('msa_depth: %s', self.msa_depth)
        logging.info('pdb version: %s', self.pdb_version)
        logging.info('frame version: %s', self.frame_version)
        logging.info('=== DatasetConfig - Finish ===')

class RNADataset(Dataset):
    """Dataset for proteins."""

    def __init__(self, config, is_train=True):
        """Constructor function."""

        super().__init__()

        # setup configurations
        assert isinstance(config, DatasetConfig)

        self.config = config
        self.config.show()

        self.is_train = is_train
        # additional configurations
        self.n_resds_max = 1000
        self.crop_size = self.config.crop_size

        # setup the pdb parser
        self.pdb_parser = PdbParser(pdb_version=config.pdb_version, frame_version=config.frame_version)

        self.n_grps_cord = self.pdb_parser.frame_info.ATOM_NUM_MAX

        self.rna_ids = None
        if self.config.is_complex:
            self.rna_ids = self.config.rna_ids
            self.num = len(self.rna_ids) * 1000 if is_train else len(self.rna_ids)
            if self.is_train:
                random.shuffle(self.rna_ids)
            else:
                self.rna_ids.sort()

        elif exists(self.config.rna_ids):
            self.rna_ids = []
            # check the existence of FASTA & NPZ files
            for rna_id in self.config.rna_ids:
                fas_fpath = os.path.join(self.config.fas_dpath, '%s.seq' % rna_id)
                if os.path.exists(fas_fpath):
                    self.rna_ids.append(rna_id)
                else:
                    logging.warning('FASTA file is missing (ID: %s)', rna_id)

            self.num = len(self.rna_ids) * 1000 if is_train else len(self.rna_ids)
            if self.is_train:
                random.shuffle(self.rna_ids)
            else:
                self.rna_ids.sort()

        self.sizes = np.asarray([self.crop_size for _ in range(self.num)])
        self.msa_mode = 'Rand' if self.is_train else 'TopN'

    def __len__(self):
        """Get the number of elements in the dataset."""

        return self.num

    def __getitem__(self, idx):
        """Get the i-th element in the dataset."""

        # build a data dict of core data
        rna_id = self.rna_ids[idx % len(self.rna_ids)]

        try:

            if self.config.is_complex:
                core_data = self.__build_core_data_cplx(rna_id)
                data_dict = self.__build_data_dict(core_data)
                # do not support drfold
                # data_dict.update(self.__build_core_data_drfold_cplx(rna_id))
            else:
                core_data = self.__build_core_data(rna_id)
                data_dict = self.__build_data_dict(core_data)
                data_dict.update(self.__build_core_data_drfold(rna_id))

            # random crop
            data_dict = self.__get_crop_data(data_dict)
            # build input tensors
            data_dict = self.__build_ss_constraint(data_dict)
            return data_dict

        except:
            print(f'sth is wrong sample {idx} {rna_id}')
            return self.__getitem__(np.random.randint(0, self.num - 1))

    def __get_pdb_fpath(self, rna_id):
        # obtain pdb_fpath
        if self.config.pdb_version == 'v0.0':
            pdb_fpath = os.path.join(self.config.pdb_dpath, '%s.pdb' % rna_id) \
                if exists(self.config.pdb_dpath) else None
            if not os.path.exists(pdb_fpath):
                pdb_fpath = os.path.join(self.config.pdb_dpath, '%s.pdb_out' % rna_id) \
                    if exists(self.config.pdb_dpath) else None
        else:
            raise NotImplementedError

        return pdb_fpath

    def __get_msa_fpath(self, rna_id):

        fas_fpath = os.path.join(self.config.fas_dpath, '%s.seq' % rna_id)

        if exists(self.config.msa_dpath) and isinstance(self.config.msa_dpath, list):
            random.shuffle(self.config.msa_dpath)
            msa_fpath = os.path.join(self.config.msa_dpath[0], '%s.a3m' % rna_id)
        elif exists(self.config.msa_dpath):
            msa_fpath = os.path.join(self.config.msa_dpath, '%s.a3m' % rna_id)

        if not os.path.exists(msa_fpath):
            msa_fpath = fas_fpath

        return fas_fpath, msa_fpath

    def __build_core_data_drfold(self, rna_id):

        drfold_dict = {}
        drfold_dict['msa'] = None
        drfold_dict['seq'] = None
        drfold_dict['hmm'] = None
        drfold_dict['ss'] = None
        drfold_dict['seq_fea'] = None
        drfold_dict['pair_fea'] = None
        drfold_dict['da_lbls_soft'] = None

        if exists(self.config.dr_ss_dpath) and exists(self.config.dr_cm_dpath):
            # only fasta for drfold inputs
            fas_fpath, msa_fpath = self.__get_msa_fpath(rna_id)
            hmm_file = f'{self.config.dr_cm_dpath}/{rna_id}.cm'
            ss_file = f'{self.config.dr_ss_dpath}/{rna_id}.txt'

            if self.is_train:
                if random.randint(0, 1) == 0:
                    hmm_file = './tmp.cm'
                if random.randint(0, 1) == 0:
                    ss_file = './tmp.txt'

            if not self.is_train and not os.path.exists(hmm_file):
                print('hmm file does not exists')
            if not self.is_train and not os.path.exists(ss_file):
                print('ss file does not exists')

            drfold_dict.update(collect_features(fas_fpath, msa_fpath, hmm_file, ss_file, self.config.msa_depth))

        if exists(self.config.dr_fea_dpath):
            data_fea = {}
            npz_file = f'{self.config.dr_fea_dpath}/{rna_id}/drfold.npz'
            data_pred = dict(np.load(npz_file))
            data_pred = {k:torch.FloatTensor(data_pred[k]) for k in data_pred.keys()}
            data_fea['da_lbls_soft'] = data_pred

            if self.config.is_dr_fea:
                npz_file = f'{self.config.dr_fea_dpath}/{rna_id}/drfold_fea.npz'
                data_fea_ = dict(np.load(npz_file))
                for k in data_fea_.keys():
                    data_fea_[k] = torch.FloatTensor(data_fea_[k])
                data_fea.update(data_fea_)

            if self.config.is_dr_pred:
                pair_preds, seq_preds = [], []

                npz_file = f'{self.config.dr_fea_dpath}/{rna_id}/drfold.npz'
                data_pred = dict(np.load(npz_file))
                pair_pred = [data_pred['dist_n'], data_pred['dist_c4'], data_pred['dist_p'],
                             data_pred['omega'], data_pred['theta'], data_pred['phi']]
                seq_pred = [data_pred['eta_bb'], data_pred['theta_bb']]

                pair_pred = [torch.FloatTensor(t) for t in pair_pred]
                seq_pred = [torch.FloatTensor(t) for t in seq_pred]
                pair_pred = torch.cat(pair_pred, dim=-1).unsqueeze(0)
                seq_pred = torch.cat(seq_pred, dim=-1).unsqueeze(0)
                pair_preds.append(pair_pred)
                seq_preds.append(seq_pred)
                pair_pred = torch.mean(torch.cat(pair_preds, dim=0), dim=0, keepdim=True)
                seq_pred = torch.mean(torch.cat(seq_preds, dim=0), dim=0, keepdim=True)

                if 'seq_fea' in data_fea:
                    data_fea['seq_fea'] = torch.cat([data_fea['seq_fea'], seq_pred], dim=-1)
                else:
                    data_fea['seq_fea'] = seq_pred

                if 'pair_fea' in data_fea:
                    data_fea['pair_fea'] = torch.cat([data_fea['pair_fea'], pair_pred], dim=-1)
                else:
                    data_fea['pair_fea'] = pair_pred

            drfold_dict.update(data_fea)

        return {f'dr_{key}': drfold_dict[key] for key in drfold_dict.keys()}

    def __build_core_data_cplx(self, rna_complex_id):
        """Build the dict of core data."""

        # print('rna_complex_id', rna_complex_id)

        def merge_split(ts):
            chains = []
            for t in ts:
                pdb = t[:4]
                chain = t[5:]
                chains.append(chain)
            _all = [pdb] + chains
            return '_'.join(_all)

        rna_ids = get_split_ids(rna_complex_id) # do not shuffle

        if self.is_train and random.randint(0,1)==0:
            rna_ids.reverse()

        rna_complex_id = merge_split(rna_ids)

        if len(rna_ids) == 1:
            return self.__build_core_data(rna_ids[0])

        core_data_dicts = [self.__build_core_data(rna_id) for rna_id in rna_ids]
        aa_seq_list = [d['seq'] for d in core_data_dicts]
        cords_list = [d['cord_gt'] for d in core_data_dicts]
        cords_mask_list = [d['cord_gt_mask'] for d in core_data_dicts]
        msa_tokens_list = [d['token'] for d in core_data_dicts]
        rna_fm_tokens_list = [d['rna_fm_token'] for d in core_data_dicts]

        # complex emsemble
        aa_seq = '-'.join(aa_seq_list)
        seq_len = len(aa_seq)
        cords = np.zeros([seq_len, self.n_grps_cord, 3])
        cords_mask = np.zeros([seq_len, self.n_grps_cord])

        complex_msa_depth = min([t.shape[0] for t in msa_tokens_list])
        msa_tokens = torch.ones([complex_msa_depth, seq_len], dtype=torch.int64) * alphabet.padding_idx
        rna_fm_tokens = torch.ones([seq_len], dtype=torch.int64) * rnafm_alphabet.padding_idx

        i_s = 0
        for i, (s, c, cmsk, m_t, r_t) in enumerate(zip(aa_seq_list,
                                                       cords_list, cords_mask_list,
                                                       msa_tokens_list, rna_fm_tokens_list)):

            c = c.reshape([len(s), self.n_grps_cord, 3]) if c is not None else None
            cmsk = cmsk.reshape([len(s), self.n_grps_cord]) if c is not None else None

            i_e = i_s + len(s)

            if c is not None:
                cords[i_s:i_e, :, :] = c
                cords_mask[i_s:i_e, :] = cmsk

            msa_tokens[:, i_s:i_e] = m_t[:complex_msa_depth]
            rna_fm_tokens[i_s:i_e] = r_t
            i_s = i_e + 1

        assert rna_fm_tokens.shape[-1] == msa_tokens.shape[-1]

        cords = cords.reshape([1, -1, 3])
        cords_mask = cords_mask.reshape([1, -1])

        assert cords.shape[1] == seq_len * self.n_grps_cord, print(cords.shape, seq_len * self.n_grps_cord)

        # preprocessed cplex da lbls
        da_dict = None
        if exists(self.config.danpz_dpath):
            danpz_fpath = os.path.join(self.config.danpz_dpath, '%s.npz' % rna_complex_id)
            da_dict = self.__build_da_lbls(danpz_fpath)

        # pack all the essential data into dict
        core_data = {
            'id': rna_complex_id,
            'seq': aa_seq,
            'token': msa_tokens,
            'rna_fm_token': rna_fm_tokens,
            't1ds': None,
            't2ds': None,
            'cmap': None,
            'cmap_msk': None,
            'da_lbls': da_dict,
            'cord_gt': cords,
            'cord_gt_mask': cords_mask,
        }
        return core_data

    def __build_core_data(self, rna_id):
        """Build the dict of core data."""

        # initialization
        t1d, t2d = None, None
        cmap, cmap_msk = None, None

        # obtain the amino-acid sequence
        fas_fpath = os.path.join(self.config.fas_dpath, '%s.seq' % rna_id)
        _, aa_seq = parse_fas_file(fas_fpath)

        # xxx
        seq_len = len(aa_seq)
        msa_tokens, rna_fm_tokens = self.__build_seq_feas(rna_id)

        # obtain ss
        if exists(self.config.ss_dpath):
            ss_fpath = os.path.join(self.config.ss_dpath, '%s.npy' % rna_id)
            cmap, cmap_msk = self.__build_ss_lbls(ss_fpath=ss_fpath)

        # obtain 3D coordinates
        cords, cords_mask = None, None
        if exists(self.config.pdb_dpath):
            pdb_fpath = self.__get_pdb_fpath(rna_id)

            if not os.path.exists(pdb_fpath):
                pdb_fpath = pdb_fpath.replace('refined', 'unrefined')

            cords, cords_mask = self.__build_pdb_lbls(pdb_fpath, fas_fpath)

            if 'unrefined' in pdb_fpath and exists(cords_mask):
                # print(f'for unrefined pdb, use only backbone info {rna_id}')
                cords_mask = cords_mask.reshape([seq_len, self.n_grps_cord])
                cords_mask[:,3:] = 0
                cords_mask = cords_mask.reshape([1, seq_len*self.n_grps_cord])

        # build distance/orientations lbls
        da_dict = None
        if exists(self.config.danpz_dpath):
            danpz_fpath = os.path.join(self.config.danpz_dpath, '%s.npz' % rna_id)
            da_dict = self.__build_da_lbls(danpz_fpath)

        assert rna_fm_tokens.shape[0] == msa_tokens.shape[-1]

        if exists(cords):
            assert cords.shape[1] == seq_len * self.n_grps_cord, \
                print(cords.shape, seq_len * self.n_grps_cord)

        if cmap is not None and seq_len != cmap.shape[0]:
            print(seq_len, cmap.shape[0], 'ss size does not match')
            cmap, cmap_msk = None, None

        # pack all the essential data into dict
        core_data = {
            'id': rna_id,
            'seq': aa_seq,
            'token': msa_tokens,
            'rna_fm_token': rna_fm_tokens,
            't1ds': t1d,
            't2ds': t2d,
            'cmap': cmap,
            'cmap_msk': cmap_msk,
            'da_lbls': da_dict,
            'cord_gt': cords,
            'cord_gt_mask': cords_mask,
        }
        return core_data

    def __build_seq_feas(self, rna_id):

        # obtain the rna msa

        fas_fpath, msa_fpath = self.__get_msa_fpath(rna_id)

        try:
            msa_tokens = get_msa_feature(msa_path=msa_fpath, msa_depth=self.config.msa_depth, mode=self.msa_mode)
        except:
            msa_tokens = get_msa_feature(msa_path=fas_fpath, msa_depth=self.config.msa_depth, mode=self.msa_mode)

        rna_fm_tokens = get_rna_fm_token(fas_fpath)

        if rna_fm_tokens.shape[0] != msa_tokens.shape[-1]:
            msa_tokens = get_msa_feature(msa_path=fas_fpath, msa_depth=self.config.msa_depth, mode=self.msa_mode)

        return msa_tokens, rna_fm_tokens

    def __build_tpl_feas(self, pkl_fpath, fas_path):
        tpl_topk = 4
        t1ds, t2ds, align_seq_tokens = \
            parse_tpl_file(pkl_fpath, fas_path, tpl_topk, self.is_train)
        return t1ds, t2ds, align_seq_tokens

    def __build_ss_lbls(self, ss_fpath):
        if ss_fpath is None:
            return None, None
        if not os.path.exists(ss_fpath):
            return None, None
        cmap = np.load(ss_fpath)
        cmap_msk = np.ones(cmap.shape)

        return cmap, cmap_msk

    def __build_ss_constraint(self, data_dict):
        ''' build ss constraint '''

        data_dict['ssc_i'] = None
        data_dict['ssc_j'] = None
        data_dict['ssc_dist'] = None

        if data_dict['cmap'] is None:
            return data_dict

        cmap = data_dict['cmap']
        aa_seq = data_dict['seq']

        index_i, index_j = [],[]
        pair_dist = []
        for i in range(cmap.shape[0]):
            for j in range(i, cmap.shape[0]):
                pair = f'{aa_seq[i]}{aa_seq[j]}'
                if cmap[i][j] == 1 and pair in BSPAIR_CONSTRAINT.keys():
                    index_i.append(i)
                    index_j.append(j)
                    pair_dist.append(BSPAIR_CONSTRAINT[pair])

        if len(pair_dist) > 0:
            data_dict['ssc_i'] = torch.LongTensor(np.array(index_j))
            data_dict['ssc_j'] = torch.LongTensor(np.array(index_i))
            data_dict['ssc_dist'] = torch.FloatTensor(np.array(pair_dist))

        return data_dict

    def __build_pdb_lbls(self, pdb_fpath, fas_fpath):

        if not exists(pdb_fpath) or not os.path.exists(pdb_fpath):
            print(f'{pdb_fpath} not exists!!!')
            return None, None

        # parse the PDB file
        aa_seq, cords, cords_mask, _, _ = self.pdb_parser.run(pdb_fpath, fas_fpath=fas_fpath)
        cords = cords.reshape([1, -1, 3])
        cords_mask = cords_mask.reshape([1, -1])
        return cords, cords_mask

    def __build_da_lbls(self, danpz_fpath):

        # if not os.path.exists(danpz_fpath):
        #     print(danpz_fpath)

        da_data_dict = parse_lbl_file(danpz_fpath,
                                      nctc_pos='first',
                                      da_bins=[self.config.dist_config.bins, 25],
                                      dist_min=self.config.dist_config.min_dist,
                                      dist_max=self.config.dist_config.max_dist)

        return da_data_dict

    def __build_data_dict(self, core_data):
        """Build a data dict from the core data (id + seq + cord + mask + dist + angl)."""

        # pack all the essential data into a dict
        data_dict = {
            'id': core_data['id'],                          # str: 1a1xA00
            'seq': core_data['seq'],                        # str: AGEDVGAPPDHLWVHQEGIYR...
            'token': torch.LongTensor(core_data['token']),  # np.ndarray: (L, N)
            'rna_fm_token': torch.LongTensor(core_data['rna_fm_token']),  # np.ndarray: (L)
            't1ds': torch.FloatTensor(core_data['t1ds'])
            if exists(core_data['t1ds']) else None,
            't2ds': torch.FloatTensor(core_data['t2ds'])
            if exists(core_data['t2ds']) else None,
            'cord_gt': torch.FloatTensor(core_data['cord_gt'])
            if exists(core_data['cord_gt']) else None,      # np.ndarray: randomized 3D coordinates (BS x L/4L x 3)
            'cord_gt_mask': torch.FloatTensor(core_data['cord_gt_mask'])
            if exists(core_data['cord_gt_mask']) else None, # np.ndarray: randomized 3D coordinates (BS x L/4L x 3)
            'cmap': torch.LongTensor(core_data['cmap'])
            if exists(core_data['cmap']) else None,
            'cmap_msk': torch.LongTensor(core_data['cmap_msk'])
            if exists(core_data['cmap_msk']) else None,
            'da_lbls': core_data['da_lbls']
            if exists(core_data['da_lbls']) else None,
            'config': self.config,
        }

        return data_dict

    def __get_crop_data(self, data_dict):

        if self.crop_size > 0 and len(data_dict['seq']) - self.crop_size > 0:
            s = np.random.randint(0, len(data_dict['seq']) - self.crop_size) if self.is_train else 0
            data_dict['seq'] = data_dict['seq'][s:s + self.crop_size]
            data_dict['token'] = data_dict['token'][:, s:s + self.crop_size]
            data_dict['rna_fm_token'] = data_dict['rna_fm_token'][s:s + self.crop_size]
            data_dict['t1ds'] = data_dict['t1ds'][:, s:s + self.crop_size] \
                if exists(data_dict['t1ds']) else None
            data_dict['t2ds'] = data_dict['t2ds'][:, s:s + self.crop_size, s:s + self.crop_size] \
                if exists(data_dict['t2ds']) else None
            data_dict['cmap'] = data_dict['cmap'][s:s + self.crop_size, s:s + self.crop_size] \
                if exists(data_dict['cmap']) else None
            data_dict['cmap_msk'] = data_dict['cmap_msk'][s:s + self.crop_size, s:s + self.crop_size] \
                if exists(data_dict['cmap_msk']) else None
            data_dict['cord_gt'] = data_dict['cord_gt'][:, s * self.n_grps_cord:(s + self.crop_size) * self.n_grps_cord,:] \
                if exists(data_dict['cord_gt']) else None
            data_dict['cord_gt_mask'] = data_dict['cord_gt_mask'][:,s * self.n_grps_cord:(s + self.crop_size) * self.n_grps_cord] \
                if exists(data_dict['cord_gt_mask']) else None

            # for drfold
            # data_dict['dr_msa'] = data_dict['dr_msa'][:, s:s + self.crop_size] \
            #     if exists(data_dict['dr_msa']) else None
            # data_dict['dr_seq'] = data_dict['dr_seq'][s:s + self.crop_size] \
            #     if exists(data_dict['dr_seq']) else None
            # data_dict['dr_hmm'] = data_dict['dr_hmm'][s:s + self.crop_size] \
            #     if exists(data_dict['dr_hmm']) else None
            # data_dict['dr_ss'] = data_dict['dr_ss'][s:s + self.crop_size, s:s + self.crop_size] \
            #     if exists(data_dict['dr_ss']) else None
            # data_dict['dr_seq_fea'] = data_dict['dr_seq_fea'][:, s:s + self.crop_size, ] \
            #     if exists(data_dict['dr_seq_fea']) else None
            # data_dict['dr_pair_fea'] = data_dict['dr_pair_fea'][:, s:s + self.crop_size, s:s + self.crop_size] \
            #     if exists(data_dict['dr_pair_fea']) else None

            if exists(data_dict['da_lbls']):
                for k in data_dict['da_lbls'].keys():
                    data_dict['da_lbls'][k] = data_dict['da_lbls'][k][s:s + self.crop_size, s:s + self.crop_size]

            # if exists(data_dict['dr_da_lbls_soft']):
            #     for k in data_dict['dr_da_lbls_soft'].keys():
            #         data_dict['dr_da_lbls_soft'][k] = data_dict['dr_da_lbls_soft'][k][s:s + self.crop_size, s:s + self.crop_size]

        return data_dict

def unitest_complex():
    msa_depth = 32
    crop_size = 128

    # data_dir = '/public/home/taoshen/data/rna/RNA3D/Ash_Total_RNA_test/'
    # clstr_list = None
    # clstr_id_list = None
    # id_fpath = f'/public/home/taoshen/data/rna/RNA3D/Ash_Total_RNA_test/exp/paper1.0/train_cplx_ids'
    # dr_fea_dpath = f'/public/home/taoshen/data/rMSA_gen/drfold_npzs'
    # msa_dpath = f'{data_dir}/rMSA/'
    # fas_dpath = f'{data_dir}/rna_seq/'
    # pdb_dpath = f'{data_dir}/rna_pdb_align_pad/'
    # dr_ss_dpath = f'{data_dir}/drfold_ss/'
    # dr_cm_dpath = f'{data_dir}/drfold_cm/'
    # danpz_dpath = f'{data_dir}/da_npzs_pad/'
    # ss_dpath =  None

    clstr_list = None
    clstr_id_list = None
    dr_ss_dpath = None
    dr_cm_dpath = None
    ss_dpath =  None
    dr_fea_dpath = None

    id_fpath = f'/public/home/taoshen/code/Projects/E2EFold3D_Dev/data_pipeline/stats/' \
               f'nrlist_3.226_4.0A/cutoff_1.0_all/train_ids_add_cplx'

    data_dir = '/public/home/taoshen/data/rna/RNA3D/mmcif2pdb_rna/'
    msa_dpath = f'{data_dir}/rMSA/'
    fas_dpath = f'{data_dir}/seq/'
    pdb_dpath = f'{data_dir}/pdb/'
    danpz_dpath = f'{data_dir}/da_npzs/'

    pdb_version = 'v0.0'
    frame_version = 'v1.1'

    ds_train_config = DatasetConfig(id_fpath=id_fpath,
                                    clstr_fpath=clstr_list,
                                    clstr_id_fpath=clstr_id_list,
                                    fas_dpath=fas_dpath,
                                    msa_dpath=msa_dpath,
                                    pdb_dpath=pdb_dpath,
                                    crop_size=crop_size,
                                    msa_depth=msa_depth,
                                    pdb_version=pdb_version,
                                    frame_version=frame_version,
                                    ss_dpath=ss_dpath,
                                    dr_fea_dpath=dr_fea_dpath,
                                    dr_ss_dpath=dr_ss_dpath,
                                    dr_cm_dpath=dr_cm_dpath,
                                    is_dr_fea=True,
                                    is_dr_pred=True,
                                    dist_config=default_config.dataset.dist,
                                    danpz_dpath=danpz_dpath,
                                    is_complex=True
                                    )

    ds_train_dataset = RNADataset(ds_train_config, is_train=True)
    for i in tqdm(range(len(ds_train_dataset))):
        data_dict = ds_train_dataset[i]
        # print(data_dict['seq'])
        # for key in data_dict.keys():
        #     if key in set(tns_keys) and exists(data_dict[key]):
        #         print(i, data_dict['id'], key, data_dict[key].shape)

def unitest_list_10cv():

    zfold_init()

    msa_depth = 32  # depth for MSA
    crop_size = 128

    id_fpath = f'/public/home/taoshen/code/Projects/E2EFold3D_Dev/data_pipeline/stats/nrlist_3.226_4.0A/cutoff_1.0_all/train_ids'
    data_dir = '/public/home/taoshen/data/rna/RNA3D/mmcif2pdb_rna/'

    msa_dpath = f'{data_dir}/rMSA/'
    fas_dpath = f'{data_dir}/seq/'  # directory path to SEQ files
    pdb_dpath = f'{data_dir}/pdb/'  # directory path to KEYPDB files
    danpz_dpath = f'{data_dir}/da_npzs/'

    clstr_list = None
    clstr_id_list = None

    pdb_version = 'v0.0'
    frame_version = 'v4.6'

    ds_train_config = DatasetConfig(id_fpath=id_fpath,
                                    clstr_fpath=clstr_list,
                                    clstr_id_fpath=clstr_id_list,
                                    fas_dpath=fas_dpath,
                                    msa_dpath=msa_dpath,
                                    pdb_dpath=pdb_dpath,
                                    crop_size=crop_size,
                                    msa_depth=msa_depth,
                                    pdb_version=pdb_version,
                                    frame_version=frame_version,
                                    dist_config=default_config.dataset.dist,
                                    danpz_dpath=danpz_dpath,
                                    )

    ds_train_dataset = RNADataset(ds_train_config, is_train=True)
    for i in tqdm(range(len(ds_train_dataset))):
        data_dict = ds_train_dataset[i]
        # print(data_dict.keys())

def move_msa():
    data_dir = '/public/home/taoshen/data/rna/RNA3D/mmcif2pdb_rna/'
    msa_dpath = f'{data_dir}/rMSA/'
    fas_dpath = f'{data_dir}/seq/'
    os.makedirs(msa_dpath, exist_ok=True)
    msa_dpath_ori = f'/public/home/taoshen/data/rna/RNA3D/Ash_Total_RNA_test/rMSA/'

    for seq in os.listdir(fas_dpath):
        tgt = seq.replace('.seq', '')
        tgt_ = tgt.replace('_','')
        src = f'{msa_dpath_ori}/{tgt_}.a3m'
        dst = f'{msa_dpath}/{tgt}.a3m'
        if os.path.exists(src):
            print(tgt)
            shutil.move(src, dst)

if __name__ == '__main__':
    unitest_complex()
    # unitest_list_10cv()