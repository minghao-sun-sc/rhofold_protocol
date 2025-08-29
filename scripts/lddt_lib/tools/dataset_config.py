import logging

class DatasetConfig():
    """Configurations for the <Dataset> class."""

    def __init__(
            self,
            id_fpath,            # IDs file
            fas_dpath,           # directory path to FASTA files
            msa_dpath,           # directory path to MSA files
            is_complex=False,
            clstr_id_fpath=None, # IDs file
            clstr_fpath=None,    # IDs file
            pdb_dpath=None,      # directory path to PDB files
            danpz_dpath=None,    # directory path to PDB files
            ss_dpath=None,       # directory path to ss npy files
            ff_dpath=None,
            tpl_dpath=None,
            dr_dpath=None,
            dr_fea_dpath=None,
            dr_ss_dpath=None,    # directory path to ss npy files
            dr_cm_dpath=None,    # directory path to ss npy files
            is_dr_fea=False,
            is_dr_pred=False,
            num=None,
            pdb_version='v0.0',
            frame_version='v4.3',
            crop_size = 64,      # crop size
            msa_depth = 64,      # msa_depth
            dist_config = None,
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
        self.noise_std_max = 10.0
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

        if self.id_fpath is not None:
            logging.info('id_fpath: %s', self.id_fpath)
            logging.info('rna num: %d', len(self.rna_ids))

        if self.clstr_id_fpath is not None and self.clstr_fpath is not None:
            logging.info('clstr_id_fpath: %s', self.clstr_id_fpath)
            logging.info('clstr_id_fpath: %s', self.clstr_fpath)
            logging.info('clstr num: %d', len(self.clstr_ids))

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

def parse_clstr_ids(file = '/public/home/taoshen/data/rna/RNA3D/data_px/cdhit_out_0.8_val_clstr.list'):

    f = open(file, 'r')
    lines = f.readlines()
    f.close()
    clstr_ids = [l.strip().split(', ')[0][2:-1] for l in lines]

    return clstr_ids


def parse_clstr(file = '/public/home/taoshen/data/rna/RNA3D/data_px/cdhit_out_0.8.clstr'):

    f = open(file, 'r')
    lines = f.readlines()
    f.close()

    items = []
    items_tmp = [lines[0]]
    for l in lines[1:]:
        if '>Cluster' in l:
            items.append(items_tmp)
            items_tmp = [l]
        else:
            items_tmp.append(l)

    items_dict = {}
    for item in items:
        clstr_id = item[0].strip()[1:]
        rna_ids = []
        for i in range(len(item)-1):
            inx = item[i+1].find('>')
            rna_id = item[i+1][inx + 1:inx + 7]
            rna_ids.append(rna_id)
        items_dict[clstr_id] = rna_ids

    return items_dict

if __name__ == '__main__':
    clstrs = parse_clstr()
    clstr_ids = parse_clstr_ids()
    for clstr_id in clstr_ids:
        print(clstr_id, clstrs[clstr_id])

