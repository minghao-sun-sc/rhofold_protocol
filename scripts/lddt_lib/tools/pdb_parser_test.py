"""Unit-tests for the <PdbParser> class."""

import logging

from rnafold.utils import zfold_init
from rnafold.tools.pdb_parser import PdbParser


def check_outputs(aa_seq, atom_cords, atom_masks, structure, error_msg):
    """Check <PdbParser>'s outputs."""

    if error_msg is not None:
        logging.warning('error message: %s', error_msg)
    else:  # no error is detected
        logging.info('rna sequence: %s (%d residues)', aa_seq, len(aa_seq))
        logging.info('key atom coordinates: %s', str(atom_cords.shape))
        logging.info('atom masks: %s', str(atom_masks.shape))
        logging.info('BioPython structure: %s', str(structure))


def main():
    """Main entry."""

    # configurations
    fas_fpath = '/data/RNA3D_data/key_seq_v2/1KH6.seq'
    pdb_fpath = '/data/RNA3D_data/key_pdb_v2/1KH6.key_pdb'
    # initialization
    zfold_init()
    parser = PdbParser(version='v2.1')

    # parse the PDB file
    logging.info('parsing the PDB file: %s', pdb_fpath)
    aa_seq, atom_cords, atom_masks, structure, error_msg = parser.run(pdb_fpath,fas_fpath=fas_fpath)
    check_outputs(aa_seq, atom_cords, atom_masks, structure, error_msg)


if __name__ == '__main__':
    main()
