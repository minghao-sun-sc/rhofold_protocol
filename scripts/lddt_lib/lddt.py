"""The built-in assessor for per-residue lDDT scores."""

import os

import torch

from .utils.math_utils import cdist
from .utils.metric_utils import KabschRMSD
from .tools.pdb_parser import PdbParser

from tqdm import tqdm
import numpy as np

class LddtAssessor():
    """The built-in assessor for lDDT scores.

    Note:
    * For simplicity, we do not consider the naming issue of partially symmetric residues.
    """

    def __init__(self):
        """Constructor function."""

        self.eps = 1e-6
        self.radius = 15.0  # inclusion radius
        self.dthr_tns = torch.tensor([0.5, 1.0, 2.0, 4.0], dtype=torch.float32).view(1, 1, -1)


    def run(self, cord_tns_ref, cord_tns_qry, cmsk_mat, atom_set='c1_'):
        """Run the built-in assessor to compute per-residue lDDT-Ca scores.

        Args:
        * cord_tns_ref: reference structure's per-atom 3D coordinates of size L x M x 3
        * cord_tns_qry: query structure's per-atom 3D coordinates of size L x M x 3
        * cmsk_mat: reference structure's per-atom 3D coordinates' validness masks of size L x M
        * atom_set: (optional) atom set over which to compute per-residue lDDT scores

        Returns:
        * plddt_vec: per-residue lDDT scores of size L
        * plmsk_vec: per-residue lDDT scores' validness masks of size L
        * clddt_val: full-chain lDDT score of size 1
        """

        # initializations
        device = cord_tns_ref.device
        n_resds, n_atoms, _ = cord_tns_ref.shape

        aa_seq = 'A' * n_resds  # does not matter

        # move <self.dthr_tns> to the correct device
        if self.dthr_tns.device != device:
            self.dthr_tns = self.dthr_tns.to(device)

        # extract 3D coordinates for CA-atom or full-atom
        if atom_set == 'c1_':
            n_atoms = 1
            cord_mat_ref = cord_tns_ref[:, 1, :] # L x 3
            cord_mat_qry = cord_tns_qry[:, 1, :] # L x 3
            cmsk_vec = cmsk_mat[:, 1]
            plmsk_vec = cmsk_mat[:, 1]
        elif atom_set == 'all':
            cord_mat_ref = cord_tns_ref.view(n_resds * n_atoms, 3)  # (L x M) x 3
            cord_mat_qry = cord_tns_qry.view(n_resds * n_atoms, 3)  # (L x M) x 3
            cmsk_vec = cmsk_mat.view(n_resds * n_atoms)
            plmsk_vec = torch.max(cmsk_mat, dim=1)[0]  # L
        else:
            raise ValueError('unrecognized atom set: ' + atom_set)

        # calculate pairwise distance matrices
        dist_mat_ref = cdist(cord_mat_ref)  # L x L or (L x M) x (L x M)
        dist_mat_qry = cdist(cord_mat_qry)
        rmsk_mat = 1 - torch.eye(n_resds, dtype=torch.int8, device=device)  # 1: same residue

        if atom_set == 'all' or atom_set == 'bb':
            rmsk_mat = rmsk_mat.repeat_interleave(n_atoms, dim=0).repeat_interleave(n_atoms, dim=1)

        dmsk_mat = rmsk_mat * \
            torch.outer(cmsk_vec, cmsk_vec) * (dist_mat_ref <= self.radius).to(torch.int8)

        # calculate the hit ratio under each distance threshold
        derr_mat = torch.abs(dist_mat_qry - dist_mat_ref)
        dhit_mat = torch.mean((derr_mat.unsqueeze(dim=2) <= self.dthr_tns).to(torch.float32), dim=2)

        if atom_set == 'all' or atom_set == 'bb':
            plddt_vec = torch.sum(dmsk_mat.view(n_resds, -1) * dhit_mat.view(n_resds, -1), dim=1) \
                        / (torch.sum(dmsk_mat.view(n_resds, -1), dim=1) + self.eps)
        else:
            plddt_vec = torch.sum(dmsk_mat * dhit_mat, dim=1) / (torch.sum(dmsk_mat, dim=1) + self.eps)

        # calculate the full-chain lDDT score
        clddt_val = torch.sum(plmsk_vec * plddt_vec) / (torch.sum(plmsk_vec) + self.eps)

        return plddt_vec, plmsk_vec, clddt_val

def main():
    """Main entry."""

    #
    pdb_version = 'v0.0'

    frame_version = 'v4.6'

    rmsd = KabschRMSD()

    parser = PdbParser(pdb_version=pdb_version, frame_version=frame_version, check_mode='lenient')
    
    tgt_all = ['R1107','R1108','R1116','R1117','R1149','R1156']
    gt_path = '/Users/ashhu/Downloads/DSSR/casp_ntv'
    af3_path = '/Users/ashhu/Downloads/DSSR/casp_af3'
    fast_path = '/Users/ashhu/Downloads/E2EfoldFM-2023-10-18/CASP15/CASP15_RNA/fasta/'

    for tgt in tgt_all:

        fasta = os.path.join(fast_path,tgt+'.fasta')
        pdb_fpath_natv =  os.path.join(gt_path,tgt+'.pdb')
        pdb_fpath_decy =   os.path.join(af3_path,tgt+'_af3.pdb')

        if not os.path.exists(fasta):
            continue

        try:
            aa_seq_ref, cord_tns_ref, cmsk_mat_ref, structure_ref, error_msg1 = \
                parser.run(pdb_fpath_natv, fas_fpath=fasta)
            aa_seq_dcy, cord_tns_qry, cmsk_mat_qry, structure_dcy, error_msg2 = \
                parser.run(pdb_fpath_decy, fas_fpath=fasta)
            cmsk_mat = cmsk_mat_ref * cmsk_mat_qry
            # # test w/ <LddtAssessor>
            assessor = LddtAssessor()
            # reconstruct CA-atom 3D coordinates
            plddt_vec, plmsk_vec, clddt_val = assessor.run(torch.FloatTensor(cord_tns_ref),
                                                           torch.FloatTensor(cord_tns_qry),
                                                           torch.LongTensor(cmsk_mat))
                                                           #atom_set='c1_')
            print(f'{tgt} ' + 'lDDT (C1-only): %.4f' % clddt_val.item())
        except:
            continue


if __name__ == '__main__':
    main()
