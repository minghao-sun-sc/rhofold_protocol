from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO

import argparse

parser = argparse.ArgumentParser(description='Add pLDDT scores as B-factors to a PDB file')
parser.add_argument('--plddt_path', type=str, help='Path to pLDDT scores')
parser.add_argument('--pdb_path', type=str, help='Path to PDB file')
args = parser.parse_args()

plddt_path = Path(args.plddt_path)
pdb_path = Path(args.pdb_path)

# plddt_path = Path('./results/rhofold/3owz_A/plddt.npy')
# pdb_path = Path('./results/rhofold/3owz_A/unrelaxed_model.pdb')

# Load and flatten the pLDDT scores
plddt_scores = np.load(plddt_path).flatten()

# Parse the structure from the PDB file
parser = PDBParser(QUIET=True)
structure = parser.get_structure('RNA', str(pdb_path))

# Update B-factors for each residue using corresponding pLDDT scores
for i, residue in enumerate(structure.get_residues()):
    if i < len(plddt_scores):
        score = plddt_scores[i]
    else:
        print(f'Warning: pLDDT score not available for residue {residue.get_id()}')
        score = 0.0
    for atom in residue:
        atom.set_bfactor(score)

# Save the updated structure to a new PDB file
output_path = pdb_path.parent / f'{pdb_path.stem}_plddt.pdb'
io = PDBIO()
io.set_structure(structure)
io.save(str(output_path))

print(f'Updated PDB file saved as {output_path}')
