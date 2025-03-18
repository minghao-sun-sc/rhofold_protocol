from pathlib import Path

import numpy as np

import argparse

parser = argparse.ArgumentParser(description='Parse pLDDT scores from npz file')
parser.add_argument('--npz_path', type=str, help='Path to npz file')
parser.add_argument('--save_plddt_path', type=str, help='Path to save pLDDT scores')

args = parser.parse_args()

npz_path = args.npz_path  # './results/rhofold/3owz_A/results.npz'
save_plddt_path = args.save_plddt_path  # './results/rhofold/3owz_A/plddt.npy'

distogram_data = np.load(npz_path)
plddt_scores = distogram_data['plddt']  # shape = (1, L)

Path(save_plddt_path).parent.mkidr(parents=True, exist_ok=True)

np.save(save_plddt_path, plddt_scores)

mean_plddt = np.mean(plddt_scores)
print(f'mean pLDDT = {mean_plddt}')
