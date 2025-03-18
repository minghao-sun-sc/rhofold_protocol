import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import argparse

parser = argparse.ArgumentParser(description='Visualize secondary structure from a CT file')
parser.add_argument('--ct_path', type=str, help='Path to CT file')
parser.add_argument('--save_db_path', type=str, help='Path to save the dot-bracket notation')
parser.add_argument('--save_npy_path', type=str, help='Path to save the .npy contact map')
parser.add_argument('--save_ss_plot_path', type=str, help='Path to save the contact map plot')
args = parser.parse_args()

ct_path = args.ct_path
db_path = args.save_db_path
npy_path = args.save_npy_path
save_ss_plot_path = args.save_ss_plot_path

# ct_path = './results/rhofold/3owz_A/ss.ct'
# db_path = './results/rhofold/3owz_A/ss.db'  # path to save the dot-bracket notation
# npy_path = './results/rhofold/3owz_A/ss.npy'  # path to save the .npy contact map
# save_ss_plot_path = './results/rhofold/3owz_A/ss.png'  # path to save the contact map plot

ct_df = pd.read_csv(ct_path, sep='\t', header=None, skiprows=[0],
                    usecols=[1, 4])  # skip the first row which is header line
# column 2 is the nucleotide, and column 5 is the one-based index for the paired nucleotide

seq = ''.join(ct_df[1].values)

# convert to dot-bracket format
db_list = ['.' for i in range(len(seq))]  # initialize all as unpaired
for i in range(len(seq)):
    pair_idx = ct_df.iloc[i, 1] - 1

    if pair_idx != 0 and pair_idx > i:
        db_list[i] = '('
        db_list[pair_idx - 1] = ')'

db = ''.join(db_list)

# save the dot-bracket notation to a file
with open(db_path, 'w') as f:
    f.write(f'>{seq}\n{db}')

# convert to contact map
contact_map = np.zeros((len(seq), len(seq)))

for i in range(len(seq)):
    pair_idx = ct_df.iloc[i, 1] - 1

    if pair_idx != 0 and pair_idx > i:
        contact_map[i, pair_idx] = 1
        contact_map[pair_idx, i] = 1

np.save(npy_path, contact_map)

plt.figure(figsize=(10, 8))

# heat map
sns.heatmap(contact_map, cmap='Purples')

# boundary box
plt.gca().add_patch(
    plt.Rectangle((0, 0), contact_map.shape[1], contact_map.shape[0], fill=False, edgecolor='black', lw=3))

# move the x ticks to the top
plt.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False)

plt.savefig(save_ss_plot_path, dpi=300, bbox_inches='tight')
