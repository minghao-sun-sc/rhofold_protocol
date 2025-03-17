import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(description='Visualize secondary structure probabilities')
parser.add_argument('--npz_path', type=str, help='Path to NPZ file with secondary structure probabilities')
parser.add_argument('--save_ss_plot_path', type=str, help='Path to save the secondary structure plot')
args = parser.parse_args()

npz_path = args.npz_path
save_ss_plot_path = args.save_ss_plot_path

# npz_path = './results/rhofold/3owz_A/results.npz'
# save_ss_plot_path = './results/rhofold/3owz_A/ss_prob_map.png'

distogram_data = np.load(npz_path)

ss_prob_map = distogram_data['ss_prob_map']

# visualize the probability map
plt.figure(figsize=(10, 8))

# draw the heatmap
sns.heatmap(ss_prob_map, cmap='Purples')

# boundary box
plt.gca().add_patch(plt.Rectangle((0, 0), ss_prob_map.shape[1], ss_prob_map.shape[0], fill=False, edgecolor='black', lw=3))

# move the x ticks to the top
plt.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False)

plt.savefig(save_ss_plot_path, dpi=300, bbox_inches='tight')
