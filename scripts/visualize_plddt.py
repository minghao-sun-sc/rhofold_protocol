from pymol import cmd

import argparse

parser = argparse.ArgumentParser(description='Visualize pLDDT scores on a protein structure')
parser.add_argument('--pdb_path', type=str, help='Path to PDB file with pLDDT scores as B-factors')
args = parser.parse_args()

pdb_path = args.pdb_path  # './results/rhofold/3owz_A/unrelaxed_model_plddt.pdb'

# Load the PDB file
cmd.load(pdb_path, 'unrelaxed_model_plddt')

# Color the model based on the B-factor using the 'blue_white_red' spectrum,
# with the color range specified from 0 to 1.
cmd.spectrum('b', 'blue_white_red', 'unrelaxed_model_plddt', minimum=0, maximum=1)

# Automatically orient the structure: centers the object and applies a rotation based on its mass distribution
cmd.orient('unrelaxed_model_plddt')

# Optionally, set the background color to white for better contrast
cmd.bg_color('white')

# Render the scene using ray tracing for a high-quality image and save as PNG.
# You can adjust the width, height, and dpi as needed.
cmd.png('unrelaxed_model_plddt.png', width=1200, height=800, dpi=300, ray=1)
