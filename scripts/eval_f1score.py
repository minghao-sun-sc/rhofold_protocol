import numpy as np

import argparse

def cal_f1score(pred_ss, target_ss):
    """
    Calculate the F1-score for predicted secondary structures compared to a target structure.

    Args:
        pred_ss (np.ndarray): Predicted secondary structure, (L, L).
        target_ss (np.ndarray): Target secondary structure, (L, L).

    Returns:
        float: F1-score.
    """
    pred_pair_num = np.sum(pred_ss)  # total predicted positive
    target_pair_num = np.sum(target_ss)  # total target positive

    tp = np.sum(pred_ss * target_ss)

    precision = tp / pred_pair_num if pred_pair_num > 0 else 0.0
    recall = tp / target_pair_num if target_pair_num > 0 else 0.0

    f1_score = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0

    return f1_score

parser = argparse.ArgumentParser(description='Evaluating the secondary structure F1-score')
parser.add_argument('--pred_ss_path', type=str, help='Path to the predicted 2D structure as .npy contact map')
parser.add_argument('--true_ss_path', type=str, help='Path to the ground truth 2D structure as .npy contact map')

args = parser.parse_args()

pred_ss_path = args.pred_ss_path  # './results/rhofold/3owz_A/ss.npy'
true_ss_path = args.true_ss_path  # './data/rhofold/3owz_A/3owz_A.npy'

pred_ss = np.load(pred_ss_path)
true_ss = np.load(true_ss_path)

f1_score = cal_f1score(pred_ss, true_ss)
print(f'F1-score = {f1_score}')
