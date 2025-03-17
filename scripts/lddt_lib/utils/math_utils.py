"""Mathematics-related utility functions."""

import numpy as np
import torch

def cdist(x1, x2=None):
    """Calculate the pairwise distance matrix.

    Args:
    * x1: input tensor of size N x D or B x N x D
    * x2: (optional) input tensor of size M x D or B x M x D

    Returns:
    * dist_mat: pairwise distance of size N x M or B x N x M

    Note:
    * If <x2> is not provided, then pairwise distance will be computed within <x1>.
    * The matrix multiplication approach will not be used to avoid the numerical stability issue.
    """

    # initialization
    compute_mode = 'donot_use_mm_for_euclid_dist'
    x2 = x1 if x2 is None else x2

    # calculate the pairwise distance matrix
    if x1.ndim == 2 and x2.ndim == 2:
        dist_mat = torch.cdist(
            x1.unsqueeze(dim=0), x2.unsqueeze(dim=0), compute_mode=compute_mode).squeeze(dim=0)
    elif x1.ndim == 3 and x2.ndim == 3:
        dist_mat = torch.cdist(x1, x2, compute_mode=compute_mode)
    else:
        raise ValueError('<x1> and <x2> must be either in the 2- or 3-dimension')

    '''
    # calculate the pairwise distance matrix
    if x1.ndim == 2 and x2.ndim == 2:
        dist_mat = torch.norm(x1.unsqueeze(dim=1) - x2.unsqueeze(dim=0), dim=2)
    elif x1.ndim == 3 and x2.ndim == 3:
        dist_mat = torch.norm(x1.unsqueeze(dim=2) - x2.unsqueeze(dim=1), dim=3)
    else:
        raise ValueError('<x1> and <x2> must be either in the 2- or 3-dimension')
    '''

    return dist_mat


def cvt_to_one_hot(arr, depth):
    """Convert an integer array into one-hot encodings.

    Args:
    * arr: integer array of size (D1, D2, ..., Dk)
    * depth: one-hot encodings's depth (denoted as C)

    Returns:
    * arr_oht: one-hot encodings of size (D1, D2, ..., Dk, C)
    """

    assert np.min(arr) >= 0 and np.max(arr) < depth
    arr_oht = np.reshape(np.eye(depth)[arr.ravel()], list(arr.shape) + [depth])

    return arr_oht


def get_rotate_mat():
    """Get a randomized 3D rotation matrix.

    Args: n/a

    Returns:
    * rotate_mat: 3D rotation matrix
    """

    # generate a randomized 3D rotation matrix
    yaw = np.random.uniform(-np.pi, np.pi)
    pitch = np.random.uniform(-np.pi, np.pi)
    roll = np.random.uniform(-np.pi, np.pi)
    sy, cy = np.sin(yaw), np.cos(yaw)
    sp, cp = np.sin(pitch), np.cos(pitch)
    sr, cr = np.sin(roll), np.cos(roll)
    rotate_mat = np.array([
        [cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
        [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
        [-sp, cp * sr, cp * cr],
    ], dtype=np.float32)

    return rotate_mat


def calc_plane_angle(cord_1, cord_2, cord_3):
    """Calculate the plane angle defined by 3 points' 3-D coordinates.

    Args:
    * cord_1: 3-D coordinate of the 1st point
    * cord_2: 3-D coordinate of the 2nd point
    * cord_3: 3-D coordinate of the 3rd point

    Returns:
    * rad: planar angle (in radian)
    """

    eps = 1e-6
    a1 = cord_1 - cord_2
    a2 = cord_3 - cord_2
    rad = np.arccos(np.clip(
        np.dot(a1, a2) / (np.linalg.norm(a1) * np.linalg.norm(a2) + eps), -1.0, 1.0))

    return rad


def calc_dihedral_angle(cord_1, cord_2, cord_3, cord_4):
    """Calculate the dihedral angle defined by 4 points' 3-D coordinates.

    Args:
    * cord_1: 3-D coordinate of the 1st point
    * cord_2: 3-D coordinate of the 2nd point
    * cord_3: 3-D coordinate of the 3rd point
    * cord_4: 3-D coordinate of the 4th point

    Returns:
    * rad: dihedral angle (in radian)
    """

    eps = 1e-6
    a1 = cord_2 - cord_1
    a2 = cord_3 - cord_2
    a3 = cord_4 - cord_3
    v1 = np.cross(a1, a2)
    v1 = v1 / np.sqrt((v1 * v1).sum(-1) + eps)
    v2 = np.cross(a2, a3)
    v2 = v2 / np.sqrt((v2 * v2).sum(-1) + eps)
    sign = np.sign((v1 * a3).sum(-1))
    rad = np.arccos(np.clip(
        (v1 * v2).sum(-1) / np.sqrt((v1 ** 2).sum(-1) * (v2 ** 2).sum(-1) + eps), -1.0, 1.0))
    if sign != 0:
        rad *= sign

    return rad


def calc_denc_tns(dist_tns, base=2.0, n_dims=11):
    """Calculate the distance encoding tensor with log-spaced thresholds.

    Args:
    * dist_tns: distance tensor of size N1 x N2 x ... x Nd
    * base: distance threshold's base
    * n_dims: number of dimensions for distance encoding

    Returns:
    * denc_tns: distance encoding tensor of size N1 x N2 x ... x Nd x D

    Note:
    Both <dist_tns> and <denc_tns> are PyTorch tensors.
    """

    dist_vals_np = np.reshape(np.power(base, np.arange(n_dims)), [1] * dist_tns.ndim + [n_dims])
    dist_vals = torch.tensor(dist_vals_np, dtype=torch.float32, device=dist_tns.device)
    denc_tns = torch.sigmoid(torch.unsqueeze(dist_tns, dim=-1) / dist_vals - 1.0)

    return denc_tns
