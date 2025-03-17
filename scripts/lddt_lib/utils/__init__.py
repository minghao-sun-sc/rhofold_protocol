"""Import all the utility functions (and constants)."""

from .comm_utils import zfold_init
from .comm_utils import get_md5sum
from .comm_utils import get_rand_str
from .comm_utils import get_nb_threads
from .comm_utils import make_config_list
from .file_utils import get_tmp_dpath
from .file_utils import clear_tmp_files
from .file_utils import find_files_by_suffix
from .file_utils import recreate_directory
from .file_utils import unpack_archive
from .file_utils import make_archive
from .math_utils import cvt_to_one_hot
from .math_utils import get_rotate_mat
from .math_utils import calc_plane_angle
from .math_utils import calc_dihedral_angle
from .math_utils import calc_denc_tns
from .torch_utils import get_tensor_size
from .torch_utils import check_tensor_size
from .torch_utils import get_peak_memory
from .torch_utils import load_pretrain
from .struct3d_utils import kabsch_numpy
from .struct3d_utils import kabsch_torch

def get_parameter_number(net):
    total_num = sum(p.numel() for p in net.parameters())
    trainable_num = sum(p.numel() for p in net.parameters() if p.requires_grad)
    return {'Total': total_num, 'Trainable': trainable_num}

__all__ = [
    'zfold_init',
    'get_md5sum',
    'get_rand_str',
    'get_nb_threads',
    'make_config_list',
    'get_tmp_dpath',
    'clear_tmp_files',
    'find_files_by_suffix',
    'recreate_directory',
    'unpack_archive',
    'make_archive',
    'cvt_to_one_hot',
    'get_rotate_mat',
    'calc_plane_angle',
    'calc_dihedral_angle',
    'calc_denc_tns',
    'get_tensor_size',
    'check_tensor_size',
    'get_peak_memory',
    'kabsch_numpy',
    'kabsch_torch',
    'load_pretrain'
]
