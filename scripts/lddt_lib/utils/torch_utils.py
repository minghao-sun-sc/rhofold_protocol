"""PyTorch-related utility functions."""

import logging

import torch


def get_tensor_size(tensor):
    """Get the PyTorch tensor's memory consumption (in MB).

    Args:
    * tensor: PyTorch tensor

    Returns:
    * mem: PyTorch tensor's memory consumption (in MB)
    """

    mem = tensor.element_size() * torch.numel(tensor) / 1024.0 / 1024.0

    return mem


def check_tensor_size(tensor, name, mem_thres=500.0):
    """Check whether the PyTorch tensor consumes more memory than the threshold.

    Args:
    * tensor: PyTorch tensor
    * name: PyTorch tensor's name
    * mem_thres: memory consumption's threshold (in MB)

    Returns: n/a
    """

    mem = get_tensor_size(tensor)
    if mem > mem_thres:
        logging.debug('tensor <%s> has consumed %.2f MB memory', name, mem)


def get_peak_memory():
    """Get the peak memory consumption (in GB).

    Args: n/a

    Returns:
    * mem: peak memory consumption (in GB)
    """

    mem = torch.cuda.memory_stats()['allocated_bytes.all.peak'] / 1024.0 / 1024.0

    return mem


def load_pretrain(model,
                  pretrain_file,
                  is_GPU=False,
                  is_fair=False,
                  ignore_mismatch=False,
                  skpi_key = None
                  ):
    print('    loading', pretrain_file)

    pretrain_state_dict = torch.load(pretrain_file) if is_GPU else torch.load(pretrain_file, map_location='cpu')

    if is_fair:
        pretrain_state_dict = pretrain_state_dict['model']

    state_dict = model.state_dict()
    for key in state_dict.keys():

        if skpi_key is not None and skpi_key in key:
            print(f'skip {key}, {skpi_key}')
            continue

        shape = state_dict[key].shape

        if key in pretrain_state_dict.keys():
            param = pretrain_state_dict[key]
        elif ('encoder.rnafold.' + key) in pretrain_state_dict.keys():
            param = pretrain_state_dict['encoder.rnafold.' + key]
        else:
            param = None
            print(f'{key} not in pretrain_state_dicts')

        if param is not None:
            if ignore_mismatch:
                if param.shape == shape:
                    state_dict[key] = param
                else:
                    print(f'skip {key}, size mismatch {param.shape} {shape}')
            else:
                state_dict[key] = param

    model.load_state_dict(state_dict)

    return model