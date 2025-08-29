import logging
import os
import sys

import numpy as np
import torch
from Bio import SeqIO
from rhofold.config import rhofold_config
from rhofold.model.rna_fm.data import get_fm_token_from_seq
from rhofold.rhofold import RhoFold
from rhofold.utils import get_device, timing
from rhofold.utils.alphabet import get_features_for_fm
from tqdm import tqdm


@torch.no_grad()
def main(config):
    """
    RNA-FM embedding generation pipeline
    """

    os.makedirs(config.output_dir, exist_ok=True)

    logger = logging.getLogger('RhoFold Inference')
    logger.setLevel(level=logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s')
    file_handler = logging.FileHandler(f'{config.output_dir}/log.txt', mode='w')
    file_handler.setLevel(level=logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    logger.info(f'Constructing RhoFold')
    model = RhoFold(rhofold_config)

    logger.info(f'Loading {config.ckpt}')
    model.load_state_dict(torch.load(config.ckpt, map_location=torch.device('cpu'))['model'])
    model.eval()

    # Input seq, MSA
    logger.info(f'Input FASTA path: {config.input_fasta}')

    with timing('RNA-FM embedding generation', logger=logger):
        config.device = get_device(config.device)
        logger.info(f'    Generating RNA-FM embedding using device {config.device}')
        model = model.to(config.device)

        records = list(SeqIO.parse(config.input_fasta, 'fasta'))

        for record in tqdm(records):
            seq_dict = {record.id: str(record.seq)}
            fm_tokens = get_fm_token_from_seq(seq_dict)

            # Use the RNA-FM tokens to generate the embeddings
            outputs = model.msa_embedder.rna_fm(fm_tokens.to(config.device),
                                                need_head_weights=False,
                                                repr_layers=[12],
                                                return_contacts=False)
            # (1, L)

            # save the embeddings
            embeddings = outputs['representations'][12].squeeze(0).cpu().numpy()
            # (1, L, 640) => (L, 640)

            output_path = os.path.join(config.output_dir, f'{record.id}.npy')

            np.save(output_path, embeddings)

        embeddings = outputs['representations'][12]
        print(embeddings.shape)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--device',
                        help='Default cpu. If GPUs are available, you can set --device cuda:<GPU_index>.',
                        default='cpu')
    parser.add_argument('--ckpt',
                        help='Path to the pretrained model, default ./checkpoints/rhofold_pretrained_params.pt.',
                        default='./checkpoints/rhofold_pretrained_params.pt')
    parser.add_argument('--input_fasta',
                        help='Path to the input fasta file. Valid nucleic acids in RNA sequence: A, U, G, C.',
                        required=True)
    parser.add_argument('--output_dir',
                        help='Path to the output dir.',
                        required=True)

    args = parser.parse_args()

    main(args)
