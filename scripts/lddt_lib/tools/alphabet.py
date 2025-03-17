from typing import Sequence, Tuple, List, Union
import re

import numpy as np
import torch
import pickle
import random
import itertools
from Bio import SeqIO
from typing import List, Tuple
import string

# rnaseq_toks = {
#     'toks': ['A','U','G','C','-']
# }

rnaseq_toks = {
    'toks': ['A','U','G','C','-','T']
}

RawMSA = Sequence[Tuple[str, str]]

deletekeys = dict.fromkeys(string.ascii_lowercase)
deletekeys["."] = None
deletekeys["*"] = None
deletekeys[" "] = None # add remove
translation = str.maketrans(deletekeys)

def read_fasta(
    path, keep_gaps=True, keep_insertions=True, to_upper=False,
):
    with open(path, "r") as f:
        for result in read_alignment_lines(
            f, keep_gaps=keep_gaps, keep_insertions=keep_insertions, to_upper=to_upper
        ):
            yield result

def read_alignment_lines(
    lines, keep_gaps=True, keep_insertions=True, to_upper=False,
):
    seq = desc = None

    def parse(s):
        if not keep_gaps:
            s = re.sub("-", "", s)
        if not keep_insertions:
            s = re.sub("[a-z]", "", s)
        return s.upper() if to_upper else s

    for line in lines:
        # Line may be empty if seq % file_line_width == 0
        if len(line) > 0 and line[0] == ">":
            if seq is not None:
                yield desc, parse(seq)
            desc = line.strip()
            seq = ""
        else:
            assert isinstance(seq, str)
            seq += line.strip()
    assert isinstance(seq, str) and isinstance(desc, str)
    yield desc, parse(seq)

class Alphabet(object):

    def __init__(
        self,
        standard_toks: Sequence[str],
        prepend_toks: Sequence[str] = ("<null_0>", "<pad>", "<eos>", "<unk>"),
        append_toks: Sequence[str] = ("<cls>", "<mask>", "<sep>"),
        prepend_bos: bool = True,
        append_eos: bool = False,
        use_msa: bool = False,
    ):

        self.standard_toks = list(standard_toks)
        self.prepend_toks = list(prepend_toks)
        self.append_toks = list(append_toks)
        self.prepend_bos = prepend_bos
        self.append_eos = append_eos
        self.use_msa = use_msa

        self.all_toks = list(self.prepend_toks)
        self.all_toks.extend(self.standard_toks)

        for i in range((8 - (len(self.all_toks) % 8)) % 8):
            self.all_toks.append(f"<null_{i  + 1}>")

        self.all_toks.extend(self.append_toks)

        self.tok_to_idx = {tok: i for i, tok in enumerate(self.all_toks)}

        self.unk_idx = self.tok_to_idx["<unk>"]
        self.padding_idx = self.get_idx("<pad>")
        self.cls_idx = self.get_idx("<cls>")
        self.mask_idx = self.get_idx("<mask>")
        self.eos_idx = self.get_idx("<eos>")
        self.sep_idx = self.get_idx("<sep>")

        self.nspecial = set([self.unk_idx, self.padding_idx, self.cls_idx, self.eos_idx, self.sep_idx])

    def __len__(self):
        return len(self.all_toks)

    def pad(self):
        return self.padding_idx

    def cls(self):
        return self.cls_idx

    def eos(self):
        return self.eos_idx

    def mask(self):
        return self.mask_idx

    def get_idx(self, tok):
        return self.tok_to_idx.get(tok, self.unk_idx)

    def get_idx_fasta(self, toks):
        return [self.get_idx(s) for s in toks]

    def get_idx_msa(self, toks_msa):
        assert isinstance(toks_msa, list)
        return [self.get_idx_fasta(toks) for toks in toks_msa]

    def get_tok(self, ind):
        return self.all_toks[ind]

    def to_dict(self):
        return {"toks": self.toks}

    def get_batch_converter(self):
        if self.use_msa:
            return MSABatchConverter(self)
        else:
            return BatchConverter(self)

    @classmethod
    def from_dict(cls, d, **kwargs):
        return cls(standard_toks=d["toks"], **kwargs)

    @classmethod
    def from_architecture(cls, name: str, ) -> "Alphabet":
        if name in ("RNA MSA Transformer", "rna_msa_transformer", "RNA"):
            standard_toks = rnaseq_toks["toks"]
            prepend_toks = ("<cls>", "<pad>", "<eos>", "<unk>")
            append_toks = ("<mask>",)
            prepend_bos = True
            append_eos = False
            use_msa = True
        else:
            raise ValueError("Unknown architecture selected")
        return cls(
            standard_toks, prepend_toks, append_toks, prepend_bos, append_eos, use_msa
        )

class BatchConverter(object):
    """Callable to convert an unprocessed (labels + strings) batch to a
    processed (labels + tensor) batch.
    """

    def __init__(self, alphabet):
        self.alphabet = alphabet

    def __call__(self, raw_batch: Sequence[Tuple[str, str]]):
        # RoBERTa uses an eos token, while ESM-1 does not.
        batch_size = len(raw_batch)

        max_len = max(len(seq_str) for _, seq_str in raw_batch)

        tokens = torch.empty(
            (
                batch_size,
                max_len
                + int(self.alphabet.prepend_bos)
                + int(self.alphabet.append_eos),
            ),
            dtype=torch.int64,
        )
        tokens.fill_(self.alphabet.padding_idx)
        labels = []
        strs = []

        for i, (label, seq_str) in enumerate(raw_batch):
            labels.append(label)
            strs.append(seq_str)
            if self.alphabet.prepend_bos:
                tokens[i, 0] = self.alphabet.cls_idx
            seq = torch.tensor(
                [self.alphabet.get_idx(s) for s in seq_str], dtype=torch.int64
            )
            tokens[
                i,
                int(self.alphabet.prepend_bos) : len(seq_str)
                + int(self.alphabet.prepend_bos),
            ] = seq
            if self.alphabet.append_eos:
                tokens[
                    i, len(seq_str) + int(self.alphabet.prepend_bos)
                ] = self.alphabet.eos_idx

        return labels, strs, tokens


class MSABatchConverter(BatchConverter):

    def __call__(self, inputs: Union[Sequence[RawMSA], RawMSA]):
        if isinstance(inputs[0][0], str):
            # Input is a single MSA
            raw_batch: Sequence[RawMSA] = [inputs]  # type: ignore
        else:
            raw_batch = inputs  # type: ignore

        batch_size = len(raw_batch)
        max_alignments = max(len(msa) for msa in raw_batch)
        max_seqlen = max(len(msa[0][1]) for msa in raw_batch)

        tokens = torch.empty(
            (
                batch_size,
                max_alignments,
                max_seqlen
                + int(self.alphabet.prepend_bos)
                + int(self.alphabet.append_eos),
            ),
            dtype=torch.int64,
        )
        tokens.fill_(self.alphabet.padding_idx)
        labels = []
        strs = []

        for i, msa in enumerate(raw_batch):
            msa_seqlens = set(len(seq) for _, seq in msa)
            if not len(msa_seqlens) == 1:
                raise RuntimeError(
                    "Received unaligned sequences for input to MSA, all sequence "
                    "lengths must be equal."
                )
            msa_labels, msa_strs, msa_tokens = super().__call__(msa)
            labels.append(msa_labels)
            strs.append(msa_strs)
            tokens[i, :msa_tokens.size(0), :msa_tokens.size(1)] = msa_tokens

        return labels, strs, tokens

def read_sequence(filename: str) -> Tuple[str, str]:
    """ Reads the first (reference) sequences from a fasta or MSA file."""
    record = next(SeqIO.parse(filename, "fasta"))
    return record.description, str(record.seq)

def remove_insertions(sequence: str) -> str:
    """ Removes any insertions into the sequence. Needed to load aligned sequences in an MSA. """
    return sequence.translate(translation)

def read_msa(filename: str, nseq: int, mode: string) -> List[Tuple[str, str]]:
    """ Reads the first nseq sequences from an MSA file, automatically removes insertions."""

    if mode == 'TopN':
        return [(record.description, remove_insertions(str(record.seq))) for record in itertools.islice(SeqIO.parse(filename, "fasta"), nseq)]

    elif mode == 'Rand':
        max_num = 10000
        all_seq = [(record.description, remove_insertions(str(record.seq))) for record in itertools.islice(SeqIO.parse(filename, "fasta"), max_num)]
        origin_seq = all_seq[:1]
        align_seq = all_seq[1: ]
        random.shuffle(align_seq)
        return origin_seq + align_seq[:nseq-1]

    else:
        raise NotImplementedError(f'Mode not exists {mode}')


def read_fas(filename: str):
    """ Reads the first nseq sequences from an MSA file, automatically removes insertions."""
    return [(record.description, remove_insertions(str(record.seq))) for record in itertools.islice(SeqIO.parse(filename, "fasta"), 1)]

def _get_msa_feature(msa_data,
                     batch_converter = Alphabet.from_architecture('RNA MSA Transformer').get_batch_converter(),
                     is_rm = True):

    if is_rm:
        # TODO
        # dirty code for removing '-' in protein seq
        rna_seq = msa_data[0][0][1]
        idx = set([i for i in range(len(rna_seq)) if rna_seq[i] == '-'])

        msa_ = []
        for id, seq in msa_data[0]:
            new_seq = ''
            for i in range(len(rna_seq)):
                if i not in idx:
                    new_seq += seq[i]
            msa_.append([id, new_seq])
        msa_data[0] = msa_

    _, _, msa_batch_tokens = batch_converter(msa_data)
    # remove [cls] token in msa_batch_tokens
    fea1d = msa_batch_tokens.squeeze(0).data.cpu().numpy().transpose((1, 0))[1:, :]

    return torch.LongTensor(fea1d[:, :].transpose((1, 0)))

def get_msa_feature(msa_path,
                    msa_depth,
                    batch_converter = Alphabet.from_architecture('RNA MSA Transformer').get_batch_converter(),
                    mode = 'TopN'):

    msa_data = [read_msa(msa_path, msa_depth, mode)]

    return _get_msa_feature(msa_data, batch_converter)


if __name__ == '__main__':
    # import os
    # msa_dir = '/data1/rna/RNA3D/data/msa'
    # fas_dir = '/data1/rna/RNA3D/data/seq'

    # i = 0
    # for tgt in os.listdir(msa_dir):
    #     tokens = get_msa_feature(msa_path=f'{fas_dir}/{tgt}'.replace('.a3m','.seq'), msa_depth=32, mode='Rand')
    #     print(tokens)
    #     i+=1
    # print(i)

    msa_path = f'/public/home/taoshen/data/rna/RNA3D/Ash_Total_RNA_test/msa/3j2cM.a3m'
    tokens = get_msa_feature(msa_path, msa_depth=32, mode='Rand')
    print(tokens.shape)



