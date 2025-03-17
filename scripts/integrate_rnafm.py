import fm
import torch
import numpy as np
from Bio import SeqIO

# Load data from a FASTA file
fasta_path = './data/rnafm/rf02684/seqs.fasta'
records = list(SeqIO.parse(fasta_path, 'fasta'))

# Prepare the input data format
data = [(record.id, str(record.seq)) for record in records]

# Get the device, CPU or GPU
device = 'cuda' if torch.cuda.is_available() else 'cpu'

# Load RNA-FM model
fm_ckpt_path = ''  # specify the path here or omit
fm_model, alphabet = fm.pretrained.rna_fm_t12(fm_ckpt_path)
batch_converter = alphabet.get_batch_converter()

fm_model.to(device)  # use GPU if available

# Begin generating the embeddings
fm_model.eval()
batch_labels, batch_strs, batch_tokens = batch_converter(data)

with torch.no_grad():
    results = fm_model(batch_tokens.to(device), repr_layers=[12])
# results is a dictionary with keys: logits, representations
emb = results['representations'][12].cpu().numpy()  # emb should have shape (B, L+2, 640)

# Truncanate the embeddings to remove the <SOS> and <EOS> tokens
emb = emb[:, 1:-1, :]  # after truncating, the shape becomes (B, L, 640)

print('Embedding has been generted, now you can use it for your downstream tasks.')
print(f'Embedding has shape: {emb.shape}')  # (B, L, 640)
