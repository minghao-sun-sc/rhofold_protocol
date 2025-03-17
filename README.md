
<p align="center">

  <h3 align="center">A language-model-based deep learning platform for predicting RNA 3D structure</h3>

  <p align="center">
    Supporting code for the paper
  </p>
</p>

## Table of contents
* [About this repository](#about-this-repository)
* [Content of this repository](#content-of-this-repository)
  * [Codebase organization](#codebase-organization)
  * [Working example dataset information](#working-example-data-information)
  * [Version of the codes in this repository](#version-of-the-codes-in-this-repository)
* [Quick start](#quick-start)
* [Contact](#contact)

<!-- ABOUT THE PROJECT -->
## About this repository

This codebase contains the relevant codes and example data associated with the paper titled "A language-model-based deep learning platform for predicting RNA 3D structure". The protocol is based on our previous research on RNA 3D structure modeling, which is published in Nature Methods (<a href='https://www.nature.com/articles/s41592-024-02487-0'>full text available</a>).

This protocol comprises four stages.
- Stage 1: prepare data
- Stage 2: generate RNA-FM embeddings
- Stage 3: perform structure inference
- Stage 4: analyze the evaluate the prediction

Detailed steps for each stage and suggestions on how to choose the configuration at each step can be found in the protocol. 

<!-- File organization -->
## Content of this repository
### Codebase organization

This codebase contains the following directories:
<ul>
<li><code>checkpoints/</code>: for downloading and keeping the model checkpoints, including the pre-trained RhoFold+ model or the RNA-FM checkpoint.</li>
<li><code>data/</code>: keeping the data as working examples for this protocol, which includes:</li>

1. a structure prediction example at <code>data/rhofold/3owz_A/</code> from Protein Data Bank (PDB) and

2. a sequence dataset example for testing embedding generation at <code>data/rnafm/rf02684/seqs.fasta</code> from Rfam.

<li><code>msa_database</code>: for downloading and keeping the MSA databases, including RNAcentral, Rfam, and nt. 

1. The <code>msa_database/bin</code> directory contains the scripts necessary for downloading and building the databases. 

2. The <code>msa_database/db</code> directory holds all the downloaded databases.</li>
<li><code>results</code>: for keeping the output from this protocol.</li>
<li><code>rhofold</code>: main module of the RhoFold+ model, which is adapted from the original RhoFold+ model to streamline the workflow of this protocol.</li>
<li><code>rmsa</code>: for keeping the rMSA tool for MSA search, which is cloned from the official rMSA2 release.</li>
<li><code>scripts</code>: for keeping the additional codes shown in the protocol.</li>
<li><code>integrate_rnafm.py</code>: example of Step 5B in the protocol for directly integrating RNA-FM into a Python script as a package.</li>

</ul>

As mentioned in the protocol, inputting relative paths can be a common source of errors

### Working example data information

This repository provides two sets of data as working examples for performing the protocol.

1. <code>data/rhofold/3owz_A/</code> for RNA 3D structure prediction. This dataset includes a structure from PDB (`.pdb` file) and its associated sequence (`.fasta` file) and secondary structure (extracted by DSSR, `.npy` file, in contact map format). There is also an `.afa` file for the constructed multiple sequence alignment (MSA).

2. <code>data/rnafm/rf02684/seqs.fasta</code> for RNA-FM embedding generation. This dataset includes a set of RNA sequences derived from Rfam, using the RF02684 Twister family's seed sequences.

### Customized scripts for this protocol

To facilitate the analysis and evaluation of the results, we have incorporated some additional scripts in this repository under the `scripts/` directory. The codes of these scripts are also shown in the protocol in their respective steps with command line argument-passed paths replaced by fixed ones.

Stage 2
- `integrate_rnafm.py` (Step 5B)

Stage 4
- `parse_plddt.py` (Step 7)
- `add_plddt_bfactor.py` (Step 8)
- `visualize_plddt.py` (Step 9)
- `eval_lddt.py` (Step 11)
- `visualize_ss_prob.py` (Step 12)
- `visualize_ss.py` (Step 13)

### Version of the codes in this repository

This codebase adapts the codes from our previous work RhoFold+ and the codes from rMSA2. The exact version of these two methods are as follows.

- `RhoFold+` from its <a href='https://github.com/ml4bio/RhoFold/tree/df930033dd40c6c3f923dcafcdc16cf50eb742c8'>official GitHub repo</a>, commit `df93003`.
- `rMSA` from its <a href='https://github.com/kad-ecoli/rMSA2'>official GitHub repo</a>, commit `3fa7c22`.

Note that the codes for RhoFold+ has been adapted from the official release, including command-line arguments, checkpoint loading, and RNA-FM embedding generation.

Users may also be aware of the <a href='https://github.com/ml4bio/RNA-FM'>RNA-FM repository</a>. In this protocol, since our main focus is on 3D structure modeling, we only use the codes from the RhoFold+ repository, which includes the same RNA-FM model and can also perform its core functions such as embedding generation.

## Quick start

The details of each step is described in the protocol. Here, we briefly outline how to download and install the dependencies for this repository, corresponding to the equipment setup section's mandatory steps in the protocol.

```
git clone https://github.com/WangJiuming/rhofold_protocol.git
cd rhofold_protocol
```
Then create the cond environment from the `.yml` configuration file.
```
conda env create -f environment.yml
```
This should install all the dependencies necessary for this protocol. To activate the environment, simply run the following.
```
conda activate rhofold_protocol
```
The pre-trained model checkpoint used in this protocol can be downloaded from the server via the following. The estimated download time is 5 minutes.
```
wget https://proj.cse.cuhk.edu.hk/aihlab/RhoFold/api/download?filename=RhoFold_pretrained.pt -O ./checkpoints/rhofold_pretrained_params.pt
```

## Contact

For questions or comments, please feel free to post an issue or reach the author at jmwang@link.cuhk.edu.hk.
