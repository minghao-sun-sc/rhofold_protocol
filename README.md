
<p align="center">

  <h3 align="center">A language-model-based deep learning platform for predicting RNA 3D structure</h3>

  <p align="center">
    Supporting code for the paper
  </p>
</p>

## Table of contents


<!-- ABOUT THE PROJECT -->
## About this repository

This codebase contains the relevant codes and example data associated with the paper titled "A language-model-based deep learning platform for predicting RNA 3D structure". The protocol is based on our previous research on RNA 3D structure modeling, which is published in Nature Methods (<a href='[www.google.com](https://www.nature.com/articles/s41592-024-02487-0)'>full text available</a>).

This protocol comprises four stages.
- Stage 1: prepare data
- Stage 2: generate RNA-FM embeddings
- Stage 3: perform structure inference
- Stage 4: analyze the evaluate the prediction

Detailed steps for each stage and suggestions on how to choose the configuration at each step can be found in the protocol. 

<!-- File organization -->
## Content of the repository
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

### Working example dataset




## Contact

For questions or comments, please feel free to post an issue or reach the author at jmwang@link.cuhk.edu.hk.
