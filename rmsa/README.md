This is the code repository for rMSA (RNA Multiple Sequence Alignment), an automated pipeline to search
and align homologs from RNAcentral and nt databases for a target RNA.

## Third party programs ##
The ``bin`` folder includes binaries precompiled for 64bit Linux for
the following programs.

* nhmmer (``bin/qnhmmer``) from HMMER 3.3
* cd-hit-est and cd-hit-est-2d from CD-HIT 4.8.1
* clustalo 1.2.4
* cmbuild, cmcalibrate, cmscan and cmsearch (``bin/qcmsearch``) from INFERNAL 1.1.3
* hhfilter and reformat.pl from HH-suite 2.0.15
* RNAfold from ViennaRNA 2.4.14
* dot2ct from RNAstructure 6.2
* plmc 2018-05-16
* PETfold 2.1

The output format of nhmmer and cmsearch are modifed from
``eslMSAFILE_STOCKHOLM`` or ``eslMSAFILE_PFAM`` to ``eslMSAFILE_A2M``.

## Run the program ##
```bash
./rMSA.pl seq.fasta
```
Run ``./rMSA.pl`` without command line argument to get full option list,
including alternative databases, temporary folders, number of CPU threads
(default is 1) and the secondary structure. Nucleotide U is converted to T.

## Citation
```bibtex
@article{zhang2023rmsa,
  title={rMSA: a sequence search and alignment algorithm to improve RNA structure modeling},
  author={Zhang, Chengxin and Zhang, Yang and Pyle, Anna Marie},
  journal={Journal of Molecular Biology},
  volume={435},
  number={14},
  pages={167904},
  year={2023},
  publisher={Elsevier}
}
```

## License ##
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.
