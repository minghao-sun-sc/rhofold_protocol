This directory is for keeping the downloaded and curated MSA databases for RhoFold+ described in the protocol.

To download and curate the MSA databases, run the following command from the root directory of this repository.
```
chmod +x ./msa_database/bin/*
bash ./msa_database/bin/builddb.sh
```
This will automatically save the checkpoint to this directory and rename it to the default name configured in RhoFold+.

This step is optional if users have already generated the MSAs or do not need to use MSAs.

For details, please consult the "Equipment setup" section in the protocol.
