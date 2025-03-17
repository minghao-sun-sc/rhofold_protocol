#!/bin/bash

# the path to this shell script
FILE=`readlink -e $0`

# the path to the directory of this script and other binaries
bindir=`dirname $FILE`

# the path to the parent directory of the binaries
rootdir=`dirname $bindir`

# navigate to the database directory
cd "$rootdir/db"

echo "downloading RNAcentral"
wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_species_specific_ids.fasta.gz -O rnacentral_species_specific_ids.fasta.gz 
wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/release_notes.txt -O release_notes.txt
wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/rfam/rfam_annotations.tsv.gz -O rfam_annotations.tsv.gz

echo "downloading rfam"
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.full_region.gz -O Rfam.full_region.gz 
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz -O Rfam.cm.gz 
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/version.txt.gz -O version.txt.gz
gzip -d -f version.txt.gz

echo "downloading nt"
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz -O nt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz.md5 -O nt.gz.md5
md5sum -c nt.gz.md5

# get the nt database's html page for extracting BLAST database download links
# note that the raw html page only contains the file names, we need to manually add the base url
wget https://ftp.ncbi.nih.gov/blast/db/ -O nt.html

base_url="https://ftp.ncbi.nih.gov/blast/db/"

for file_name in `grep -Po 'href="\K(nt\.[^"]+\.tar\.gz)(?=")' nt.html | uniq`; do

    # for url in `grep -ohP "https://ftp.ncbi.nih.gov[:\d]*/blast/db/nt.\d+.tar.gz" nt.html |uniq`;do
    # filename=`basename $url`

    url="${base_url}${file_name}"
    echo "downloading nt: $url"
    wget $url -O "$file_name"
    wget $url.md5 -O "${file_name}.md5"
    md5sum -c "${file_name}.md5"

done

rm *md5

echo "download complete"

