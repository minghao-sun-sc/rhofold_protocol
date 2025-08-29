#/bash/bin

# Reference:https://academic.oup.com/bioinformatics/article/35/21/4459/5480133?login=false
# Download: https://seq2fun.dcmb.med.umich.edu//RNA-align/download.html

if [ $# -lt 2 ]; then
  echo "usage: tmscore_rna.sh inp.pdb ntv.pdb tag"
  echo "example: tmscore_rna.sh farfar2S_000001.pdb native_RNA/2RD2.pdb 2RD2"
  echo "========================================="
  echo "Results: tag TMscore RMSD"
  echo "Example: 4FB0 0.18916   36.44"
  exit 0;
fi

if [ $# -eq 2 ]; then
  tag="tag"
else
  tag=$3
fi

inp=$1
ntv=$2

tm=`/public/home/taoshen/code/Projects/StructurePrediction/rnafold/tools/bin/RNAalign/RNAalign $inp $ntv -TMscore 5 | grep "TM-score=" | awk '{print $2}' | tail -n 1`
rmsd=`/public/home/taoshen/code/Projects/StructurePrediction/rnafold/tools/bin/RNAalign/RNAalign $inp $ntv -TMscore 5 | grep "RMSD=" | cut -d "," -f 2| cut -d "=" -f 2 | tail -n 1`

echo "$tag $tm $rmsd"