#!/bin/bash

# this is the path to this shell script itself
FILE=`readlink -e $0`

# this is the directory of this shell script, along with other binaries
bin_dir=`dirname $FILE`

# this is the directory of the database <=> parent directory of the binary directory 
root_dir=`dirname $bin_dir`

echo "path to this script: $FILE"
echo "path to the binaries' directory: $bin_dir"
echo "path to the database directory: $root_dir"

# navigate to the database's directory
cd $root_dir

$bin_dir/download.sh
$bin_dir/curate.sh

