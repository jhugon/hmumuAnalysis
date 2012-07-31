#!/bin/bash

# Can be run like: ./run.sh [bucket] [outdataset] [indataset1] [indataset2] ...

echo "Setting-up Environment"
startdir=`pwd`
cd ~/work
source setupEnv.sh
cd $startdir

echo "Building Analysis Package"
scons

echo "Downloading Data Files"
./downloadMyDatasetsFromS3.py "$@"

for i in ${@:3}; do
echo "analyzing $i"
./analyzer $i.root $i/*.root
done

echo "Uploading Results to S3..."
./uploadResultsToS3.py "$@"
echo "done"
