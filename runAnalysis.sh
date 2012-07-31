#!/bin/bash

echo "Setting-up Environment"
startdir=`pwd`
cd ~/work
source setupEnv.sh
cd $startdir

echo "Building Analysis Package"
scons

echo "Downloading Data Files"
./downloadMyDatasetsFromS3.py

for i in WHmumu ZHmumu vbfHmumu3 Zmumujets ZmumujetsMgt100 ggHmumu; do
echo "analyzing $i"
./analyzer $i.root $i/*.root
done
echo "done"
