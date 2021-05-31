#!/bin/bash

cd src
tar xzf ARACNE.src.tar.gz
cd ARACNE
make clean
make
cd ../../
pwd
ln -s ../src/ARACNE bin/ARACNE

[[ -f "bin/ARACNE/usage.txt" ]] && echo "Finished successfully." && exit 0
echo "Something was wrong"
exit 15
