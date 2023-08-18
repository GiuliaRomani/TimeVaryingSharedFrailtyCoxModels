#!/bin/bash

cd ..

make distclean

make

cd Src

clear 

./main ../Data/DataTool/DataToolFile.txt ../Data/DataIndividuals/DataIndividualsFile.txt
