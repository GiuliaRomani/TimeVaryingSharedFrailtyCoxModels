#!/bin/bash

cd ..

make distclean

make

cd Src

clear 

./main ../Data/DataTool/DataToolFile2018.txt ../Data/DataIndividuals/DataIndividualsFile2018.txt
