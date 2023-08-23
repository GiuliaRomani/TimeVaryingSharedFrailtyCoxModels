#!/bin/bash

cd ..

make distclean

make

cd Src

clear 

./main ../Data/DataTool/DataToolFile2010.txt ../Data/DataIndividuals/DataIndividualsFile2010.txt
