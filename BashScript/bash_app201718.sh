#!/bin/bash

cd ..

make distclean

make

cd Src

clear 

./main ../Data/DataTool/DataToolFile201718.txt ../Data/DataIndividuals/DataIndividualsFile201718.txt
