#!/bin/bash

make distclean

make

cd Src

./main ../Data/DataTool/DataToolFile2018.txt ../Data/DataIndividuals/DataIndividualsFile2018.txt
