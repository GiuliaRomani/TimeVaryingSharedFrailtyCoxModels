#!/bin/bash

make distclean

make

cd Src

./main ../Data/DataTool/DataToolFile.txt ../Data/DataIndividuals/DataIndividualsFile.txt
