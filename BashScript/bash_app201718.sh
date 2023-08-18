#!/bin/bash

make distclean

make

cd Src

./main ../Data/DataTool/DataToolFile201718.txt ../Data/DataIndividuals/DataIndividualsFile201718.txt
