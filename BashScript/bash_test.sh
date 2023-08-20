#!/bin/bash

cd ..

make distclean

make

cd Src

clear 

./main ../Data/DataTool/DataToolFileTest.txt ../Data/DataIndividuals/DataIndividualsFileTest.txt

