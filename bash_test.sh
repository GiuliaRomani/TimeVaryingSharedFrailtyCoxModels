#!/bin/bash

make distclean

make

cd Src

./main ../Data/DataTool/DataToolFileTest.txt ../Data/DataIndividuals/DataIndividualsFileTest.txt
