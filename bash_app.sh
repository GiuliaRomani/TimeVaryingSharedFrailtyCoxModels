#!/bin/bash

make distclean

make

./main Data/DataTool/DataToolFile.txt Data/DataIndividuals/DataIndividualsFile.txt
