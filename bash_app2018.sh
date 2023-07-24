#!/bin/bash

make distclean

make

./main Data/DataToolFile2018.txt Data/DataIndividualsFile2018.txt
