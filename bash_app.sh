#!/bin/bash

make distclean

make

./main Data/DataToolFile.txt Data/DataIndividualsFile.txt
