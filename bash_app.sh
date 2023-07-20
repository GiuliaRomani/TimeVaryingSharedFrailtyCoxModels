#!/bin/bash

make distclean

make

./main DataToolFile.txt DataIndividualsFile.txt
