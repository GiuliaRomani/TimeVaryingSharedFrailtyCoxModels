#!/bin/bash

make distclean

make

./main Data/DataToolFile201718.txt Data/DataIndividualsFile201718.txt
