#!/bin/bash

make distclean

make

./main Data/DataToolFileTest.txt Data/DataIndividualsFileTest.txt
