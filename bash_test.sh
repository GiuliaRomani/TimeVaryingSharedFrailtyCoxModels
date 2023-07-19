#!/bin/bash

make distclean

make

./main DataToolFileTest.txt DataIndividualsFileTest.txt
