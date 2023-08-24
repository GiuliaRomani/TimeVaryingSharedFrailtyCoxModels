#!/bin/bash

#------------------------------------------------------------------
# Bash script for executing a model with test data and test 
# information provided in the .txt files.
#
# If you have compiled the codes and you just want to play with the 
# models, comment the following lines: make distclean, make docs, 
# make, make docsclean and execute the script on terminal.
#------------------------------------------------------------------

# Change directory to go into the main one, where Makefile is contained
cd ..

# Clean the sub-directories
make distclean

# Create doxygen documentation. It directly opens the index.html file
make docs

# Compile 
make

# Change directory and go in Src, where the executable ./main is contained
cd Src

# Clear the terminal
clear 

# Execute
./main ../Data/DataTool/DataToolFileTest.txt ../Data/DataIndividuals/DataIndividualsFileTest.txt

# Remove all doxygen documentation (Doc/doxygen folder)
#make docsclean

