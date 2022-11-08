#!/bin/bash

# Set Up Example
cd poisson_example/src;
ffc -l dolfin forms.ufl
ln -s ../../TDF_file.cpp .
ln -s ../../TDF_file.hpp .

# Build Example
cd ../
mkdir build
cd build;
cmake ..;
make

./poisson_example
ls -lah *.tdf *.h5 meta*
