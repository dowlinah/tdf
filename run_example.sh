#!/bin/bash

# Set Up Example
cd poisson_example/src;
ffc -l dolfin forms.ufl
ln -s ../../TLV_file.cpp .
ln -s ../../TLV_file.hpp .

# Build Example
cd ../
mkdir build
cd build;
cmake ..;
make

./poisson_example
ls -lah *.tlv *.h5 meta*
