# FEniCS TDF File

Authors: Anthony Dowling, Lin Jiang, Dr. Ming-Cheng Cheng, Dr. Yu Liu

Department of Electrical and Computer Engineering, Clarkson University.

Emails: {dowlinah,jiangl2,mcheng,yuliu}@clarkson.edu

TLV-based File format implementation used for storing FEniCS Function objects.

## Using the TDF File in a FEniCS C++ Program

In FEniCS, the HDF5 file format can be used as follows:

```
Function u;
HDF5File in_file = HDF5File(mesh->mpi_comm(), "file.h5", "r");
input_file.read(u);
HDF5File out_file = HDF5File(mesh->mpi_comm(), "outfile.h5", "w");
output_file.write(u);
```

Our TDF format can be used similarly:

```
Function u;
TDFFile in_file = TDFFile(mesh->mpi_comm(), "file.tdf", "meta");
input_file.cache_metadata = true; // enables caching of metadata
input_file.read(u);
TDFFile out_file = TDFFile(mesh->mpi_comm(), "out.tdf", "out_meta");
output_file.save_metadata = true; // only necessary for first write
output_file.write(u);
```

The arguments are very similar to dolfin's HDF5
implementation, except for the third.  The third argument to
the constructor is the metadata prefix. When data is saved
using the TDF format, three metadata vectors are saved in
separate files using this prefix. The suffixes are added
automatically. The same metadata prefix that was used to
write the data should be used during loading to ensure
consistency.

When using the file format, make sure to include the header
file, and compile the CPP file into your program.

### Options

The `save_metadata` option can be set to `false` after
writing the first time-step's solution to save time on later
writes and only write the solution data vector instead of
also writing the metadata vectors.

The `cache_metadata` option enables the use of cached
metadata when reading multiple solution files. Setting it to
true causes the TDFFile object to read the metadata vectors
when it reads the first solution, then it uses the
previously read metadata to order subsequently read solution
vectors into the mesh.

## Running Poisson Equation Example

Included is an example FEniCS program that computes the
Poisson Equation. Assuming all dependencies are installed on
your system, running `bash run_example.sh` will execute the
example program and display the sizes of all of the metadata
files, the solution TDF file, and the HDF5 file that stores
the same information. The included `reset_example.sh` script
can be used to remove the files generates by
`run_example.sh`

## Generating Documentation

To generate the documentation for the user-level API, make
sure that doxygen, make, and latex are installed on your
system, then, in the root folder, run:

```
doxygen
cd docs/latex
make
```

and that will generate refman.pdf, which details the use of
the TDF file class.
