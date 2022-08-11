# FEniCS TLV File

Authors: Anthony Dowling, Lin Jiang, Dr. Ming-Cheng Cheng, Dr. Yu Liu

TLV File format implementation used for storing FEniCS Function objects.

## Using the TLV File in a FEniCS C++ Program

In FEniCS, the HDF5 file format can be used as follows:

```
Function u;
HDF5File in_file = HDF5File(mesh->mpi_comm(), "file.h5", "r");
input_file.read(u);
HDF5File out_file = HDF5File(mesh->mpi_comm(), "outfile.h5", "w");
output_file.write(u);
```

Our TLV format can be used similarly:

```
Function u;
TLVFile in_file = TLVFile(mesh->mpi_comm(), "file.tlv", "meta");
input_file.cache_metadata = true; // enables caching of metadata
input_file.read(u);
TLVFile out_file = TLVFile(mesh->mpi_comm(), "out.tlv", "out_meta");
output_file.save_metadata = true; // only necessary for first write
output_file.write(u);
```

The arguments are very similar to dolfin's HDF5
implementation, except for the third.  The third argument to
the constructor is the metadata prefix. When data is saved
using the TLV format, three metadata vectors are saved in
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
true causes the TLVFile object to read the metadata vectors
when it reads the first solution, then it uses the
previously read metadata to order subsequently read solution
vectors into the mesh.

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
the TLV file class.
