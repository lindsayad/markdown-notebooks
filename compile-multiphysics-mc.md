# 10/12/17

So there are warnings with making on Mac about ranpath no symbols. This occurs
regardless of whether MC finds hdf5 or not.

If I use FC=h5pfc, then I get to 43% build and fail with the error:
```
f951: Fatal Error: Can't rename module file ‘include/source_header.mod0’ to
‘include/source_header.mod’: No such file or directory
```
If I use FC=mpif90 and split hdf5 and mpi in cmake, then I fail at
