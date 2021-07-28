# recon_code

Code to perform density field reconstruction given an object catalog
and two random files.  The goal of density field reconstruction is to
sharpen the baryon acoustic oscillation feature for measuring distances.
This code implements the standard, lowest-order algorithm presented in
Eisenstein et al. (2007) with periodic boundary conditions using a
multigrid relaxation technique with a full multigrid V-cycle based on
damped Jacobi iteration.
Inputs are a catalog of objects and two random catalogs which serve to
specify the selection function and act as the "shifted" fields. Output
are 'shifted' versions of the object and second random catalogs.

# Use

This version of the code accepts and writes out binary "files" in the
filehandler format [https://github.com/martinjameswhite/filehandler].
The syntax is:

`
    recon <data-file> <random-file> <random-file> <output-file> \
    <output-shifted-randoms-file> <bias> <f-growth> <R-filter> \
    <Omega m> <RSD Convention (RecIso|RecSym)>
`

There is also a python wrapper available, which runs the cpp code as
a subprocess and takes care of input/output convertion from fits and
hdf5 formats. The syntax in this case is:

`
    recon.py <YAML configuration file> \
    [--data_path <fits, hdf5 or fh file> [<fits, hdf5 or fh file> ...]] \
    [--randoms1_path <fits, hdf5 or fh file> [<fits, hdf5 or fh file> ...]] \
    [--randoms2_path <fits, hdf5 or fh file> [<fits, hdf5 or fh file> ...]] \
    [--output_data_path [<fits, hdf5 or fh file>]] \
    [--output_randoms_path [<fits, hdf5 or fh file>]] \
`

Where the YAML configuration file holds configuration regarding
cosmology, RSD Conventions, IO formats and basic catalog pre-
processing (check `example_config.yml` for the options available).
It is possible to use multiple split files for a catalog, simply by
specifying multiple space-separated paths in the call.

