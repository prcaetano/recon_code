#!/usr/bin/env python
import sys
import os
import yaml
import subprocess
import tempfile
import numpy as np
import ndfilehandler as fh
import h5py
from astropy.table import Table
from astropy.io import fits
#from ctypes import *


def read_fits_file(fits_fname, fh_fname):
    """
    Reads relevant columns from fits file and converts to a filehandler file
    to be passed to the C code
    """
    data = {}
    hdul = fits.open(fits_fname)
    columns = ["RA", "DEC", "Z"]
    for column in columns:
        try:
            data[column.lower()] = hdul[1].data.field(column)
        except ValueError:
            print("ERROR: ", fits_fname, " should contain an ", column, " column.",
                  file=sys.stderr)
            exit(1)
    fh.write_file(fh_fname, data)


def write_fits_file(fh_fname, output_fits_fname, input_fits_fname):
    """
    Writes back to original fits file the additional columns computed.
    """
    data = fh.read_file(fh_fname,None)
    tab = Table.read(input_fits_fname, memmap=True)
    for column in ["RA", "DEC", "Z", "X_COORD", "Y_COORD", "Z_COORD"]:
        tab[column+'_RECON'] = data[column.lower()]
    tab.write(output_fits_fname)


def read_hdf5_file(hdf5_fname, hdf5_path, fh_fname):
    """
    Reads relevant columns from hdf5 file and converts to a filehandler file
    to be passed to the C code
    """
    data = {}
    f = h5py.File(hdf5_fname, "r")
    columns = ["RA", "DEC", "Z"]
    for column in columns:
        dataset = os.path.join(hdf5_path, column.lower())
        try:
            data[column.lower()] = np.array(f[dataset])
        except KeyError:
            print("ERROR: ", hdf5_fname, " should contain an dataset ", column.lower(),
                  " under the path ", hdf5_path, ".", file=sys.stderr)
            exit(1)
    fh.write_file(fh_fname, data)
    f.close()


def write_hdf5_file(fh_fname, output_hdf5_fname, input_hdf5_fname, hdf5_path):
    """
    Writes back to original hdf5 file the additional columns computed.
    """
    if input_hdf5_fname != output_hdf5_fname:
        subprocess.run(["cp " + input_hdf5_fname + " " + output_hdf5_fname], shell=True)

    f = h5py.File(output_hdf5_fname, "r+")
    data = fh.read_file(fh_fname, None)
    columns = ["RA", "DEC", "Z", "X_COORD", "Y_COORD", "Z_COORD"]
    for column in columns:
        full_path = os.path.join(hdf5_path, column.lower() + "_recon")
        f[full_path] = data[column.lower()]
    f.close()


if __name__ == "__main__":
    try:
        config_file = sys.argv[1]
    except IndexError:
        print("Usage: recon.py <YAML configuration file>")
        exit()

    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    input_path = config["input"]["input_path"]
    data_file = os.path.join(input_path,
                             config["input"]["data_fname"])
    randoms1_file = os.path.join(input_path,
                                 config["input"]["randoms1_fname"])
    randoms2_file = os.path.join(input_path,
                                 config["input"]["randoms2_fname"])

    output_path = config["output"]["output_path"]
    output_data_file = os.path.join(output_path,
                                    config["output"]["displaced_catalog_fname"])
    shifted_randoms_file = os.path.join(output_path,
                                        config["output"]["shifted_random_catalog_fname"])

    b = config["reconstruction_config"]["bias"]
    f = config["reconstruction_config"]["f"]
    Rf = config["reconstruction_config"]["R_smooth"]
    Om = config["reconstruction_config"]["Omega_m"]
    rsd_convention = config["reconstruction_config"]["RSD_convention"]
    try:
        hdf_preffix = config["reconstruction_config"]["hdf5_preffix"]
    except KeyError:
        hdf_preffix = "/"

    data_fmt = os.path.splitext(data_file)[-1][1:]
    randoms1_fmt = os.path.splitext(randoms1_file)[-1][1:]
    randoms2_fmt = os.path.splitext(randoms2_file)[-1][1:]

    hdf_exts = ['hdf', 'h4', 'hdf4', 'he2', 'h5', 'hdf5', 'he5']

    with tempfile.TemporaryDirectory(dir=output_path) as tmpdirname:
        if data_fmt == "fits":
            new_data_file = os.path.join(tmpdirname, "data.fh")
            read_fits_file(data_file, new_data_file)
            data_file, old_data_file = new_data_file, data_file

            new_output_data_file = os.path.join(tmpdirname, "data_out.fh")
            old_output_data_file, output_data_file = output_data_file, new_output_data_file
        elif data_fmt in hdf_exts:
            new_data_file = os.path.join(tmpdirname, "data.fh")
            read_hdf5_file(data_file, hdf_preffix, new_data_file)
            data_file, old_data_file = new_data_file, data_file

            new_output_data_file = os.path.join(tmpdirname, "data_out.fh")
            old_output_data_file, output_data_file = output_data_file, new_output_data_file
        elif data_fmt == "fh":
            pass
        else:
            raise RuntimeError("Unsupported format ", data_fmt, " for input data.")

        if randoms1_fmt == "fits":
            new_randoms1_file = os.path.join(tmpdirname, "randoms1.fh")
            read_fits_file(randoms1_file, new_randoms1_file)
            randoms1_file = new_randoms1_file
        elif randoms1_fmt in hdf_exts:
            new_randoms1_file = os.path.join(tmpdirname, "randoms1.fh")
            read_hdf5_file(randoms1_file, hdf_preffix, new_randoms1_file)
            randoms1_file = new_randoms1_file
        elif randoms1_fmt == "fh":
            pass
        else:
            raise RuntimeError("Unsupported format ", randoms1_fmt, " for the randoms.")

        if randoms2_fmt == "fits":
            new_randoms2_file = os.path.join(tmpdirname, "randoms2.fh")
            read_fits_file(randoms2_file, new_randoms2_file)
            randoms2_file, old_randoms2_file = new_randoms2_file, randoms2_file

            new_shifted_randoms_file = os.path.join(tmpdirname, "randoms_out.fh")
            old_shifted_randoms_file, shifted_randoms_file = shifted_randoms_file, new_shifted_randoms_file
        elif randoms2_fmt in hdf_exts:
            new_randoms2_file = os.path.join(tmpdirname, "randoms2.fh")
            read_hdf5_file(randoms2_file, hdf_preffix, new_randoms2_file)
            randoms2_file, old_randoms2_file = new_randoms2_file, randoms2_file

            new_shifted_randoms_file = os.path.join(tmpdirname, "randoms_out.fh")
            old_shifted_randoms_file, shifted_randoms_file = shifted_randoms_file, new_shifted_randoms_file
        elif randoms2_fmt == "fh":
            pass
        else:
            raise RuntimeError("Unsupported format ", randoms2_fmt, " for the randoms.")

        ## Import external C code
        #try:
        #    lib = CDLL("recon.so")
        #except OSError:
        #    raise RuntimeError("Couldn't load the library recon.so. Make sure it was "
        #          "compiled and LD_LIBRARY_PATH contains the directory holding it.")
        #recon = lib.recon
        #recon.argtypes = [c_char_p, c_char_p, c_char_p,
        #                  c_char_p, c_char_p,
        #                  c_float, c_float, c_float, c_float, c_char_p]
        #recon.restype = c_int

        #recon(data_file.encode(), randoms1_file.encode(), randoms2_file.encode(),
        #      output_data_file.encode(), shifted_randoms_file.encode(),
        #      c_float(b), c_float(f), c_float(Rf), c_float(Om), reciso.encode())

        cmd = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), "recon")
        args = "{} {} {} {} {} {} {} {} {} {}".format(data_file, randoms1_file,
                                                      randoms2_file, output_data_file,
                                                      shifted_randoms_file, b, f, Rf, Om,
                                                      rsd_convention)
        subprocess.run([cmd + " " + args], shell=True)

        if data_fmt == "fits":
            write_fits_file(fh_fname=output_data_file,
                            output_fits_fname=old_output_data_file,
                            input_fits_fname=old_data_file)
        elif data_fmt in hdf_exts:
            write_hdf5_file(fh_fname=output_data_file,
                            output_hdf5_fname=old_output_data_file,
                            input_hdf5_fname=old_data_file,
                            hdf5_path=hdf_preffix)


        if randoms2_fmt == "fits":
            write_fits_file(fh_fname=shifted_randoms_file,
                            output_fits_fname=old_shifted_randoms_file,
                            input_fits_fname=old_randoms2_file)
        elif randoms2_fmt in hdf_exts:
            write_hdf5_file(fh_fname=shifted_randoms_file,
                            output_hdf5_fname=old_shifted_randoms_file,
                            input_hdf5_fname=old_randoms2_file,
                            hdf5_path=hdf_preffix)


