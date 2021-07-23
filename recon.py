#!/usr/bin/env python
import sys
import os
import yaml
import subprocess
import tempfile
import ndfilehandler as fh
from astropy.table import Table
from astropy.io import fits
#from ctypes import *


def read_fits_file(fits_fname, fh_fname):
    """
    Reads relevant colums from fits file and converts to a filehandler file
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

    data_fmt = None
    randoms_fmt = None

    with tempfile.TemporaryDirectory(dir=output_path) as tmpdirname:
        if os.path.splitext(data_file)[-1] == '.fits':
            new_data_file = os.path.join(tmpdirname, "data.fh")
            read_fits_file(data_file, new_data_file)
            data_file, old_data_file = new_data_file, data_file
            data_fmt = 'fits'

            new_output_data_file = os.path.join(tmpdirname, "data_out.fh")
            old_output_data_file, output_data_file = output_data_file, new_output_data_file

        if os.path.splitext(randoms1_file)[-1] == '.fits':
            new_randoms1_file = os.path.join(tmpdirname, "randoms1.fh")
            read_fits_file(randoms1_file, new_randoms1_file)
            randoms1_file = new_randoms1_file

        if os.path.splitext(randoms2_file)[-1] == '.fits':
            new_randoms2_file = os.path.join(tmpdirname, "randoms2.fh")
            read_fits_file(randoms2_file, new_randoms2_file)
            randoms2_file, old_randoms2_file = new_randoms2_file, randoms2_file
            randoms_fmt = 'fits'

            new_shifted_randoms_file = os.path.join(tmpdirname, "randoms_out.fh")
            old_shifted_randoms_file, shifted_randoms_file = shifted_randoms_file, new_shifted_randoms_file

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

        if randoms_fmt == "fits":
            write_fits_file(fh_fname=shifted_randoms_file,
                            output_fits_fname=old_shifted_randoms_file,
                            input_fits_fname=old_randoms2_file)

