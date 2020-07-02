#!/usr/bin/env python
import sys
import os
import yaml
import subprocess
from ctypes import *


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

    is_sim_data = config["reconstruction_config"]["is_simulation"]
    is_redshift_space = config["reconstruction_config"]["is_redshift_space"]
    try:
        read_weight = config["reconstruction_config"]["read_weight"]
    except KeyError:
        read_weight = False

    rsd_convention_enum = {"RecSym": 0, "RecIso": 1}
    try:
        rsd_convention = config["reconstruction_config"]["rsd_convention"]
    except KeyError:
        rsd_convention = "RecSym"
    rsd_convention_int = rsd_convention_enum[rsd_convention]

    # Import external C code
    try:
        lib = CDLL("recon.so")
    except OSError:
        raise RuntimeError("Couldn't load the library recon.so. Make sure it was "
              "compiled and LD_LIBRARY_PATH contains the directory holding it.")
    recon = lib.recon
    recon.argtypes = [c_char_p, c_char_p, c_char_p,
                      c_char_p, c_char_p,
                      c_float, c_float, c_float, c_float, c_bool,
                      c_bool, c_bool, c_int]
    recon.restype = c_int

    recon(data_file.encode(), randoms1_file.encode(), randoms2_file.encode(),
          output_data_file.encode(), shifted_randoms_file.encode(),
          c_float(b), c_float(f), c_float(Rf), c_float(Om), c_bool(is_sim_data),
          c_bool(is_redshift_space), c_bool(read_weight), c_int(rsd_convention_int))

