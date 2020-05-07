#!/usr/bin/env python
import sys
import os
import yaml
from ctypes import *

# Import external C code
try:
    lib = CDLL("recon.so")
except OSError:
    raise RuntimeError("Couldn't load the library recon.so. Make sure it was "
          "compiled and LD_LIBRARY_PATH contains the directory holding it.")
recon = lib.recon
recon.argtypes = [c_char_p, c_char_p, c_char_p,
                  c_char_p, c_char_p,
                  c_float, c_float, c_float, c_float, c_bool]
recon.restype = c_int


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
                             config["input"]["data_fname"]).encode()
    randoms1_file = os.path.join(input_path,
                                 config["input"]["randoms1_fname"]).encode()
    randoms2_file = os.path.join(input_path,
                                 config["input"]["randoms2_fname"]).encode()

    output_path = config["output"]["output_path"]
    output_data_file = os.path.join(output_path,
                                    config["output"]["displaced_catalog_fname"]).encode()
    shifted_randoms_file = os.path.join(output_path,
                                        config["output"]["shifted_random_catalog_fname"]).encode()

    b = c_float(config["reconstruction_config"]["bias"])
    f = c_float(config["reconstruction_config"]["f"])
    Rf = c_float(config["reconstruction_config"]["R_smooth"])
    Om = c_float(config["reconstruction_config"]["Omega_m"])

    is_sim_data = c_bool(config["reconstruction_config"]["is_simulation"])

    recon(data_file, randoms1_file, randoms2_file,
          output_data_file, shifted_randoms_file,
          b, f, Rf, Om, is_sim_data)


