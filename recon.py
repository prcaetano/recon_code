#!/usr/bin/env python
import sys
import os
import yaml
from ctypes import *

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

    input_path = config["input_reconstruction"]["input_path"]
    data_file = os.path.join(input_path,
                             config["input_reconstruction"]["data_fname"]).encode()
    randoms1_file = os.path.join(input_path,
                                 config["input_reconstruction"]["randoms1_fname"]).encode()
    randoms2_file = os.path.join(input_path,
                                 config["input_reconstruction"]["randoms2_fname"]).encode()

    output_path = config["output_reconstruction"]["output_path"]
    output_data_file = os.path.join(output_path,
                                    config["output_reconstruction"]["data_fname"]).encode()
    shifted_randoms_file = os.path.join(output_path,
                                        config["output_reconstruction"]["shifted_randoms_fname"]).encode()

    b = c_float(config["input_reconstruction"]["bias"])
    f = c_float(config["input_reconstruction"]["f"])
    Rf = c_float(config["input_reconstruction"]["R_smooth"])
    Om = c_float(config["input_reconstruction"]["Omega_m"])

    is_sim_data = c_bool(config["input_reconstruction"]["is_simulation"])

    recon(data_file, randoms1_file, randoms2_file,
          output_data_file, shifted_randoms_file,
          b, f, Rf, Om, is_sim_data)


