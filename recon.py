#!/usr/bin/env python
import argparse
import sys
import os
import yaml
import subprocess
import tempfile
import numpy as np
import ndfilehandler as fh
import h5py
from glob import glob
from astropy.table import Table
from astropy.io import fits
#from ctypes import *


def read_fits_files(fits_fnames, fh_fname, columns, row_mask):
    """
    Reads relevant columns from fits file and converts to a filehandler file
    to be passed to the C code
    """
    data = {}

    # Processing row_mask string to apply row selection
    if row_mask is not None:
        while (row_mask.find("${") != -1):
            col_label = row_mask[row_mask.find("${"):].split("}")[0][2:]
            sub_expression = "hdul[1].data.field('{}')".format(col_label)
            row_mask = row_mask.replace("${"+col_label+"}", sub_expression)

    masks = []
    sizes = np.zeros(len(fits_fnames)+1, dtype=int)
    for i, fits_fname in enumerate(fits_fnames):
        hdul = fits.open(fits_fname)
        if row_mask is not None:
            try:
                mask = eval(row_mask)
            except:
                raise RuntimeError("Invalid row selection specified.")
            sizes[i+1] = sizes[i] + mask.astype(int).sum()
            masks.append(mask)
        else:
            sizes[i+1] = sizes[i] + len(hdul[1].data)
        hdul.close()

    for field in columns:
        data[field] = np.zeros(sizes[-1], dtype='f8')

    for i, fits_fname in enumerate(fits_fnames):
        hdul = fits.open(fits_fname)
        for field in columns:
            col = columns[field]
            try:
                col_data = hdul[1].data.field(col)
                if row_mask is not None:
                    col_data = col_data[masks[i]]
                data[field][sizes[i]:sizes[i+1]] = col_data
            except ValueError:
                raise RuntimeError(fits_fname, " should contain an ", col, " column.")
        hdul.close()
    fh.write_file(fh_fname, data)


def write_fits_file(input_fh_fname, output_fh_fname, output_fits_fname,
                    input_columns, output_columns):
    """
    Converts output fh file to fits, keeping only the specified columns
    """
    output_data = fh.read_file(output_fh_fname, list(output_columns))
    input_data = fh.read_file(input_fh_fname, list(input_columns))

    data = {}
    for col in input_columns:
        out_col = input_columns[col]
        data[out_col] = input_data.pop(col)
    for col in output_columns:
        out_col = output_columns[col]
        data[out_col] = output_data.pop(col)
    Table(data).write(output_fits_fname)


def read_hdf5_files(hdf5_fnames, hdf5_path, fh_fname, columns, row_mask):
    """
    Reads relevant columns from hdf5 file and converts to a filehandler file
    to be passed to the C code
    """
    data = {}

    # Processing row_mask string to apply row selection
    if row_mask is not None:
        while (row_mask.find("${") != -1):
            col_label = row_mask[row_mask.find("${"):].split("}")[0][2:]
            sub_expression = "f[os.path.join(hdf5_path, '{}')][:]".format(col_label)
            row_mask = row_mask.replace("${"+col_label+"}", sub_expression)

    masks = []
    sizes = np.zeros(len(hdf5_fnames)+1, dtype=int)
    for i, hdf5_fname in enumerate(hdf5_fnames):
        f = h5py.File(hdf5_fname, "r")
        if row_mask is not None:
            try:
                mask = eval(row_mask)
            except:
                raise RuntimeError("Invalid row selection specified.")
            sizes[i+1] = sizes[i] + mask.astype(int).sum()
            masks.append(mask)
        else:
            # Assuming hdf5 file well formed: all columns w/ same size
            sizes[i+1] = sizes[i] + f[os.path.join(hdf5_path, columns['ra'])].len()
        f.close()

    for field in columns:
        data[field] = np.zeros(sizes[-1], dtype='f8')

    for i, hdf5_fname in enumerate(hdf5_fnames):
        f = h5py.File(hdf5_fname, "r")
        for field in columns:
            col = columns[field]
            dataset = os.path.join(hdf5_path, col)
            try:
                if row_mask is not None:
                    data[field][sizes[i]:sizes[i+1]] = f[dataset][masks[i]]
                else:
                    data[field][sizes[i]:sizes[i+1]] = f[dataset][:]
            except KeyError:
                raise RuntimeError(hdf5_fname, " should contain an dataset ", col,
                                   " under the path ", hdf5_path, ".")
        f.close()
    fh.write_file(fh_fname, data)


def write_hdf5_file(input_fh_fname, output_fh_fname, output_hdf5_fname,
                    hdf5_path, input_columns, output_columns):
    """
    Converts output fh file to hdf5, keeping only the specified columns
    """
    output_data = fh.read_file(output_fh_fname, list(output_columns))
    input_data = fh.read_file(input_fh_fname, list(input_columns))

    with h5py.File(output_hdf5_fname, "w") as f:
        for col in input_columns:
            out_col = input_columns[col]
            f[os.path.join(hdf5_path, out_col)] = input_data[col]
        for col in output_columns:
            out_col = output_columns[col]
            f[os.path.join(hdf5_path, col)] = output_data[col]


def write_fh_file(input_fh_fname, output_fh_fname, input_columns,
                  output_columns):

    def get_cols(fname):
        cols = glob(os.path.join(fname, "*.nd*"))
        cols = [col.split("/")[-1].split(".nd")[0] for col in cols]
        return cols

    def remove_column(fname, col):
        glob_expr = os.path.join(fname, col + ".nd.*")
        subprocess.run(["rm -f {}".format(glob_expr)], shell=True, check=True)

    def move_column(input_fname, in_col, dest_fname, out_col):
        full_input_fname = glob(os.path.join(input_fname, in_col + ".nd.*"))[0]
        extension = ".nd." + full_input_fname.split(".nd.")[1]
        full_dest_fname = os.path.join(dest_fname, out_col + extension)
        subprocess.run(["mv {} {}".format(full_input_fname, full_dest_fname)],
                       shell=True, check=True)

    def rename_column(fname, in_col, out_col):
        move_column(fname, in_col, fname, out_col)

    input_file_cols = get_cols(input_fh_fname)
    output_file_cols = get_cols(output_fh_fname)

    for col in output_file_cols:
        if col not in output_columns:
            remove_column(output_fh_fname, col)

    for col in output_columns:
        dest_col = output_columns[col]
        rename_column(output_fh_fname, col, dest_col)

    for col in input_columns:
        dest_col = input_columns[col]
        move_column(input_fh_fname, col, output_fh_fname, dest_col)


def read_config(config, section, key, default=None):
    "Reads (section, key) pair from config file, setting to default if not available."
    try:
        ret = config[section][key]
    except KeyError:
        ret = default
    return ret


if __name__ == "__main__":
    #try:
    #    config_file = sys.argv[1]
    #except IndexError:
    #    print("Usage: recon.py <YAML configuration file>")
    #    exit()
    parser = argparse.ArgumentParser(prog='recon.py', description="",)
    parser.add_argument('config_file', metavar='<YAML configuration file>', action='store', type=str,
                        help='yaml file holding configuration options')
    parser.add_argument('--data_path', metavar='<fits, hdf5 or fh file>', nargs='+',
                        help='path(s) to input data file (overrides configuration file)')
    parser.add_argument('--randoms1_path', nargs='+',
                        metavar='<fits, hdf5 or fh file>',
                        help='path(s) to first set of randoms (idem)')
    parser.add_argument('--randoms2_path', nargs='+',
                        metavar='<fits, hdf5 or fh file>',
                        help='path(s) to second set of randoms (idem)')
    parser.add_argument('--output_data_path', nargs='?', metavar='<fits, hdf5 or fh file>',
                        help='path to output data file (idem)')
    parser.add_argument('--output_randoms_path', metavar='<fits, hdf5 or fh file>', nargs='?',
                        help='path to output randoms file (idem)')

    args = parser.parse_args()
    config_file = args.config_file

    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    data = args.data_path
    randoms1 = args.randoms1_path
    randoms2 = args.randoms2_path
    output_data_file = args.output_data_path
    output_randoms_file = args.output_randoms_path


    if (data is None) and (randoms1 is None) and (randoms2 is None):
        input_path = config["input"]["input_path"]
    if data is None:
        data = config["input"]["data_fname"]
        if not isinstance(data, list):
            data = [data]
        data_files = list(map(lambda x: os.path.join(input_path, x), data))
    else:
        data_files = data

    if randoms1 is None:
        randoms1 = config["input"]["randoms1_fname"]
        if not isinstance(randoms1, list):
            randoms1 = [randoms1]
        randoms1_files = list(map(lambda x: os.path.join(input_path, x), randoms1))
    else:
        randoms1_files = randoms1

    if randoms2 is None:
        randoms2 = config["input"]["randoms2_fname"]
        if not isinstance(randoms2, list):
            randoms2 = [randoms2]
        randoms2_files = list(map(lambda x: os.path.join(input_path, x), randoms2))
    else:
        randoms2_files = randoms2

    hdf_preffix = read_config(config, "input", "hdf5_preffix", "/")
    ra_col = read_config(config, "input", "ra_col", "ra")
    dec_col = read_config(config, "input", "dec_col", "dec")
    redshift_col = read_config(config, "input", "redshift_col", "z")
    ra_col_randoms = read_config(config, "input", "ra_col_randoms", ra_col)
    dec_col_randoms = read_config(config, "input", "dec_col_randoms", dec_col)
    redshift_col_randoms = read_config(config, "input", "redshift_col_randoms", redshift_col)

    input_columns = {'ra': ra_col, 'dec': dec_col, 'z': redshift_col}
    input_columns_randoms = {'ra': ra_col_randoms, 'dec': dec_col_randoms, 'z': redshift_col_randoms}
    row_mask = read_config(config, "input", "row_mask", "")
    row_mask_randoms = read_config(config, "input", "row_mask_randoms", row_mask)
    if row_mask == "":
        row_mask = None
    if row_mask_randoms == "":
        row_mask_randoms = None


    if output_data_file is not None:
        tentative_output_path = '/'.join(output_data_file.split("/")[:-1])
    output_path = read_config(config, "output", "output_path", tentative_output_path)
    if output_data_file is None:
        output_data_file = os.path.join(output_path,
                                        config["output"]["output_data_fname"])
    if output_randoms_file is None:
        output_randoms_file = os.path.join(output_path,
                                           config["output"]["output_randoms_fname"])
    recon_suffix = read_config(config, "output", "columns_suffix", "_recon").lower()

    output_columns_str = read_config(config, "output", "output_coordinates", "both")
    output_columns = {}
    if (output_columns_str == "both") or (output_columns_str == "spherical"):
        output_columns['ra'] = (ra_col+recon_suffix).lower()
        output_columns['dec'] = (dec_col+recon_suffix).lower()
        output_columns['z'] = (redshift_col+recon_suffix).lower()
    if (output_columns_str == "both") or (output_columns_str == "cartesian"):
        output_columns['x_coord'] = 'x_coord'+recon_suffix
        output_columns['y_coord'] = 'y_coord'+recon_suffix
        output_columns['z_coord'] = 'z_coord'+recon_suffix

    keep_cols = read_config(config, "output", "keep_cols", "")
    keep_cols_randoms = read_config(config, "output", "keep_cols_randoms", keep_cols)
    keep_cols_dict = {}
    keep_cols_dict_randoms = {}
    if keep_cols != "":
        for col in keep_cols:
            input_columns[col.lower()] = col
            keep_cols_dict[col.lower()] = col.lower()
    if keep_cols != "":
        for col in keep_cols_randoms:
            input_columns_randoms[col.lower()] = col
            keep_cols_dict_randoms[col.lower()] = col.lower()



    b = config["reconstruction_config"]["bias"]
    f = config["reconstruction_config"]["f"]
    Rf = config["reconstruction_config"]["R_smooth"]
    Om = config["reconstruction_config"]["Omega_m"]
    rsd_convention = config["reconstruction_config"]["RSD_convention"]

    data_fmt = os.path.splitext(data_files[0])[-1][1:]
    randoms1_fmt = os.path.splitext(randoms1_files[0])[-1][1:]
    randoms2_fmt = os.path.splitext(randoms2_files[0])[-1][1:]
    output_data_fmt = os.path.splitext(output_data_file)[-1][1:]
    output_randoms_fmt = os.path.splitext(output_randoms_file)[-1][1:]

    hdf_exts = ['hdf', 'h4', 'hdf4', 'he2', 'h5', 'hdf5', 'he5', 'h5py']

    with tempfile.TemporaryDirectory(dir=output_path) as tmpdirname:
        fh_data_file = os.path.join(tmpdirname, "data.fh")
        fh_randoms1_file = os.path.join(tmpdirname, "randoms1.fh")
        fh_randoms2_file = os.path.join(tmpdirname, "randoms2.fh")
        fh_output_data_file = os.path.join(tmpdirname, "data_out.fh")
        fh_output_randoms_file = os.path.join(tmpdirname, "randoms_out.fh")

        if output_data_fmt == "fh":
            fh_output_data_file = output_data_file
        if output_randoms_fmt == "fh":
            fh_output_randoms_file = output_randoms_file

        if data_fmt == "fits":
            read_fits_files(fits_fnames=data_files, fh_fname=fh_data_file,
                            columns=input_columns, row_mask=row_mask)
        elif data_fmt in hdf_exts:
            read_hdf5_files(hdf5_fnames=data_files, hdf5_path=hdf_preffix,
                            fh_fname=fh_data_file, columns=input_columns,
                            row_mask=row_mask)
        elif data_fmt == "fh":
            fh_data_file = data_files[0]
            fh_output_data_file = output_data_file
        else:
            raise RuntimeError("Unsupported format ", data_fmt, " for input data.")

        if randoms1_fmt == "fits":
            read_fits_files(fits_fnames=randoms1_files, fh_fname=fh_randoms1_file,
                            columns=input_columns_randoms, row_mask=row_mask_randoms)
        elif randoms1_fmt in hdf_exts:
            read_hdf5_files(hdf5_fnames=randoms1_files, hdf5_path=hdf_preffix,
                            fh_fname=fh_randoms1_file, columns=input_columns_randoms,
                            row_mask=row_mask_randoms)
        elif randoms1_fmt == "fh":
            fh_randoms1_file = randoms1_files[0]
        else:
            raise RuntimeError("Unsupported format ", randoms1_fmt, " for the randoms.")

        if randoms2_fmt == "fits":
            read_fits_files(fits_fnames=randoms2_files, fh_fname=fh_randoms2_file,
                            columns=input_columns_randoms, row_mask=row_mask_randoms)
        elif randoms2_fmt in hdf_exts:
            read_hdf5_files(hdf5_fnames=randoms2_files, hdf5_path=hdf_preffix,
                            fh_fname=fh_randoms2_file, columns=input_columns_randoms,
                            row_mask=row_mask_randoms)
        elif randoms2_fmt == "fh":
            fh_randoms2_file = randoms2_files[0]
            fh_output_randoms_file = output_randoms_file
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
        #      output_data_file.encode(), output_randoms_file.encode(),
        #      c_float(b), c_float(f), c_float(Rf), c_float(Om), reciso.encode())

        cmd = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), "recon")
        args = "{} {} {} {} {} {} {} {} {} {}".format(fh_data_file, fh_randoms1_file,
                                                      fh_randoms2_file, fh_output_data_file,
                                                      fh_output_randoms_file, b, f, Rf, Om,
                                                      rsd_convention)
        subprocess.run([cmd + " " + args], shell=True, check=True)


        if output_data_fmt == "fits":
            write_fits_file(input_fh_fname=fh_data_file,
                            output_fh_fname=fh_output_data_file,
                            output_fits_fname=output_data_file,
                            input_columns=keep_cols_dict,
                            output_columns=output_columns)
        elif output_data_fmt in hdf_exts:
            write_hdf5_file(input_fh_fname=fh_data_file,
                            output_fh_fname=fh_output_data_file,
                            output_hdf5_fname=output_data_file,
                            hdf5_path=hdf_preffix,
                            input_columns=keep_cols_dict,
                            output_columns=output_columns)
        elif output_data_fmt == "fh":
            write_fh_file(input_fh_fname=fh_data_file,
                          output_fh_fname=output_data_file,
                          input_columns=keep_cols_dict,
                          output_columns=output_columns)
        else:
            raise RuntimeError("Unsupported format ", output_data_fmt, " for output data.")

        if output_randoms_fmt == "fits":
            write_fits_file(input_fh_fname=fh_randoms2_file,
                            output_fh_fname=fh_output_randoms_file,
                            output_fits_fname=output_randoms_file,
                            input_columns=keep_cols_dict_randoms,
                            output_columns=output_columns)
        elif output_randoms_fmt in hdf_exts:
            write_hdf5_file(input_fh_fname=fh_randoms2_file,
                            output_fh_fname=fh_output_randoms_file,
                            output_hdf5_fname=output_randoms_file,
                            hdf5_path=hdf_preffix,
                            input_columns=keep_cols_dict_randoms,
                            output_columns=output_columns)
        elif output_randoms_fmt == "fh":
            write_fh_file(input_fh_fname=fh_randoms2_file,
                          output_fh_fname=output_randoms_file,
                          input_columns=keep_cols_dict_randoms,
                          output_columns=output_columns)
        else:
            raise RuntimeError("Unsupported format ", output_randoms_fmt, " for output randoms.")

