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


fits_exts = ['fits']
hdf_exts = ['hdf', 'h4', 'hdf4', 'he2', 'h5', 'hdf5', 'he5', 'h5py']
ascii_exts = ['dat', 'txt', 'tsv', 'csv', 'xyz', 'xyzw']


def read_fits_files(fits_fnames, fh_fname, columns, row_mask, **kwargs):
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
                    input_columns, output_columns, **kwargs):
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

def read_hdf5_files(hdf5_fnames, hdf5_path, fh_fname, columns, row_mask, **kwargs):
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
            sizes[i+1] = sizes[i] + f[os.path.join(hdf5_path, columns[columns.keys()[0]])].len()
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
                    hdf5_path, input_columns, output_columns, **kwargs):
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

def read_ascii_files(ascii_fnames, fh_fname, columns, row_mask, **kwargs):
    """
    Reads relevant columns from ascii file and converts to a filehandler file
    to be passed to the C code
    """
    data = {}

    data = np.concatenate([np.loadtxt(ascii_fname) for ascii_fname in ascii_fnames])

    # Processing row_mask string to apply row selection
    if row_mask is not None:
        while (row_mask.find("${") != -1):
            col_label = row_mask[row_mask.find("${"):].split("}")[0][2:]
            sub_expression = "data_arr[:,{}]".format(int(col_label))
            row_mask = row_mask.replace("${"+col_label+"}", sub_expression)

    data_cols = {field: [] for field in columns}
    for i, ascii_fname in enumerate(ascii_fnames):
        data_arr = np.loadtxt(ascii_fname)
        if row_mask is not None:
            try:
                mask = eval(row_mask)
            except:
                raise RuntimeError("Invalid row selection specified.")
            data_arr = data_arr[mask]
        for field in columns:
            col = columns[field]
            try:
                col_data = data_arr[:,col]
            except KeyError:
                raise RuntimeError(ascii_fname, " should contain an column at key ", col)
            data_cols[field].append(col_data)
    data = {field: np.concatenate(data_cols[field]).astype('f8') for field in columns}
    fh.write_file(fh_fname, data)

def write_ascii_file(input_fh_fname, output_fh_fname, output_ascii_fname,
                     input_columns, output_columns, **kwargs):
    """
    Converts output fh file to ascii, keeping only the specified columns
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
    columns = [output_columns[col] for col in output_columns]
    columns.extend([input_columns[col] for col in input_columns])
    tab = Table(data)[columns]
    tab.write(output_ascii_fname, format='ascii.fast_commented_header')

def read_input_files(fnames, **kwargs):
    data_fmts = [os.path.splitext(fname)[-1][1:] for fname in fnames]
    for fmt in data_fmts:
        if fmt != data_fmts[0]:
            raise RuntimeError("Input formats should be uniform when multiple files provided.")
    data_fmt = data_fmts[0]

    if data_fmt in fits_exts:
        read_fits_files(fits_fnames=fnames, **kwargs)
    elif data_fmt in hdf_exts:
        read_hdf5_files(hdf5_fnames=fnames, **kwargs)
    elif data_fmt in ascii_exts:
        read_ascii_files(ascii_fnames=fnames, **kwargs)
    else:
        raise RuntimeError("Unsupported format ", data_fmt, " for input file ", fnames[0], ".")

def write_output_file(output_fname, **kwargs):
    output_data_fmt = os.path.splitext(output_fname)[-1][1:]

    if output_data_fmt in fits_exts:
        write_fits_file(output_fits_fname=output_fname, **kwargs)
    elif output_data_fmt in hdf_exts:
        write_hdf5_file(output_hdf5_fname=output_fname, **kwargs)
    elif output_data_fmt in ascii_exts:
        write_ascii_file(output_ascii_fname=output_fname, **kwargs)
    else:
        raise RuntimeError("Unsupported format ", output_data_fmt, " for output file ", output_fname, ".")

def load_config(parser):
    # namespace to return
    ns = argparse.Namespace

    # parse input and load YAML file
    args = parser.parse_args()
    with open(args.config_file, "r") as f:
        config = yaml.safe_load(f)

    # Small function to read from config, optionally setting defaults and overwrites
    def read_config(config, section, key, default=None, overwrite=None):
        try:
            ret = config[section][key]
        except KeyError:
            ret = default
        if overwrite is not None:
            ret = overwrite
        return ret

    # Small function to handle joining base_paths and filenames from config files
    # or load directly from cli
    def load_fname_paths(section, base_path_key, key, overwrite=None):
        if overwrite is not None:
            paths = overwrite
        else:
            input_path  = read_config(config, section, base_path_key, default="/")
            fnames = read_config(config, section, key)
            if not isinstance(fnames, list):
                fnames = [fnames]
            paths = [os.path.join(input_path, file_) for file_ in fnames]

        if paths is None:
            raise RuntimeError("If input/{} is not specified in the configuration file, "
                               "it should be provided as argument "
                               "to the application call.".format(config_key))
        return paths

    # Reads input/output files
    ns.data_files          = load_fname_paths("input", "input_path", "data_fname",
                                              overwrite=args.data_path)
    ns.randoms1_files      = load_fname_paths("input", "input_path", "randoms1_fname",
                                              overwrite=args.randoms1_path)
    ns.randoms2_files      = load_fname_paths("input", "input_path", "randoms2_fname",
                                              overwrite=args.randoms2_path)
    ns.output_data_file    = load_fname_paths("output", "output_path", "output_data_fname",
                                              overwrite=args.output_data_path)[0]
    ns.output_randoms_file = load_fname_paths("output", "output_path", "output_randoms_fname",
                                              overwrite=args.output_randoms_path)[0]

    ns.output_path = read_config(config, "output", "output_path",
                                 default=os.path.dirname(ns.output_data_file))

    # hdf_preffix of dataset if hdf5 format is used
    ns.hdf_preffix = read_config(config, "input", "hdf5_preffix", "/")

    # Reads input columns
    ra_col               = read_config(config, "input", "ra_col", "ra")
    dec_col              = read_config(config, "input", "dec_col", "dec")
    redshift_col         = read_config(config, "input", "redshift_col", "z")
    ra_col_randoms       = read_config(config, "input", "ra_col_randoms", ra_col)
    dec_col_randoms      = read_config(config, "input", "dec_col_randoms", dec_col)
    redshift_col_randoms = read_config(config, "input", "redshift_col_randoms", redshift_col)
    x_col                = read_config(config, "input", "x_col", "x_coord")
    y_col                = read_config(config, "input", "y_col", "y_coord")
    z_col                = read_config(config, "input", "z_col", "z_coord")
    x_col_randoms        = read_config(config, "input", "x_col_randoms", x_col)
    y_col_randoms        = read_config(config, "input", "y_col_randoms", y_col)
    z_col_randoms        = read_config(config, "input", "z_col_ranndoms", z_col)

    ns.is_cartesian = read_config(config, "input", "input_coordinates", "spherical").lower() == "cartesian"
    if not ns.is_cartesian:
        ns.input_columns = {'ra': ra_col,
                            'dec': dec_col,
                            'z': redshift_col}
        ns.input_columns_randoms = {'ra': ra_col_randoms,
                                    'dec': dec_col_randoms,
                                    'z': redshift_col_randoms}
    else:
        ns.input_columns = {'x_coord': x_col,
                            'y_coord': y_col,
                            'z_coord': z_col}
        ns.input_columns_randoms = {'x_coord': x_col_randoms,
                                    'y_coord': y_col_randoms,
                                    'z_coord': z_col_randoms}

    # Reads row mask
    ns.row_mask         = read_config(config, "input", "row_mask", "")
    ns.row_mask_randoms = read_config(config, "input", "row_mask_randoms", ns.row_mask)
    if ns.row_mask == "":
        ns.row_mask = None
    if ns.row_mask_randoms == "":
        ns.row_mask_randoms = None

    # Reads output columns
    recon_suffix       = read_config(config, "output", "columns_suffix", "_recon").lower()
    output_columns_str = read_config(config, "output", "output_coordinates", "both")
    ns.output_columns = {}
    if (output_columns_str == "both") or (output_columns_str == "spherical"):
        if ns.is_cartesian:
            raise RuntimeError("Not implemented spherical output for cartesian input.")
        ns.output_columns['ra']  = (ra_col+recon_suffix).lower()
        ns.output_columns['dec'] = (dec_col+recon_suffix).lower()
        ns.output_columns['z']   = (redshift_col+recon_suffix).lower()
    if (output_columns_str == "both") or (output_columns_str == "cartesian"):
        ns.output_columns['x_coord'] = 'x_coord'+recon_suffix
        ns.output_columns['y_coord'] = 'y_coord'+recon_suffix
        ns.output_columns['z_coord'] = 'z_coord'+recon_suffix

    # Reads columns to keep from input and add to both input and output columns
    keep_cols         = read_config(config, "output", "keep_cols", "")
    keep_cols_randoms = read_config(config, "output", "keep_cols_randoms", keep_cols)
    ns.keep_cols_dict = {}
    ns.keep_cols_dict_randoms = {}
    if keep_cols != "":
        for col in keep_cols:
            ns.input_columns[str(col).lower()] = col
            ns.keep_cols_dict[str(col).lower()] = str(col).lower()
    if keep_cols != "":
        for col in keep_cols_randoms:
            ns.input_columns_randoms[str(col).lower()] = col
            ns.keep_cols_dict_randoms[str(col).lower()] = str(col).lower()

    # Reads parameters to pass to recon code
    ns.b              = read_config(config, "reconstruction_config", "bias",
                                    overwrite=args.bias)
    ns.f              = read_config(config, "reconstruction_config", "f",
                                    overwrite=args.f)
    ns.Rf             = read_config(config, "reconstruction_config", "R_smooth",
                                    overwrite=args.R_smooth)
    ns.Om             = read_config(config, "reconstruction_config", "Omega_m",
                                    overwrite=args.Omega_m)
    ns.Ngrid          = read_config(config, "reconstruction_config", "Ngrid",
                                    default=512, overwrite=args.Ngrid)
    ns.rsd_convention = read_config(config, "reconstruction_config", "RSD_convention",
                                    overwrite=args.RSD_convention)

    return ns


if __name__ == "__main__":
    #CLI parsing
    parser = argparse.ArgumentParser(prog='recon.py', description="",
                                     epilog='All optional parameters, if specified, overwrite corresponding values set at the specified configuration file. Accepted input file formats are ascii, fits and hdf5.')
    parser.add_argument('config_file', metavar='<YAML configuration file>', action='store', type=str,
                        help='yaml file holding configuration options')
    parser.add_argument('--data_path', metavar='<fname>', nargs='+',
                        help='path(s) to input data file')
    parser.add_argument('--randoms1_path', nargs='+',
                        metavar='<fname>',
                        help='path(s) to first set of randoms')
    parser.add_argument('--randoms2_path', nargs='+',
                        metavar='<fname>',
                        help='path(s) to second set of randoms')
    parser.add_argument('--output_data_path', nargs=1, metavar='<fname>',
                        help='path to output data file')
    parser.add_argument('--output_randoms_path', metavar='<fname>', nargs=1,
                        help='path to output randoms file')
    parser.add_argument('--bias', metavar='<bias>', help='value of bias')
    parser.add_argument('--f', metavar='<f>', help='value of growth factor')
    parser.add_argument('--R_smooth', metavar='<Rs>', help='smoothing scale to use (in Mpc/h)')
    parser.add_argument('--Omega_m', metavar='<Om>', help='Omega m value to assume')
    parser.add_argument('--Ngrid', metavar='<Ng>', help='Grid size')
    parser.add_argument('--RSD_convention', metavar='<RSD_convention>',
                        choices=['RecIso', 'RecSym'], help='RSD convention (RecSym or RecIso)')

    #Process configuration using both CLI and provided YAML config file
    config = load_config(parser)

    with tempfile.TemporaryDirectory(dir=config.output_path) as tmpdirname:
        fh_data_file = os.path.join(tmpdirname, "data.fh")
        fh_randoms1_file = os.path.join(tmpdirname, "randoms1.fh")
        fh_randoms2_file = os.path.join(tmpdirname, "randoms2.fh")
        fh_output_data_file = os.path.join(tmpdirname, "data_out.fh")
        fh_output_randoms_file = os.path.join(tmpdirname, "randoms_out.fh")

        # Read input files and convert to fh format
        read_input_files(fnames=config.data_files, fh_fname=fh_data_file,
                         columns=config.input_columns, row_mask=config.row_mask,
                         hdf5_path=config.hdf_preffix)
        read_input_files(fnames=config.randoms1_files, fh_fname=fh_randoms1_file,
                         columns=config.input_columns_randoms, row_mask=config.row_mask_randoms,
                         hdf5_path=config.hdf_preffix)
        read_input_files(fnames=config.randoms2_files, fh_fname=fh_randoms2_file,
                         columns=config.input_columns_randoms, row_mask=config.row_mask_randoms,
                         hdf5_path=config.hdf_preffix)

        # Call C subprocess
        cmd = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), "recon")
        args = " ".join(str(x) for x in [fh_data_file, fh_randoms1_file, fh_randoms2_file,
                                         fh_output_data_file, fh_output_randoms_file,
                                         config.b, config.f, config.Rf, config.Om, config.Ngrid,
                                         config.rsd_convention, int(config.is_cartesian)])
        subprocess.run([cmd + " " + args], shell=True, check=True)


        # Convert back from fh and write to definitive output files
        write_output_file(input_fh_fname=fh_data_file,
                          output_fh_fname=fh_output_data_file,
                          output_fname=config.output_data_file,
                          hdf5_path=config.hdf_preffix,
                          input_columns=config.keep_cols_dict,
                          output_columns=config.output_columns)

        write_output_file(input_fh_fname=fh_randoms2_file,
                          output_fh_fname=fh_output_randoms_file,
                          output_fname=config.output_randoms_file,
                          hdf5_path=config.hdf_preffix,
                          input_columns=config.keep_cols_dict_randoms,
                          output_columns=config.output_columns)

