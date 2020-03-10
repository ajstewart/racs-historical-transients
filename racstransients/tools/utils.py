#!/usr/bin/env python

import os
import subprocess
import logging
import glob
import sys
import pandas as pd
import whichcraft

try:
    import colorlog
except ImportError:
    pass


logger = logging.getLogger(__name__)

def setup_logging(logfilename,level,use_colorlog):
    levels={'CRITICAL' : logging.CRITICAL,
        'ERROR' : logging.ERROR,
        'WARNING' : logging.WARNING,
        'INFO' : logging.INFO,
        'DEBUG' : logging.DEBUG
    }

    root_logger = logging.getLogger()
    root_logger.setLevel(levels[level])

    # create a file handler
    handler = logging.FileHandler('{}.log'.format(logfilename), mode='w')
    handler.setLevel(levels[level])

    # create a logging format
    logformat='[%(asctime)s] - %(name)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(logformat, datefmt="%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)

    # add the handlers to the logger
    root_logger.addHandler(handler)
    if use_colorlog:
        formatter = colorlog.ColoredFormatter(
            # "%(log_color)s%(levelname)-8s%(reset)s %(blue)s%(message)s",
            "%(log_color)s[%(asctime)s] - %(name)s - %(levelname)s - %(blue)s%(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            reset=True,
            log_colors={
            'DEBUG':    'cyan',
            'INFO':     'green',
            'WARNING':  'yellow',
            'ERROR':    'red',
            'CRITICAL': 'red,bg_white',
            },
            secondary_log_colors={},
            style='%'
        )
    term=logging.StreamHandler()
    term.setLevel(levels[level])
    term.setFormatter(formatter)
    root_logger.addHandler(term)

    logger.info("Job Started")

def createdir(name, clobber=False, cont=False):
    if os.path.isdir(name):
        if cont==True:
            logger.warning("Continue = True!")
            logger.warning("Working on directory without refreshing!")
            return True
        if clobber==True:
            logger.warning("Clobber = True!")
            logger.warning("Previous output directory '{}' will be overwritten".format(name))
            subprocess.call(["rm", "-rf", name])
        else:
            logger.error("Previous output directory '{0}' already exists and clobber set to {1}".format(name, clobber))
            return False
    else:
        if cont==True:
            logger.error("Directory chosen to continue run {} does not exist!".format(name))
            sys.exit()
    os.mkdir(name)
    return True

def createdir_print(name, clobber=False, cont=False):
    if os.path.isdir(name):
        if cont==True:
            print("Continue = True!")
            print("Working on directory without refreshing!")
            return True
        if clobber==True:
            print("Clobber = True!")
            print("Previous output directory '{}' will be overwritten".format(name))
            subprocess.call(["rm", "-rf", name])
        else:
            print("Previous output directory '{0}' already exists and clobber set to {1}".format(name, clobber))
            return False
    else:
        if cont==True:
            print("Directory chosen to continue run {} does not exist!".format(name))
            sys.exit()
    os.mkdir(name)
    return True

def create_sub_dirs(dirs, cont=False):
    if not cont:
        for i in dirs:
            os.mkdir(i)
    return

def checkdir(files):
    if isinstance(files,list):
        for i in files:
            if not os.path.isdir(i):
                logger.error("{} cannot be found".format(i))
                return False
    else:
        if not os.path.isdir(files):
            logger.error("{} cannot be found".format(files))
            return False
    return True


def checkfile(files):
    if isinstance(files,list):
        for i in files:
            if not os.path.isfile(i):
                logger.error("{} cannot be found".format(i))
                return False
    else:
        if not os.path.isfile(files):
            logger.error("{} cannot be found".format(files))
            return False
    return True

def checkfile_print(files):
    if isinstance(files,list):
        for i in files:
            if not os.path.isfile(i):
                print("{} cannot be found".format(i))
                return False
    else:
        if not os.path.isfile(files):
            print("{} cannot be found".format(files))
            return False
    return True

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    # from shutil import which

    return whichcraft.which(name) is not None

def pybdsf2aegean(pybdsf_output, new_ouput):
    pybdsf=pd.read_csv(pybdsf_output, skiprows=5, sep=", ", engine="python")

    pybdsf_new=pybdsf.filter(["Isl_id", "Source_id", "Isl_rms", "RA", "E_RA", "DEC", "E_DEC", "Peak_flux", "E_Peak_flux", "Total_flux", "E_Total_flux", "Maj", "E_Maj",
        "Min", "E_Min", "PA", "E_PA", "Maj_img_plane", "Min_img_plane", "PA_img_plane"])

    pybdsf_new["Maj"]=pybdsf_new["Maj"]*3600.
    pybdsf_new["E_Maj"]=pybdsf_new["E_Maj"]*3600.
    pybdsf_new["Min"]=pybdsf_new["Min"]*3600.
    pybdsf_new["E_Min"]=pybdsf_new["E_Min"]*3600.
    pybdsf_new["PA"]=pybdsf_new["PA"]*3600.
    pybdsf_new["E_PA"]=pybdsf_new["E_PA"]*3600.
    pybdsf_new["Maj_img_plane"]=pybdsf_new["Maj_img_plane"]*3600.
    pybdsf_new["Min_img_plane"]=pybdsf_new["Min_img_plane"]*3600.
    pybdsf_new["PA_img_plane"]=pybdsf_new["PA_img_plane"]*3600.

    pybdsf_new.columns=["island","source","local_rms","ra","err_ra","dec","err_dec","peak_flux","err_peak_flux","int_flux","err_int_flux","a","err_a","b","err_b","pa","err_pa",
    "psf_a", "psf_b", "psf_pa"]

    pybdsf_new["flags"]=0

    pybdsf_new.to_csv(new_ouput, index=False, sep=",")

def selavy2aegean(selavy_output, new_ouput):
    with open(selavy_output, "r") as f:
        lines=f.readlines()

    columns=lines[0].split()[1:-1]
    data=[i.split() for i in lines[2:]]

    selavy=pd.DataFrame(data, columns=columns)

#   island_id    component_id component_name ra_hms_cont dec_dms_cont ra_deg_cont dec_deg_cont   ra_err  dec_err  freq  flux_peak flux_peak_err flux_int flux_int_err maj_axis min_axis pos_ang maj_axis_err min_axis_err pos_ang_err maj_axis_deconv min_axis_deconv pos_ang_deconv maj_axis_deconv_err min_axis_deconv_err pos_ang_deconv_err chi_squared_fit rms_fit_gauss spectral_index spectral_curvature spectral_index_err spectral_curvature_err  rms_image has_siblings fit_is_estimate spectral_index_from_TT flag_c4

    selavy_new=selavy.filter(["island_id", "component_id", "rms_image", "ra_deg_cont", "ra_err", "dec_deg_cont", "dec_err", "flux_peak", "flux_peak_err", "flux_int", "flux_int_err", "maj_axis",
        "maj_axis_err", "min_axis", "min_axis_err", "pos_ang", "pos_ang_err", "maj_axis_deconv", "min_axis_deconv", "pos_ang_deconv"])

    selavy_new["flux_peak"]=selavy_new["flux_peak"].astype(float)/1.e3
    selavy_new["flux_peak_err"]=selavy_new["flux_peak_err"].astype(float)/1.e3
    selavy_new["flux_int_err"]=selavy_new["flux_int_err"].astype(float)/1.e3
    selavy_new["flux_int"]=selavy_new["flux_int"].astype(float)/1.e3
    selavy_new["rms_image"]=selavy_new["rms_image"].astype(float)/1.e3

    selavy_new.columns=["island","source","local_rms","ra","err_ra","dec","err_dec","peak_flux","err_peak_flux","int_flux","err_int_flux","a","err_a","b","err_b","pa","err_pa",
        "psf_a", "psf_b", "psf_pa"]

    selavy_new["flags"]=0

    selavy_new.to_csv(new_ouput, index=False, sep=",")
