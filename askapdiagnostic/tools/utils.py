#!/usr/bin/env python

import os
import subprocess
import logging
import glob
import sys

try:
    import colorlog
except ImportError:
    pass
    

logger = logging.getLogger(__name__)

def setup_logging(logfilename,use_colorlog):
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # create a file handler
    handler = logging.FileHandler('{}.log'.format(logfilename), mode='w')
    handler.setLevel(logging.INFO)

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
    term.setLevel(logging.INFO)
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
            print "Continue = True!"
            print "Working on directory without refreshing!"
            return True
        if clobber==True:
            print "Clobber = True!"
            print "Previous output directory '{}' will be overwritten".format(name)
            subprocess.call(["rm", "-rf", name])
        else:
            print "Previous output directory '{0}' already exists and clobber set to {1}".format(name, clobber)
            return False
    else:
        if cont==True:
            print "Directory chosen to continue run {} does not exist!".format(name)
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
                print "{} cannot be found".format(i)
                return False
    else:
        if not os.path.isfile(files):
            print "{} cannot be found".format(files)
            return False
    return True