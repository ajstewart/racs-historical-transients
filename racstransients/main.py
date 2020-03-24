#!/usr/bin/env python

import argparse
import logging
from racstransients.tools import utils
from racstransients.tools.fitsimage import Askapimage
from racstransients.tools.crossmatch import crossmatch
from racstransients.tools.catalog import Catalog
from racstransients.plotting import plots
import multiprocessing
import datetime
import os
import pandas as pd
import sys
import subprocess
import configparser
import getpass
# from multiprocessing import Pool


#Steps are:
# 1. Read ASKAP image and get centre and size.
# 2. Perform source finding on image / or read already performed crossmatching.
# 3. Fetch SUMSS from Vizier / or read in SUMSS already performed.
# 4. Mask SUMSS to leave only sources that should be in the image.
# 5. Perform cross-matching using above outputs.
# 6. Produce plots and perform other analysis.

def exit(logger):
    logger.info("Exiting...")
    sys.exit()

def source_finding(askapimg, sf, logger, options=None, save_diag_images=False):
    if options!=None:
        config = configparser.SafeConfigParser()
        config.optionxform = str
        config.read([options])
    if sf=="aegean":
        if options==None:
            aegean_sf_options={
                "cores":1,
                "maxsummits":5,
                "seedclip":5,
                "floodclip":4,
                "nocov":""
                }
        else:
            aegean_sf_options = dict(config.items("aegean"))
            todel=[]
            for i in aegean_sf_options:
                if aegean_sf_options[i].lower()=="true":
                    aegean_sf_options[i]=""
                #account for if someone puts false
                elif aegean_sf_options[i].lower()=="false":
                    todel.append(i)
            if len(todel)>0:
                for d in todel:
                    del aegean_sf_options[d]

        askapimg.find_sources(sf=sf, options=aegean_sf_options)
        if save_diag_images:
            aegean_sf_options["save"]=""
            askapimg.find_sources(options=aegean_sf_options)
        if "floodclip" in aegean_sf_options:
            aegean_fill_sigma=aegean_sf_options["floodclip"]
        else:
            aegean_sigma=99
        if "seedclip" in aegean_sf_options:
            aegean_det_sigma=aegean_sf_options["seedclip"]
        else:
            aegean_det_sigma=99
        return [float(aegean_det_sigma), float(aegean_fill_sigma)]

    elif sf=="pybdsf":
        if options==None:
            pybdsf_sf_options={
                "frequency":askapimg.freq,
                }
        else:
            pybdsf_sf_options = dict(config.items("pybdsf"))

        askapimg.find_sources(sf=sf, options=pybdsf_sf_options)

        return [5.0, 3.0]

    elif sf=="selavy":
        base_options=theimg._get_default_selavy_options()
        if options!=None:
            selavy_sf_options = dict(config.items("selavy"))
            for i in base_options:
                if i not in selavy_sf_options:
                    selavy_sf_options[i]=base_options[i]

        askapimg.find_sources(sf=sf, options=selavy_sf_options)

        return [5.0, 3.0]

def convert_str_to_list(liststring):
    return liststring.split(",")

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def str2int(v):
    try:
        v = int(v)
        return v
    except:
        raise argparse.ArgumentTypeError('Integer value expected.')

def str2float(v):
    try:
        v = float(v)
        return v
    except:
        raise argparse.ArgumentTypeError('Float value expected.')

def str2upper(v):
    try:
        v = v.upper()
        return v
    except:
        raise argparse.ArgumentTypeError('String value expected.')

def get_catalogue(image, catalogue, boundary_value, logger):
    source_cat=image.imagename.replace(".fits", "_{}_comp.csv".format(catalogue.lower()))
    cat_df=image.get_catalogue(catalogue, boundary_value=boundary_value)
    if len(cat_df.index)<1:
        logger.error("{} catalog fetching seems to have failed.".format(catalogue))
        exit(logger)
    #Write the filtered catalog to disk
    if catalogue=="NVSS":
        image.write_nvss_sources()
    else:
        image.write_sumss_sources()
    logger.info("Fetching {} catalogue complete.".format(catalogue))
    logger.info("Written to disk as {}.".format(source_cat))
    return source_cat, cat_df

def determine_catalogues(image_dec, sumss_only, nvss_only, logger):
    #Check if NVSS or SUMSS is required (or both)
    if sumss_only and nvss_only:
        logger.critical("'sumss-only' and 'nvss-only' has been selected, please choose one!")
        sys.exit()
    if sumss_only:
        logger.warning("SUMSS only is selected! Will not use NVSS even if available.")
    elif nvss_only:
        logger.warning("NVSS only is selected! Will not use SUMSS even if available.")
    if image_dec <= -30.0:
        if nvss_only:
            sumss=False
            basecat="nvss"
        else:
            sumss=True
            basecat="sumss"
        if image_dec <= -43:
            nvss=False
        else:
            if sumss_only:
                nvss=False
            else:
                nvss=True
    else:
        if image_dec >= -27.0:
            sumss=False
        else:
            if not nvss_only:
                sumss=True
        if not sumss_only:
            basecat="nvss"
            nvss=True
        else:
            basecat="sumss"
            nvss=False

    if sumss==False and nvss==False:
        logger.error("No catalogues chosen for analysis (Dec = {}, SUMSS only is {})".format(image_dec, sumss_only))
        exit(logger)

    if sumss and nvss:
        dualmode=True
    else:
        dualmode=False

    if dualmode:
        logger.info("Both SUMSS and NVSS will be used.")
    logger.info("{} is the base catalog.".format(basecat.upper()))

    if dualmode and basecat=="sumss":
        matched_to_tag="SUMSS & NVSS"
    elif dualmode and basecat=="nvss":
        matched_to_tag="NVSS & SUMSS"
    elif sumss:
        matched_to_tag="SUMSS"
    else:
        matched_to_tag="NVSS"

    return sumss, nvss, dualmode, basecat, matched_to_tag


def main():
    launchtime=datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    logger = logging.getLogger(__name__)

    try:
        import colorlog
        use_colorlog=True
    except ImportError:
        use_colorlog=False

    conf_parser = argparse.ArgumentParser(
        # Turn off help, so we print all options in response to -h
            add_help=False
            )
    conf_parser.add_argument("-c", "--conf_file", help="Specify config file", metavar="FILE")
    args, remaining_argv = conf_parser.parse_known_args()
    defaults = {
        "output_tag":"",
        "log_level":"INFO",
        "nice":10,
        "clobber":False,
        "sumss_only":False,
        "nvss_only":False,
        "weight_crop":False,
        "weight_crop_value":0.04,
        "weight_crop_image":"weights.fits",
        "convolve":False,
        "convolved_image":"None",
        "convolved_non_conv_askap_csv":None,
        "convolved_non_conv_askap_islands_csv":None,
        "sourcefinder":"aegean",
        "frequency":99,
        "askap_csv":None,
        "askap_islands_csv":None,
        "sumss_csv":None,
        "nvss_csv":None,
        "askap_csv_format":"aegean",
        "remove_extended":False,
        "askap_ext_thresh":1.2,
        "sumss_ext_thresh":1.2,
        "nvss_ext_thresh":1.2,
        "use_all_fits":False,
        "write_ann":False,
        "produce_overlays":True,
        "boundary_value":"nan",
        "askap_flux_error":0.0,
        "diagnostic_max_separation":5.0,
        "transient_max_separation":45.0,
        "postage_stamps":False,
        "postage_stamp_selection":"all",
        "postage_stamp_ncores":int(multiprocessing.cpu_count()/2),
        "postage_stamp_radius":13.0,
        "postage_stamp_zscale_contrast":0.2,
        "sumss_mosaic_dir":"None",
        "nvss_mosaic_dir":"None",
        "aegean_settings_config":None,
        "pybdsf_settings_config":None,
        "selavy_settings_config":None,
        "transients":False,
        "transients_askap_snr_thresh":5.0,
        "transients_large_flux_ratio_thresh":3.0,
        "db_inject":True,
        "db_engine":"postgresql",
        "db_username":"postgres",
        "db_password":"postgres",
        "db_host":"localhost",
        "db_port":"5432",
        "db_database":"postgres",
        "db_tag":"RACS Analysis",
        "website_media_dir":"none"
        }

    if args.conf_file:
        if not utils.checkfile(args.conf_file):
            print("Config file {} not found!".format(args.conf_file))
            sys.exit()
        config = configparser.SafeConfigParser()
        config.optionxform = str
        config.read([args.conf_file])
        for section in config.sections():
            defaults.update(dict(config.items(section)))
    # Don't surpress add_help here so it will handle -h
    parser = argparse.ArgumentParser(
        # Inherit options from config_parser
        parents=[conf_parser],
        # print script description with -h/--help
        description=__doc__,
        # Don't mess with format of description
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

    parser.set_defaults(**defaults)
    # parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("images", type=str, nargs="+", help="Define the images to process")
    parser.add_argument("--output-tag", type=str, help="Add a tag to the output name.")
    parser.add_argument("--log-level", type=str2upper, choices=["WARNING", "INFO", "DEBUG"], help="Set the logging level.")
    parser.add_argument("--nice", type=str2int, help="Set the 'nice' level of processes.")
    parser.add_argument("--clobber", type=str2bool, help="Overwrite output if already exists.")
    parser.add_argument("--sumss-only", type=str2bool, help="Only use SUMSS in the image analysis.")
    parser.add_argument("--nvss-only", type=str2bool, help="Only use NVSS in the image analysis.")
    parser.add_argument("--weight-crop", type=str2bool, help="Crop image using the weights image.")
    parser.add_argument("--weight-crop-value", type=str2float, help="Define the minimum normalised value from the weights image to crop to.")
    parser.add_argument("--weight-crop-image", type=str, help="Define the weights image to use.")
    parser.add_argument("--convolve", type=str2bool, help="Convolve the image using CASA to SUMSS resolution for crossmatching.")
    parser.add_argument("--convolved-image", type=str, help="Define a convolved image that has already been produced.")
    parser.add_argument("--convolved-non-conv-askap-csv", type=str, help="Define the unconvolved catalogue to use when using convolved mode, \
        otherwise it will be generated automatically (if aegaen or pybdsf)")
    parser.add_argument("--convolved-non-conv-askap-islands-csv", type=str, help="Define the unconvolved island catalogue to use when using convolved mode, \
        otherwise it will be generated automatically (if aegaen or pybdsf)")
    parser.add_argument("--sourcefinder", type=str, choices=["aegean", "pybdsf", "selavy"], help="Select which sourcefinder to use")
    parser.add_argument("--frequency", type=str2float, help="Provide the frequency of the image in Hz. Use if 'RESTFRQ' is not in the header")
    # parser.add_argument("--fetch-sumss", type=str2bool, help="Fetch the SUMSS catalogue for the image area from Vizier.")
    parser.add_argument("--askap-csv", type=str, help="Manually define a aegean format csv file containing the extracted sources to use for the ASKAP image.")
    parser.add_argument("--askap-islands-csv", type=str, help="Manually define a csv file containing the extracted islands to use for the ASKAP image.")
    parser.add_argument("--sumss-csv", type=str, help="Manually provide the SUMSS catalog csv.")
    parser.add_argument("--nvss-csv", type=str, help="Manually provide the NVSS catalog csv.")
    parser.add_argument("--askap-csv-format", type=str, choices=["aegean", "selavy"], help="Define which source finder provided the ASKAP catalog (currently only supports aegean).")
    parser.add_argument("--remove-extended", type=str2bool,
        help="Remove perceived extended sources from the catalogues. Uses the following arguments 'askap-ext-thresh' and 'sumss-ext-thresh' to set the threshold.")
    parser.add_argument("--askap-ext-thresh", type=str2float,
        help="Define the maximum scaling threshold of the size of the ASKAP source compared to the PSF. Used to exclude extended sources. Only 1 axis has to exceed.")
    parser.add_argument("--sumss-ext-thresh", type=str2float,
        help="Define the maximum scaling threshold of the size of the SUMSS source compared to the PSF. Use to exclude extended sources. Only 1 axis has to exceed.")
    parser.add_argument("--nvss-ext-thresh", type=str2float,
        help="Define the maximum scaling threshold of the size of the NVSS source compared to the PSF. Use to exclude extended sources. Only 1 axis has to exceed.")
    #Option below provides framework for future source finder support, only aegaen for now.
    parser.add_argument("--use-all-fits", type=str2bool, help="Use all the fits from Aegean ignoring all flags. Default only those with flag '0' are used.")
    parser.add_argument("--write-ann", type=str2bool, help="Create kvis annotation files of the catalogues.")
    parser.add_argument("--produce-overlays", type=str2bool, help="Create overlay figures of the sources on the ASKAP image.")
    parser.add_argument("--boundary-value", type=str, choices=["nan", "zero"], help="Define whether the out-of-bounds value in the ASKAP FITS is 'nan' or 'zero'.")
    parser.add_argument("--askap-flux-error", type=str2float, help="Percentage error to apply to flux errors.")
    # parser.add_argument("--crossmatch-base", type=str, help="Define the base catalogue in the cross matching (currently not supported).", default="sumss")
    parser.add_argument("--diagnostic-max-separation", type=str2float, help="Maximum crossmatch distance (in arcsec) to be consdiered when creating the diagnostic plots.")
    parser.add_argument("--transient-max-separation", type=str2float, help="Maximum crossmatch distance (in arcsec) to be consdiered when searching for transients.")
    parser.add_argument("--postage-stamps", type=str2bool, help="Produce postage stamp plots of the cross matched sources within the max separation.")
    parser.add_argument("--postage-stamp-selection", type=str, choices=["all", "transients"], help="Select which postage stamps to create.")
    parser.add_argument("--postage-stamp-ncores", type=int, help="Select how many cores to use when creating the postage stamps.")
    parser.add_argument("--postage-stamp-radius", type=float, help="Select the radius of the postage stamp cutouts (arcmin).")
    parser.add_argument("--postage-stamp-zscale-contrast", type=float, help="Select the ZScale contrast to use in the postage stamps.")
    # parser.add_argument("--nprocs", type=str2int, help="Number of simulataneous SUMSS images at once when producing postage stamps.", default=1)
    parser.add_argument("--sumss-mosaic-dir", type=str, help="Directory containing the SUMSS survey mosaic image files.")
    parser.add_argument("--nvss-mosaic-dir", type=str, help="Directory containing the NVSS survey mosaic image files.")
    parser.add_argument("--aegean-settings-config", type=str, help="Select a config file containing the Aegean settings to be used (instead of defaults if none provided).")
    parser.add_argument("--pybdsf-settings-config", type=str, help="Select a config file containing the PyBDSF settings to be used (instead of defaults if none provided).")
    parser.add_argument("--selavy-settings-config", type=str, help="Select a config file containing the Selavy settings to be used (instead of defaults if none provided).")
    parser.add_argument("--transients", type=str2bool, help="Perform a transient search analysis using the crossmatch data. Requires '--max-separation' to be defined.")
    parser.add_argument("--transients-askap-snr-thresh", type=str2float, help="Define the threshold for which ASKAP sources are considered to not have a SUMSS match baseed\
     upon the estimated SUMSS SNR if the source was placed in the SUMSS image.")
    parser.add_argument("--transients-large-flux-ratio-thresh", type=str2float, help="Define the threshold for which sources are considered to have a large flux ratio.\
     Median value +/- threshold x std.")
    parser.add_argument("--db-inject", type=str2bool, help="Turn databse injection on or off.")
    parser.add_argument("--db-engine", type=str, help="Define the database engine.")
    parser.add_argument("--db-username", type=str, help="Define the username to use for the database")
    parser.add_argument("--db-password", type=str, help="Define the password to use for the database")
    parser.add_argument("--db-host", type=str, help="Define the host for the databse.")
    parser.add_argument("--db-port", type=str, help="Define the port for the databse.")
    parser.add_argument("--db-database", type=str, help="Define the name of the database.")
    parser.add_argument("--db-tag", type=str, help="The description field in the databased attached to the image.")
    parser.add_argument("--website-media-dir", type=str, help="Copy the image directory directly to the static media directory of the website.")
    # parser.add_argument("--dont-mask-sumss", type=str2bool, help="Do not filter the SUMSS catalogue such that only sources that should be in the ASKAP image remain.")
    args = parser.parse_args(remaining_argv)

    os.nice(args.nice)

    sf=args.sourcefinder

    logname="askapdiagnostic-{}".format(launchtime)
    utils.setup_logging(logname, args.log_level, use_colorlog)
    logger.debug("Debug on")
    #get user running
    username=getpass.getuser()

    if args.transients:
        if args.transient_max_separation==None:
            logger.error("Transients option selected but 'transient-max-separation' is not defined.")
            logger.error("Define the max-separation using '--transient-max-separation' and run again.")
            exit(logger)
        askap_snr_thresh=args.transients_askap_snr_thresh
        large_flux_ratio_thresh=args.transients_large_flux_ratio_thresh
        # if args.postage_stamps!=True:
        #     logger.warning("Transients option selected - turning on postage stamp production.")
        #     args.postage_stamps=True

    if args.weight_crop:
        if args.weight_crop_image == "None":
            logger.error("When using the weight crop feature, '--weight-crop-image' must be defined.")
            exit(logger)
        weight_image = os.path.abspath(args.weight_crop_image)
        if not utils.checkfile(weight_image):
            exit(logger)

    if args.convolve:
        if args.convolved_image=="None":
            if not utils.is_tool("casa"):
                logger.error("CASA not found in current PATH.")
                logger.error("To convolve an image CASA needs to be available. Please rectify and try again.")
                exit(logger)
        else:
            convolved_image = os.path.abspath(args.convolved_image)
            if not utils.checkfile(convolved_image):
                exit(logger)
        if args.convolved_non_conv_askap_csv != None:
            convolved_non_conv_askap_csv=os.path.abspath(args.convolved_non_conv_askap_csv)
            if not utils.checkfile(convolved_non_conv_askap_csv):
                exit(logger)
            else:
                non_convolved_src_cat=convolved_non_conv_askap_csv
                conv_source_find=False
        else:
            non_convolved_src_cat=None
            conv_source_find=True

        if args.convolved_non_conv_askap_islands_csv != None:
            convolved_non_conv_askap_islands_csv=os.path.abspath(args.convolved_non_conv_askap_islands_csv)
            if not utils.checkfile(convolved_non_conv_askap_islands_csv):
                exit(logger)
            else:
                non_convolved_isl_cat=convolved_non_conv_askap_islands_csv
        else:
            non_convolved_isl_cat=None

    else:
        conv_source_find = False

    for image in args.images:
        logger.info("Beginning processing of {}".format(image))
        #First check that image exists and get abs path
        image = os.path.abspath(image)
        if not utils.checkfile(image):
            exit(logger)

        #Load in the image and read in the information ready.
        theimg = Askapimage(image, readinfo=True)

        if theimg.freq == None:
            logger.info("Setting image frequency to {} Hz.".format(args.frequency))
            theimg.freq = args.frequency

        #Also get the rms
        base_image_rms=theimg.get_rms_clipping()

        #Create output name and check if output already exists (as now there is easy access to the image name).
        output="{}".format(theimg.imagename.replace(".fits", "_results"))
        if args.output_tag!="":
            output=output.replace("_results", "_results_{}".format(args.output_tag))

        #Get full path
        full_output=os.path.abspath(output)

        if not utils.createdir(output, clobber=args.clobber):
            exit(logger)

        #Fetch what catalogues we will be using
        sumss, nvss, dualmode, basecat, matched_to_tag = determine_catalogues(theimg.centre.dec.degree, args.sumss_only, args.nvss_only, logger)

        theimg.matched_to = matched_to_tag

        #And load the SUMSS beam info
        if sumss:
            theimg.calculate_sumss_beam()

        if args.postage_stamps:
            if sumss:
                if args.sumss_mosaic_dir == "None":
                    logger.error("SUMSS mosaic directory has not been defined for generating the postage stamps!")
                    logger.error("Define the SUMSS mosaic directory using '--sumss-mosaic-dir' and run again.")
                    exit(logger)
            if nvss:
                if args.nvss_mosaic_dir == None:
                    logger.error("NVSS mosaic directory has not been defined for generating the postage stamps!")
                    logger.error("Define the NVSS mosaic directory using '--nvss-mosaic-dir' and run again.")
                    exit(logger)

            logger.info("Will create the following postage stamps:")
            postage_options=args.postage_stamp_selection
            if postage_options=="all":
                logger.info("all")
                postage_options=["all"]
                if args.transients:
                    postage_options.append("transients")
            else:
                postage_options = [postage_options]

        #Check if catalogues have been provided
        if args.askap_csv != None:
            askap_cat_file = os.path.abspath(args.askap_csv)
            if not utils.checkfile(askap_cat_file):
                exit(logger)
            source_find=False
        else:
            logger.info("No ASKAP catalog provided - will perform source finding on the image using {}.".format(sf))
            source_find=True
            if sf == "selavy":
                if not utils.is_tool("selavy"):
                    logger.error("selavy cannot be found in the current path!. Please rectify and try again.")
                    exit(logger)

        if args.askap_islands_csv != None:
            askap_cat_islands_file = os.path.abspath(args.askap_islands_csv)
            if not utils.checkfile(askap_cat_islands_file):
                exit(logger)
        else:
            logger.warning("No ASKAP island catalog provided - will not perform island transient checks")
            askap_cat_islands_file = None

        if sumss:
            if args.sumss_csv != None:
                sumss_source_cat=os.path.abspath(args.sumss_csv)
                if not utils.checkfile(sumss_source_cat):
                    exit(logger)
                fetch_sumss=False
            else:
                logger.info("No SUMSS catalog provided - will perform fetch of SUMSS sources.")
                fetch_sumss=True
        else:
            sumss_source_cat=None

        if nvss:
            if args.nvss_csv != None:
                nvss_source_cat=os.path.abspath(args.nvss_csv)
                if not utils.checkfile(nvss_source_cat):
                    exit(logger)
                fetch_nvss=False
            else:
                logger.info("No NVSS catalog provided - will perform fetch of NVSS sources.")
                fetch_nvss=True
        else:
            nvss_source_cat=None

        #If aegean settings provided need to check that the file is present
        if source_find or conv_source_find:
            sf_option_files={"aegean":args.aegean_settings_config,
                        "pybdsf":args.pybdsf_settings_config,
                        "selavy":args.selavy_settings_config}
            sf_option_file = sf_option_files[sf]
            if sf_option_file != None:
                sf_option_file=os.path.abspath(sf_option_file)
                if not utils.checkfile(sf_option_file):
                    exit(logger)
                subprocess.call(["cp", sf_option_file, os.path.join(output, "{}_settings_used.cfg".format(sf))])
            else:
                sf_option_file = None

        #Now change directory
        os.chdir(output)

        if args.weight_crop:
            logger.info("Will crop image using weights image: {}. To a value of {}.".format(args.weight_crop_image, args.weight_crop_value))
            weight_cropped_image=theimg.weight_crop(weight_image, args.weight_crop_value)
            original_theimg = theimg
            theimg = Askapimage(weight_cropped_image, readinfo=True)
            theimg.original_name = original_theimg.imagename
            theimg.get_rms_clipping()

        if args.convolve:
            logger.info("Will convolve image to {} resolution".format(basecat.upper()))
            if not args.weight_crop:
                original_theimg = theimg
                non_convolved = theimg.image
            else:
                non_convolved = theimg.imagename
                original_theimg = None
            if non_convolved_src_cat != None:
                logger.info("Loading provided pre-convolved source catalogue: {}".format(non_convolved_src_cat))
                subprocess.call(["cp", non_convolved_src_cat, "."])
                sf_sigmas = [99.,99.]
                need_converting = True
            else:
                logger.info("Running sourcefinder on pre-convolved image...")
                need_converting = False
                non_convolved_src_cat=theimg.imagename.replace(".fits", "_comp.csv")
                sf_sigmas=source_finding(theimg, sf, logger, options=sf_option_file)
                if not os.path.isfile(non_convolved_src_cat):
                    logger.error("Source finding failed.")
                    exit(logger)

            if need_converting:
                non_conv_askap_cat = pd.read_fwf(non_convolved_src_cat, engine="python", skiprows=[1,])
            else:
                non_conv_askap_cat = pd.read_csv(non_convolved_src_cat, engine="python", delimiter=",")
            non_conv_askap_cat = Catalog(non_conv_askap_cat, "askap", "{}".format(theimg.imagename.replace(".fits", "_askap")), ra_col="ra",
                dec_col="dec", flux_col="int_flux", frequency=theimg.freq, add_name_col=True, convert_from_selavy=need_converting)
            non_conv_askap_cat._nearest_neighbour_d2d()
            non_conv_askap_cat._add_askap_sn()
            non_conv_askap_cat.add_error_uncertainty(args.askap_flux_error)
            if sumss:
                non_conv_askap_cat.calculate_scaled_flux("sumss", to_catalog="sumss")
                non_conv_askap_cat.calculate_scaled_flux("sumss_err", to_catalog="sumss", flux_col="err_int_flux")
                non_conv_askap_cat.add_sumss_sn(flux_col="askap_scaled_to_sumss")
            if nvss:
                non_conv_askap_cat.calculate_scaled_flux("nvss", to_catalog="nvss")
                non_conv_askap_cat.calculate_scaled_flux("nvss_err", to_catalog="nvss", flux_col="err_int_flux")
                non_conv_askap_cat.add_nvss_sn(flux_col="askap_scaled_to_nvss")

            if args.convolved_image=="None":
                if basecat=="sumss":
                    convolved_image = theimg.convolve2sumss()
                else:
                    convolved_image = theimg.convolve2sumss(nvss=True)
                if convolved_image == "Failed":
                    logger.error("Something went wrong with CASA convolving. Check the logs.")
                    exit(logger)
            else:
                logger.info("Using already supplied convolved image: {}.".format(convolved_image))

            loaded_catalog_data = {}

            loaded_catalog_data = non_conv_askap_cat.find_matching_catalog_image(
                sumss_mosaic_dir=args.sumss_mosaic_dir, nvss_mosaic_dir=args.nvss_mosaic_dir, sumss=sumss, nvss=nvss, dualmode=dualmode,
                loaded_data=loaded_catalog_data
            )
            non_conv_askap_cat.add_telescope_beam_columns("manual", manual_bmaj=theimg.bmaj*3600., manual_bmin=theimg.bmin*3600.)

            #Load up the islands if provided - only need to load this, no other processing needed.
            if non_convolved_isl_cat != None:
                non_convolved_isl_cat_df = pd.read_fwf(non_convolved_isl_cat, skiprows=[1,])
            else:
                non_convolved_isl_cat_df = pd.DataFrame([])

            theimg = Askapimage(convolved_image, readinfo=True)
            theimg.original_name = original_theimg.imagename
            theimg.non_convolved = non_convolved
            theimg.get_rms_clipping()
            theimg.calculate_sumss_beam()
            if theimg.freq == None:
                logger.warning("Frequency not found in convolved image.")
                logger.warning("Setting to {} Hz".format(original_theimg.freq))
                theimg.freq = original_theimg.freq

        else:
            non_conv_askap_cat = None
            non_convolved_isl_cat_df = None
            original_theimg = None

        #Perform source finding if needed
        if source_find:
            askap_cat_file=theimg.imagename.replace(".fits", "_comp.csv")
            sf_sigmas=source_finding(theimg, sf, logger, options=sf_option_file)
            #Check that the source finding was successful
            if not os.path.isfile(askap_cat_file):
                logger.error("Source finding failed.")
                exit(logger)
        else:
            logger.info("Using provided catalog for {}:".format(theimg.imagename))
            logger.info(askap_cat_file)
            subprocess.call(["cp", askap_cat_file, "."])
            sf_sigmas=[99.,99.]

        #Load the askap catalog into a new catalog object (but load it first to filter good sources)
        if source_find:
            askap_catalog=pd.read_csv(askap_cat_file, engine="python")
            need_converting = False
        else:
            askap_catalog=pd.read_fwf(askap_cat_file, engine="python", skiprows=[1,])
            need_converting = True
        aegean_total = len(askap_catalog.index)
        logger.info("Total number of ASKAP sources: {}".format(aegean_total))
        if not args.use_all_fits:
            askap_catalog=askap_catalog[askap_catalog["flag"]==0].reset_index(drop=True)
            aegean_good_fits=len(askap_catalog.index)
            logger.info("Total number of good fit ASKAP sources found: {} ({} sources removed).".format(aegean_good_fits, aegean_total - aegean_good_fits))
        theimg.total_askap_sources = len(askap_catalog.index)
        askap_catalog=Catalog(askap_catalog, "askap", "{}".format(theimg.imagename.replace(".fits", "_askap")), ra_col="ra", dec_col="dec",
            flux_col="int_flux", frequency=theimg.freq, add_name_col=True, convert_from_selavy=need_converting)
        askap_catalog._nearest_neighbour_d2d()
        askap_catalog._add_askap_sn()
        askap_catalog.add_error_uncertainty(args.askap_flux_error)
        askap_catalog.add_distance_from_pos(theimg.centre)
        askap_catalog.add_single_val_col("rms", theimg.rms)
        askap_catalog.add_telescope_beam_columns("manual", manual_bmaj=theimg.bmaj*3600., manual_bmin=theimg.bmin*3600.)
        #Below we also add the same columns to the convolved askap catalog, while the wrong way round this is to keep consistency with the postage stamps later.
        if args.convolve:
            askap_catalog.add_single_val_col("rms_preconv", original_theimg.rms)
            non_conv_askap_cat.add_single_val_col("rms", theimg.rms)
            non_conv_askap_cat.add_single_val_col("rms_preconv", original_theimg.rms)
            non_conv_askap_cat.add_distance_from_pos(theimg.centre)
            askap_catalog._merge_askap_non_convolved_catalogue(non_conv_askap_cat, sumss=sumss, nvss=nvss)

        #Load the askap convolved island file if available
        if askap_cat_islands_file != None:
            askap_cat_islands_df = pd.read_fwf(askap_cat_islands_file, skiprows=[1,])
        else:
            askap_cat_islands_df = pd.DataFrame([])

        #Find the respective catalog images for the askap sources
        #This will put results in 'catalog_Mosaic', 'catalog_Mosaic_path', 'catalog_Mosaic_rms'
        loaded_catalog_data = askap_catalog.find_matching_catalog_image(
            sumss_mosaic_dir=args.sumss_mosaic_dir, nvss_mosaic_dir=args.nvss_mosaic_dir, sumss=sumss, nvss=nvss, dualmode=dualmode, loaded_data=loaded_catalog_data
        )

        del loaded_catalog_data

        #Now check for SUMSS file or fetch if none
        if sumss:
            askap_catalog.calculate_scaled_flux("sumss", to_catalog="sumss")
            askap_catalog.calculate_scaled_flux("sumss_err", to_catalog="sumss", flux_col="err_int_flux")
            # askap_catalog.calculate_scaled_flux("sumss", flux_col="St", reverse=True)
            if fetch_sumss:
                sumss_source_cat, sumss_cat_df=get_catalogue(theimg, "SUMSS", args.boundary_value, logger)
                # sumss_source_cat=theimg.imagename.replace(".fits", "_sumss_comp.csv")
                # sumss_cat_df=theimg.get_catalogue("SUMSS", boundary_value=args.boundary_value)
                # if len(sumss_cat_df.index)<1:
                #     logger.error("SUMSS catalog fetching seems to have failed.")
                #     exit(logger)
                # #Write the filtered catalog to disk
                # theimg.write_sumss_sources()
                # logger.info("Fetching SUMSS catalogue complete.")
                # logger.info("Written to disk as {}.".format(sumss_source_cat))
            else:
                logger.info("Using provided SUMSS catalog for {}:".format(theimg.imagename))
                subprocess.call(["cp", sumss_source_cat, "."])
                logger.info(sumss_source_cat)
                sumss_cat_df=pd.read_csv(sumss_source_cat, delimiter=",", engine="python")

        if nvss:
            askap_catalog.calculate_scaled_flux("nvss", to_catalog="nvss")
            askap_catalog.calculate_scaled_flux("nvss_err", to_catalog="nvss", flux_col = "err_int_flux")
            # askap_catalog.calculate_scaled_flux("nvss", flux_col="S1.4", reverse=True)
            if fetch_nvss:
                nvss_source_cat, nvss_cat_df=get_catalogue(theimg, "NVSS", args.boundary_value, logger)
            else:
                logger.info("Using provided NVSS catalog for {}:".format(theimg.imagename))
                subprocess.call(["cp", nvss_source_cat, "."])
                logger.info(sumss_source_cat)
                nvss_cat_df=pd.read_csv(nvss_source_cat, delimiter=",", engine="python")

        #Initialise the SUMSS catalogue object
        if sumss:
            if fetch_sumss:
                raw_sumss=Catalog(theimg.raw_sumss_sources, "sumss", "{}".format(theimg.imagename.replace(".fits", "_raw_sumss")))
            sumss_catalog=Catalog(sumss_cat_df, "sumss", "{}".format(theimg.imagename.replace(".fits", "_sumss")), frequency=843.0e6, add_name_col=True)
            theimg.total_sumss_sources = len(sumss_cat_df.index)
            sumss_catalog.add_telescope_beam_columns("sumss")
            sumss_catalog._nearest_neighbour_d2d()
            sumss_catalog._add_sumss_mosaic_info(args.sumss_mosaic_dir)
            sumss_catalog.calculate_scaled_flux("askap", frequency=askap_catalog.frequency, flux_col="St")
            sumss_catalog.calculate_scaled_flux("askap_err", frequency=askap_catalog.frequency, flux_col="e_St")
            sumss_catalog.add_manual_sn(base_image_rms, "sumss_askap_snr", flux_col="sumss_scaled_to_askap", flux_scaling=1.e-3)
            sumss_catalog.add_sumss_sn(flux_col="St", dec_col="_DEJ2000", flux_scaling=1.e-3, use_image_rms=True)
            askap_catalog.add_sumss_sn(flux_col="askap_scaled_to_sumss", use_image_rms=True)
            #Filter the ASKAP catalogues if SUMSS only is used and we are near the 30 dec boundary
            if args.sumss_only:
                max_sumss = sumss_catalog.df["_DEJ2000"].max()
                if max_sumss >-31.0:
                    logger.warning("ASKAP image near the SUMSS border. Will check for ASKAP sources that are out of range when searching for transients.")
                    clean_for_sumss=True
                else:
                    clean_for_sumss=False
            else:
                clean_for_sumss=False
                max_sumss=0.0
        else:
            #Define empty ones for diagnostic plots
            sumss_catalog=[]
            theimg.total_sumss_sources = 0
            clean_for_sumss=False
            max_sumss = 0.0

        if nvss:
            if fetch_nvss:
                raw_nvss=Catalog(theimg.raw_nvss_sources, "nvss", "{}".format(theimg.imagename.replace(".fits", "_raw_nvss")))
            nvss_catalog=Catalog(nvss_cat_df, "nvss", "{}".format(theimg.imagename.replace(".fits", "_nvss")), flux_col="S1.4", frequency=1.4e9, add_name_col=True)
            theimg.total_nvss_sources = len(nvss_cat_df.index)
            nvss_catalog.add_telescope_beam_columns("nvss")
            nvss_catalog._nearest_neighbour_d2d()
            nvss_catalog._find_nvss_mosaic_info(args.nvss_mosaic_dir)
            nvss_catalog.calculate_scaled_flux("askap", frequency=askap_catalog.frequency, flux_col="S1.4")
            nvss_catalog.calculate_scaled_flux("askap_err", frequency=askap_catalog.frequency, flux_col="e_S1.4")
            nvss_catalog.add_manual_sn(base_image_rms, "nvss_askap_snr", flux_col="nvss_scaled_to_askap", flux_scaling=1.e-3)
            nvss_catalog.add_nvss_sn(flux_scaling=1.e-3, use_image_rms=True)
            askap_catalog.add_nvss_sn(flux_col="askap_scaled_to_nvss", use_image_rms=True)
            if args.nvss_only:
                min_nvss = nvss_catalog.df["_DEJ2000"].min()
                if min_nvss <-39.0:
                    logger.warning("ASKAP image near the NVSS border. Will check for ASKAP sources that are out of range when searching for transients.")
                    clean_for_nvss=True
                else:
                    clean_for_nvss=False
            else:
                clean_for_nvss=False
                min_nvss=0.0
        else:
            #Define empty ones for diagnostic plots
            nvss_catalog=[]
            theimg.total_nvss_sources = 0
            clean_for_nvss = False
            min_nvss = 0.0

        #Filter out extended sources if requested for diagnostics only
        if args.remove_extended:
            if sumss:
                noext_sumss_catalog=Catalog(sumss_catalog.remove_extended(threshold=args.sumss_ext_thresh,sumss_psf=True),
                     "sumss", "{}".format(theimg.imagename.replace(".fits", "_sumss_noext")), frequency=843.0e6)
            if nvss:
                noext_nvss_catalog=Catalog(nvss_catalog.remove_extended(threshold=args.nvss_ext_thresh,nvss_psf=True),
                     "nvss", "{}".format(theimg.imagename.replace(".fits", "_nvss_noext")), flux_col="S1.4", frequency=1.4e9)
                logger.info("NVSS No Ext Cat: {}".format(len(noext_nvss_catalog.df.index)))
            noext_askap_catalog=Catalog(askap_catalog.remove_extended(threshold=args.askap_ext_thresh,ellipse_a="a", ellipse_b="b",
                 beam_a=theimg.bmaj*3600., beam_b=theimg.bmin*3600.),
                 "askap", "{}".format(theimg.imagename.replace(".fits", "_askap_noext")),ra_col="ra", dec_col="dec", flux_col="int_flux", frequency=theimg.freq)

            if sumss:
                sumss_touse = noext_sumss_catalog
            if nvss:
                nvss_touse = noext_nvss_catalog
            askap_touse = noext_askap_catalog
        else:
            if sumss:
                sumss_touse = sumss_catalog
            if nvss:
                nvss_touse = nvss_catalog
            askap_touse = askap_catalog


        #Produce two plots of the image with overlay of ASKAP and SUMSS sources
        theimg.plots={}
        if args.produce_overlays:
            if args.remove_extended:
                if sumss:
                    theimg.plots["sumss_overlay"]=theimg.create_overlay_plot(sumss_touse.df, overlay_cat_label="SUMSS Sources", overlay_cat_2=sumss_catalog.sumss_ext_cat, overlay_cat_label_2="SUMSS Extended Sources", sumss=True)
                else:
                    theimg.plots["sumss_overlay"]="N/A"
                if nvss:
                    theimg.plots["nvss_overlay"]=theimg.create_overlay_plot(nvss_touse.df, overlay_cat_label="NVSS Sources", overlay_cat_2=nvss_catalog.nvss_ext_cat, overlay_cat_label_2="NVSS Extended Sources", sumss=False, nvss=True)
                else:
                    theimg.plots["nvss_overlay"]="N/A"
                theimg.plots["askap_overlay"]=theimg.create_overlay_plot(askap_touse.df, overlay_cat_label="ASKAP Extracted Sources", overlay_cat_2=askap_catalog.askap_ext_cat, overlay_cat_label_2="ASKAP Extended Extracted Sources")
            else:
                if sumss:
                    theimg.plots["sumss_overlay"]=theimg.create_overlay_plot(sumss_touse.df, overlay_cat_label="SUMSS Sources", sumss=True)
                else:
                    theimg.plots["sumss_overlay"]="N/A"
                if nvss:
                    theimg.plots["nvss_overlay"]=theimg.create_overlay_plot(nvss_touse.df, overlay_cat_label="NVSS Sources", sumss=False, nvss=True)
                else:
                    theimg.plots["nvss_overlay"]="N/A"
                theimg.plots["askap_overlay"]=theimg.create_overlay_plot(askap_touse.df, overlay_cat_label="ASKAP Extracted Sources")
        else:
            theimg.plots["sumss_overlay"]="N/A"
            theimg.plots["nvss_overlay"]="N/A"
            theimg.plots["askap_overlay"]="N/A"



        #Add the respective image to the ASKAP catalog for later
        askap_touse.add_single_val_col("image", theimg.image)
        askap_catalog.add_single_val_col("image", theimg.image)

        #Calculate distance from centre for each catalog for later
        # askap_touse.add_distance_from_pos(theimg.centre)

        if sumss:
            sumss_touse.add_distance_from_pos(theimg.centre)
            sumss_catalog.add_distance_from_pos(theimg.centre)
        if nvss:
            nvss_touse.add_distance_from_pos(theimg.centre)
            nvss_catalog.add_distance_from_pos(theimg.centre)

        #Annotation files
        if args.write_ann:
            if sumss:
                if fetch_sumss:
                    raw_sumss.write_ann(color="RED", name="sumss_raw.ann")
                sumss_catalog.write_ann(name="sumss.ann")
            if nvss:
                if fetch_nvss:
                    raw_nvss.write_ann(color="RED", name="nvss_raw.ann")
                nvss_catalog.write_ann(name="nvss.ann")
            askap_catalog.write_ann(color="BLUE", ellipse_a="a", ellipse_b="b", ellipse_pa="pa")

        #Start crossmatching section

        #First we do a crossmatch for diagnostics only
        logger.info("Performing diagnostic cross-match")
        # basecat=args.crossmatch_base
        # logger.info("Base catalogue: {}".format(basecat))

        if sumss:
            logger.debug("sumss_touse = {}".format(sumss_touse._crossmatch_catalog))
        if nvss:
            logger.debug("nvss_touse = {}".format(nvss_touse._crossmatch_catalog))
        logger.debug("askap_touse = {}".format(askap_touse._crossmatch_catalog))

        #Create new crossmatch object

        if dualmode or sumss:
            sumss_askap_crossmatch_diag=crossmatch(sumss_touse, askap_touse)
            sumss_askap_crossmatch_diag.perform_crossmatch()
            if dualmode:
                cross_name_tag="_sumss_nvss_askap_crossmatch_diagnostic.csv"
            else:
                cross_name_tag="_sumss_askap_crossmatch_diagnostic.csv"

        if dualmode or nvss:
            nvss_askap_crossmatch_diag=crossmatch(nvss_touse, askap_touse, base_catalogue_name="nvss")
            nvss_askap_crossmatch_diag.perform_crossmatch()
            if not dualmode:
                cross_name_tag="_nvss_askap_crossmatch_diagnostic.csv"

        #If dual mode the diagnostic match needs to be merged
        if dualmode:
            if basecat=="sumss":
                sumss_askap_crossmatch_diag.merge(nvss_askap_crossmatch_diag, basecat)
            else:
                nvss_askap_crossmatch_diag.merge(sumss_askap_crossmatch_diag, basecat)

        #Now we have it merged we can ignore the other one
        if dualmode:
            if basecat=="sumss":
                diag_crossmatch=sumss_askap_crossmatch_diag
            else:
                diag_crossmatch=nvss_askap_crossmatch_diag
        elif sumss:
            diag_crossmatch=sumss_askap_crossmatch_diag
        else:
            diag_crossmatch=nvss_askap_crossmatch_diag

        crossmatch_name=theimg.imagename.replace(".fits", cross_name_tag)

        #Calculate flux ratios and RA and Dec differences for plots
        logger.info("Calculating flux ratios and separations of crossmatches.")
        if dualmode:
            tag="sumss_nvss"
            # if basecat=="sumss":
                # askap_flux_col = "askap_int_flux_scale_sumss"
                # other_flux_col = "sumss_St"
                # other_ra_col = "sumss__RAJ2000"
                # other_dec_col = "sumss__DEJ2000"
            # else:
                # askap_flux_col = "askap_int_flux_scale_nvss"
                # other_flux_col = "nvss_S1.4"
                # other_ra_col = "nvss__RAJ2000"
                # other_dec_col = "nvss__DEJ2000"
        elif sumss:
            # askap_flux_col = "askap_int_flux_scale_sumss"
            # other_flux_col = "sumss_St"
            # other_ra_col = "sumss__RAJ2000"
            # other_dec_col = "sumss__DEJ2000"
            tag="sumss"
        else:
            # askap_flux_col = "askap_int_flux_scale_nvss"
            # other_flux_col = "nvss_S1.4"
            # other_ra_col = "nvss__RAJ2000"
            # other_dec_col = "nvss__DEJ2000"
            tag="nvss"

        diag_crossmatch.calculate_ratio("master_scaled_askap_iflux", "master_catflux", "askap_{}_int_flux_ratio".format(tag), col2_scaling=1.e-3, dualmode=dualmode, basecat=basecat)
        diag_crossmatch.calculate_ratio("master_catflux", "master_scaled_askap_iflux", "{}_askap_int_flux_ratio".format(tag), col1_scaling=1.e-3, dualmode=dualmode, basecat=basecat)
        diag_crossmatch.calculate_offsets(tag=tag)
        # diag_crossmatch.calculate_diff("askap_ra", "master_ra", "askap_{}_ra_offset".format(tag), dualmode=dualmode, basecat=basecat)

        # diag_crossmatch.calculate_diff("askap_dec", "master_dec", "askap_{}_dec_offset".format(tag), dualmode=dualmode, basecat=basecat)

        #Write out crossmatch df to file
        diag_crossmatch.write_crossmatch(crossmatch_name)

        #Get acceptable seperation df for plots
        if args.diagnostic_max_separation!=None:
            plotting_df = diag_crossmatch.crossmatch_df[diag_crossmatch.crossmatch_df["d2d"]<= args.diagnostic_max_separation]
        else:
            plotting_df = diag_crossmatch.crossmatch_df

        plotting_len=len(plotting_df.index)
        logger.debug("Length of plotting_df: {}".format(plotting_len))

        #Plots
        logger.info("Producing plots.")
        #Ratio view plot
        plot_titles={"source_counts":theimg.imagename+" source counts",
        }
        # if dualmode:
        #     plot_titles["flux_ratio_image_view"]=theimg.imagename+" ASKAP / (SUMSS & NVSS) flux ratio"
        #     plot_titles["position_offset"]=theimg.imagename+" SUMSS & NVSS Position Offset"
        #     plot_titles["flux_ratios_from_centre"]=theimg.imagename+" ASKAP / SUMSS & NVSS int. flux ratio vs. Distance from Image Centre"
        #     plot_titles["flux_ratios"]=theimg.imagename+" ASKAP / SUMSS & NVSS int. flux ratio vs. ASKAP Int Flux"
        # elif sumss:
        #     plot_titles["flux_ratio_image_view"]=theimg.imagename+" ASKAP / SUMSS flux ratio"
        #     plot_titles["position_offset"]=theimg.imagename+" SUMSS Position Offset"
        #     plot_titles["flux_ratios_from_centre"]=theimg.imagename+" ASKAP / SUMSS int. flux ratio vs. Distance from Image Centre"
        #     plot_titles["flux_ratios"]=theimg.imagename+" ASKAP / SUMSS int. flux ratio vs. ASKAP Int Flux"
        # else:
        #     plot_titles["flux_ratio_image_view"]=theimg.imagename+" ASKAP / NVSS flux ratio"
        #     plot_titles["position_offset"]=theimg.imagename+" NVSS Position Offset"
        #     plot_titles["flux_ratios_from_centre"]=theimg.imagename+" ASKAP / NVSS int. flux ratio vs. Distance from Image Centre"
        #     plot_titles["flux_ratios"]=theimg.imagename+" ASKAP / NVSS int. flux ratio vs. ASKAP Int Flux"

        plot_titles["flux_ratio_image_view"]=theimg.imagename
        plot_titles["position_offset"]=theimg.imagename
        plot_titles["flux_ratios_from_centre"]=theimg.imagename
        plot_titles["flux_ratios"]=theimg.imagename

        theimg.plots["flux_ratio_image_view"]=plots.flux_ratio_image_view_astropy(plotting_df, theimg.image, title=plot_titles["flux_ratio_image_view"], base_filename=theimg.imagename.replace(".fits", ""), basecat=basecat,
            ratio_col="askap_{}_int_flux_ratio".format(tag))
        theimg.plots["position_offset"]=plots.position_offset(plotting_df, title=plot_titles["position_offset"], base_filename=theimg.imagename.replace(".fits", ""),
            bmaj=theimg.bmaj*3600., bmin=theimg.bmin*3600., pa=theimg.bpa, basecat=basecat, ra_offset_col="askap_{}_ra_offset".format(tag),
            dec_offset_col="askap_{}_dec_offset".format(tag), dualmode=dualmode)
        # if basecat=="sumss":
        #     theimg.plots["source_counts"]=plots.source_counts(sumss_touse.df, askap_touse.df, diag_crossmatch.crossmatch_df, args.diagnostic_max_separation,
        #         title=plot_titles["source_counts"], base_filename=theimg.imagename.replace(".fits", ""))
        # else:
        theimg.plots["source_counts"]=plots.source_counts(askap_catalog.df, diag_crossmatch.crossmatch_df, args.diagnostic_max_separation,
            sumss_cat=sumss_catalog, nvss_cat=nvss_catalog, title=plot_titles["source_counts"], base_filename=theimg.imagename.replace(".fits", ""), basecat=basecat, dualmode=dualmode)
        theimg.plots["flux_ratios_from_centre"]=plots.flux_ratios_distance_from_centre(plotting_df, args.diagnostic_max_separation,
            title=plot_titles["flux_ratios_from_centre"], base_filename=theimg.imagename.replace(".fits", ""), basecat=basecat, ratio_col="askap_{}_int_flux_ratio".format(tag), dualmode=dualmode)
        theimg.plots["flux_ratios"]=plots.flux_ratios_askap_flux(plotting_df, args.diagnostic_max_separation,
            title=plot_titles["flux_ratios"], base_filename=theimg.imagename.replace(".fits", ""), basecat=basecat, ratio_col="askap_{}_int_flux_ratio".format(tag), dualmode=dualmode)

        #Now create crossmatch for the purpose of transient searching (and a crossmatch before convolving if required)
        logger.info("Performing transient cross-match")
        if dualmode or sumss:
            sumss_askap_crossmatch_transient=crossmatch(sumss_catalog, askap_catalog)
            sumss_askap_crossmatch_transient.perform_crossmatch()
        if dualmode or nvss:
            nvss_askap_crossmatch_transient=crossmatch(nvss_catalog, askap_catalog, base_catalogue_name="nvss")
            nvss_askap_crossmatch_transient.perform_crossmatch()


        if dualmode:
            if basecat=="sumss":
                sumss_askap_crossmatch_transient.merge(nvss_askap_crossmatch_transient, basecat)
            else:
                nvss_askap_crossmatch_transient.merge(sumss_askap_crossmatch_transient, basecat)

        if dualmode:
            if basecat=="sumss":
                transient_crossmatch=sumss_askap_crossmatch_transient
            else:
                transient_crossmatch=nvss_askap_crossmatch_transient
        elif sumss:
            transient_crossmatch=sumss_askap_crossmatch_transient
        else:
            transient_crossmatch=nvss_askap_crossmatch_transient

        if non_conv_askap_cat != None:
            if sumss:
                sumss_askap_preconv_crossmatch_transient=crossmatch(sumss_catalog, non_conv_askap_cat)
                sumss_askap_preconv_crossmatch_transient.perform_crossmatch()
            if nvss:
                nvss_askap_preconv_crossmatch_transient=crossmatch(nvss_catalog, non_conv_askap_cat, base_catalogue_name="nvss")
                nvss_askap_preconv_crossmatch_transient.perform_crossmatch()
            if dualmode:
                sumss_askap_preconv_crossmatch_transient.merge(nvss_askap_preconv_crossmatch_transient, basecat, raw=True)
                preconv_crossmatch_transient = sumss_askap_preconv_crossmatch_transient
            else:
                if sumss:
                    preconv_crossmatch_transient=sumss_askap_preconv_crossmatch_transient
                else:
                    preconv_crossmatch_transient=nvss_askap_preconv_crossmatch_transient

            #Also need to calculate all the ratios
            preconv_crossmatch_transient.calculate_ratio("master_scaled_askap_iflux", "master_catflux", "askap_cat_int_flux_ratio".format(tag), col2_scaling=1.e-3, dualmode=dualmode, basecat=basecat)
            preconv_crossmatch_transient.calculate_ratio("master_catflux", "master_scaled_askap_iflux", "cat_askap_int_flux_ratio".format(tag), col1_scaling=1.e-3, dualmode=dualmode, basecat=basecat)
            preconv_crossmatch_transient.calculate_offsets(tag=tag)
            # preconv_crossmatch_transient.calculate_diff("askap_ra", "master_ra", "askap_{}_ra_offset".format(tag), dualmode=dualmode, basecat=basecat)
            # preconv_crossmatch_transient.calculate_diff("askap_dec", "master_dec", "askap_{}_dec_offset".format(tag), dualmode=dualmode, basecat=basecat)
        else:
            preconv_crossmatch_transient = None
            preconv_crossmatch_transient = None

        crossmatch_name=theimg.imagename.replace(".fits", "_{}_askap_crossmatch_transient.csv".format(tag))
        logger.info("Calculating flux ratios and separations of crossmatches.")
        transient_crossmatch.calculate_ratio("master_scaled_askap_iflux", "master_catflux", "askap_cat_int_flux_ratio".format(tag), col2_scaling=1.e-3, dualmode=dualmode, basecat=basecat)
        transient_crossmatch.calculate_ratio("master_catflux", "master_scaled_askap_iflux", "cat_askap_int_flux_ratio".format(tag), col1_scaling=1.e-3, dualmode=dualmode, basecat=basecat)
        transient_crossmatch.calculate_offsets(tag=tag)
        # transient_crossmatch.calculate_diff("askap_ra", "master_ra", "askap_{}_ra_offset".format(tag), dualmode=dualmode, basecat=basecat)
        # transient_crossmatch.calculate_diff("askap_dec", "master_dec", "askap_{}_dec_offset".format(tag), dualmode=dualmode, basecat=basecat)

        transient_crossmatch.get_good_matches(maxsep=args.transient_max_separation)
        transient_crossmatch.get_bad_matches(maxsep=args.transient_max_separation)


        if args.transients:
            #Need to check if convolving and pass empty values if convolving is not used.
            if args.convolve:
                preconv_askap_img_wcs=original_theimg.wcs
                preconv_askap_img_data=original_theimg.data
                preconv_askap_img_header=original_theimg.header
            else:
                preconv_askap_img_wcs="None"
                preconv_askap_img_data=[]
                preconv_askap_img_header={}
            transient_crossmatch.transient_search(max_separation=args.transient_max_separation, askap_snr_thresh=askap_snr_thresh, large_flux_thresh=large_flux_ratio_thresh,
                pre_conv_crossmatch=preconv_crossmatch_transient, image_beam_maj=theimg.bmaj*3600., image_beam_min=theimg.bmin*3600., image_beam_pa=theimg.bpa, dualmode=dualmode, sumss=sumss, nvss=nvss,
                askap_img_wcs=theimg.wcs, askap_img_header=theimg.header, askap_img_data=theimg.data, preconv_askap_img_wcs=preconv_askap_img_wcs, preconv_askap_img_header=preconv_askap_img_header,
                preconv_askap_img_data=preconv_askap_img_data,askap_cat_islands_df=askap_cat_islands_df, non_convolved_isl_cat_df=non_convolved_isl_cat_df,
                askap_image=theimg, preconv_askap_image=original_theimg, clean_for_sumss=clean_for_sumss, max_sumss=max_sumss, clean_for_nvss=clean_for_nvss, min_nvss=min_nvss)
            os.makedirs("transients/no-match")
            os.makedirs("transients/large-ratio")
            os.makedirs("transients/askap-notseen")
            subprocess.call("mv transients*.csv transients/", shell=True)

        #Write out crossmatch df to file
        transient_crossmatch.write_crossmatch(crossmatch_name)

        if args.postage_stamps:
            logger.info("Starting postage stamp production.")
            transient_crossmatch.produce_postage_stamps(theimg.data, theimg.wcs, "all", args.postage_stamp_ncores, radius=args.postage_stamp_radius, contrast=args.postage_stamp_zscale_contrast,
            convolve=args.convolve, askap_nonconv_image=original_theimg, askap_pre_convolve_catalog=non_conv_askap_cat, dualmode=dualmode,
            basecat=basecat)
            os.mkdir("postage-stamps")
 #            os.makedirs("postage-stamps/bad")
 #            if args.transients:
 #                try:
 #                    subprocess.check_output("mv transient*NOMATCH*_sidebyside.jpg transients/no-match/", shell=True, stderr=subprocess.STDOUT)
 #                except subprocess.CalledProcessError as e:
 #                    logger.warning("No transient 'NO MATCH' images to move.")
 #                try:
 #                    subprocess.check_output("mv transient*LARGERATIO*_sidebyside.jpg transients/large-ratio/", shell=True, stderr=subprocess.STDOUT)
 #                except subprocess.CalledProcessError as e:
 #                    logger.warning("No transient 'LARGE RATIO' images to move.")
 #                try:
 #                    subprocess.check_output("mv transient_askapnotseen*_sidebyside.jpg transients/askap-notseen/", shell=True, stderr=subprocess.STDOUT)
 #                except subprocess.CalledProcessError as e:
 #                    logger.warning("No transient 'ASKAP NOT SEEN' images to move.")

            subprocess.call("mv *_postagestamps.jpg postage-stamps/", shell=True)
            # subprocess.call("mv *BAD*_sidebyside.jpg postage-stamps/bad/", shell=True)

        #Hack fix for now CHANGE!
        # theimg.rms=base_image_rms
        # Database Entry
        theimg.matched_to = matched_to_tag
        # image_id=theimg.inject_db(datestamp=launchtime)
        #Processing settings
        if args.transients:
            transients_noaskapmatchtocatalog_total=transient_crossmatch.transients_noaskapmatchtocatalog_total
            transients_noaskapmatchtocatalog_candidates=transient_crossmatch.transients_noaskapmatchtocatalog_candidates
            transients_nocatalogmatchtoaskap_total=transient_crossmatch.transients_nocatalogmatchtoaskap_total
            transients_nocatalogmatchtoaskap_candidates=transient_crossmatch.transients_nocatalogmatchtoaskap_candidates
            # transients_largeratio_total=transient_crossmatch.transients_largeratio_total
            # transients_largeratio_candidates=transient_crossmatch.transients_largeratio_candidates
            transients_goodmatches_total=transient_crossmatch.transients_goodmatches_total
            transients_master_total=transient_crossmatch.transients_master_total
            transients_master_candidates_total = transient_crossmatch.transients_master_candidates_total
            transients_master_flagged_total = transient_crossmatch.transients_master_flagged_total
        else:
            transients_noaskapmatchtocatalog_total=0
            transients_noaskapmatchtocatalog_candidates=0
            transients_nocatalogmatchtoaskap_total=0
            transients_nocatalogmatchtoaskap_candidates=0
            # transients_largeratio_total=0
            # transients_largeratio_candidates=0
            transients_goodmatches_total=0
            transients_master_total=0
            transients_master_candidates_total = 0
            transients_master_flagged_total = 0
        if args.db_inject:
            if args.convolve:
                image_2 = original_theimg.image
            else:
                image_2 = "N/A"
            image_id=theimg.inject_db(basecat=basecat, datestamp=launchtime, user=username, description=args.db_tag, db_engine=args.db_engine, db_username=args.db_username,
                db_host=args.db_host, db_password=args.db_password,
                db_port=args.db_port, db_database=args.db_database, transients_noaskapmatchtocatalog_total=transients_noaskapmatchtocatalog_total,
                transients_noaskapmatchtocatalog_candidates=transients_noaskapmatchtocatalog_candidates,
                transients_nocatalogmatchtoaskap_total=transients_nocatalogmatchtoaskap_total,
                transients_nocatalogmatchtoaskap_candidates=transients_nocatalogmatchtoaskap_candidates,
                # transients_largeratio_total=transients_largeratio_total,
                # transients_largeratio_candidates=transients_largeratio_candidates,
                transients_goodmatches_total=transients_goodmatches_total,
                transients_master_total=transients_master_total,transients_master_candidates_total=transients_master_candidates_total, transients_master_flagged_total=transients_master_flagged_total, image_2 = image_2)
            theimg.inject_processing_db(image_id, full_output, askap_cat_file, sumss_source_cat, nvss_source_cat, args.askap_ext_thresh,
                args.sumss_ext_thresh, args.nvss_ext_thresh, args.transient_max_separation, sf_sigmas, db_engine=args.db_engine, db_username=args.db_username, db_password=args.db_password, db_host=args.db_host,
                db_port=args.db_port, db_database=args.db_database)
            if args.transients:
                transient_crossmatch.inject_transients_db(image_id, sumss, nvss, db_engine=args.db_engine, db_username=args.db_username, db_password=args.db_password, db_host=args.db_host,
                    db_port=args.db_port, db_database=args.db_database, max_separation=args.transient_max_separation, dualmode=dualmode, basecat=basecat, askap_image = theimg.image,
                    askap_image_2 = image_2)
        else:
            image_id = None
        # if args.transients:
            # transient_crossmatch.inject_transients_db(image_id, db_engine=args.db_engine, db_username=args.db_username, db_host=args.db_host,
            # db_port=args.db_port, db_database=args.db_database, dualmode=dualmode, basecat=basecat)

        if args.website_media_dir!="none" and image_id is not None:
            media_dir=os.path.join(args.website_media_dir, str(image_id))
        # else:
            # media_dir=os.path.join("..", "static", "media", "{}".format(image_id))
            stamp_media_dir=os.path.join(media_dir, "stamps")
        # os.makedirs(media_dir)
            os.makedirs(stamp_media_dir)
            subprocess.call("for i in `ls postage-stamps/*.jpg`; do cp $i {}/ ; done".format(stamp_media_dir), shell=True)
        # subprocess.call("cp postage-stamps/bad/*.jpg {}/".format(stamp_media_dir), shell=True)
        # if args.transients:
        #     subprocess.call("cp transients/askap-notseen/*.jpg {}/".format(stamp_media_dir), shell=True)
            subprocess.call("cp *.png {}/".format(media_dir), shell=True)

        os.chdir("..")

        logger.info("Analysis complete for {}.".format(theimg.imagename))


    logger.info("All images done. Log will be placed in {}.".format(output))
    subprocess.call(["cp", logname+".log", output])
    subprocess.call(["rm", logname+".log"])
