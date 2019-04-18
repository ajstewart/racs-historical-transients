#!/usr/bin/env python

import argparse
import logging
from askapdiagnostic.tools import utils
from askapdiagnostic.tools.fitsimage import askapimage
from askapdiagnostic.tools.crossmatch import crossmatch
from askapdiagnostic.tools.catalog import Catalog
from askapdiagnostic.plotting import plots
import datetime
import os
import pandas as pd
import sys
import subprocess
import ConfigParser
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
        config = ConfigParser.SafeConfigParser()
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
        "weight_crop":False,
        "weight_crop_value":0.04,
        "weight_crop_image":"weights.fits",
        "convolve":False,
        "convolved_image":"convolved.fits",
        "convolved_non_conv_askap_csv":None,
        "sourcefinder":"aegean",
        "frequency":864e6,
        "askap_csv":None,
        "sumss_csv":None,
        "nvss_csv":None,
        "askap_csv_format":"aegean",
        "remove_extended":False,
        "askap_ext_thresh":1.2,
        "sumss_ext_thresh":1.2,
        "nvss_ext_thresh":1.2,
        "use_all_fits":False,
        "write_ann":False,
        "boundary_value":"nan",
        "diagnostic_max_separation":5.0,
        "transient_max_separation":45.0,
        "postage_stamps":False,
        "postage_stamp_selection":["all",],
        "sumss_mosaic_dir":None,
        "nvss_mosaic_dir":None,
        "aegean_settings_config":None,
        "pybdsf_settings_config":None,
        "selavy_settings_config":None,
        "transients":False,
        "transients_askap_sumss_snr_thresh":5.0,
        "transients_large_flux_ratio_thresh":3.0,
        "db_engine":"postgresql",
        "db_username":"postgres",
        "db_host":"localhost",
        "db_port":"5432",
        "db_database":"postgres",
        "db_tag":"RACS Analysis",
        "website_media_dir":"none"
        }
        
    if args.conf_file:
        config = ConfigParser.SafeConfigParser()
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
    parser.add_argument("--log-level", type=str, choices=["WARNING", "INFO", "DEBUG"], help="Set the logging level.")
    parser.add_argument("--nice", type=str2int, help="Set the 'nice' level of processes.")
    parser.add_argument("--clobber", type=str2bool, help="Overwrite output if already exists.")
    parser.add_argument("--weight-crop", type=str2bool, help="Crop image using the weights image.")
    parser.add_argument("--weight-crop-value", type=str2float, help="Define the minimum normalised value from the weights image to crop to.")
    parser.add_argument("--weight-crop-image", type=str, help="Define the weights image to use.")
    parser.add_argument("--convolve", type=str2bool, help="Convolve the image using CASA to SUMSS resolution for crossmatching.")
    parser.add_argument("--convolved-image", type=str, help="Define a convolved image that has already been produced.")
    parser.add_argument("--convolved-non-conv-askap-csv", type=str, help="Define the unconvolved catalogue to use when using convolved mode, \
        otherwise it will be generated automatically (if aegaen or pybdsf)")
    parser.add_argument("--sourcefinder", type=str, choices=["aegean", "pybdsf", "selavy"], help="Select which sourcefinder to use")
    parser.add_argument("--frequency", type=str2float, help="Provide the frequency of the image in Hz. Use if 'RESTFRQ' is not in the header")
    # parser.add_argument("--fetch-sumss", type=str2bool, help="Fetch the SUMSS catalogue for the image area from Vizier.")
    parser.add_argument("--askap-csv", type=str, help="Manually define a aegean csv file containing the extracted sources to use for the ASKAP image.")
    parser.add_argument("--sumss-csv", type=str, help="Manually provide the SUMSS catalog csv.")
    parser.add_argument("--nvss-csv", type=str, help="Manually provide the NVSS catalog csv.")
    parser.add_argument("--askap-csv-format", type=str, choices=["aegean"], help="Define which source finder provided the ASKAP catalog (currently only supports aegean).")
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
    parser.add_argument("--boundary-value", type=str, choices=["nan", "zero"], help="Define whether the out-of-bounds value in the ASKAP FITS is 'nan' or 'zero'.")
    # parser.add_argument("--crossmatch-base", type=str, help="Define the base catalogue in the cross matching (currently not supported).", default="sumss")
    parser.add_argument("--diagnostic-max-separation", type=str2float, help="Maximum crossmatch distance (in arcsec) to be consdiered when creating the diagnostic plots.")
    parser.add_argument("--transient-max-separation", type=str2float, help="Maximum crossmatch distance (in arcsec) to be consdiered when searching for transients.")
    parser.add_argument("--postage-stamps", type=str2bool, help="Produce postage stamp plots of the cross matched sources within the max separation.")
    parser.add_argument("--postage-stamp-selection", type=convert_str_to_list, choices=["all", "good", "bad", "transients"], help="Select which postage stamps to create.")
    # parser.add_argument("--nprocs", type=str2int, help="Number of simulataneous SUMSS images at once when producing postage stamps.", default=1)
    parser.add_argument("--sumss-mosaic-dir", type=str, help="Directory containing the SUMSS survey mosaic image files.")
    parser.add_argument("--nvss-mosaic-dir", type=str, help="Directory containing the NVSS survey mosaic image files.")
    parser.add_argument("--aegean-settings-config", type=str, help="Select a config file containing the Aegean settings to be used (instead of defaults if none provided).")
    parser.add_argument("--pybdsf-settings-config", type=str, help="Select a config file containing the PyBDSF settings to be used (instead of defaults if none provided).")
    parser.add_argument("--selavy-settings-config", type=str, help="Select a config file containing the Selavy settings to be used (instead of defaults if none provided).")
    parser.add_argument("--transients", type=str2bool, help="Perform a transient search analysis using the crossmatch data. Requires '--max-separation' to be defined.")
    parser.add_argument("--transients-askap-sumss-snr-thresh", type=str2float, help="Define the threshold for which ASKAP sources are considered to not have a SUMSS match baseed\
     upon the estimated SUMSS SNR if the source was placed in the SUMSS image.")
    parser.add_argument("--transients-large-flux-ratio-thresh", type=str2float, help="Define the threshold for which sources are considered to have a large flux ratio.\
     Median value +/- threshold x std.")
    parser.add_argument("--db-engine", type=str, help="Define the database engine.")
    parser.add_argument("--db-username", type=str, help="Define the username to use for the database")
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
    
    #get user running
    username=getpass.getuser()
    
    if args.postage_stamps:
        if args.sumss_mosaic_dir == None:
            logger.error("SUMSS mosaic directory has not been defined for generating the postage stamps!")
            logger.error("Define the SUMSS mosaic directory using '--sumss-mosaic-dir' and run again.")
            exit(logger)
        logger.info("Will create the following postage stamps:")
        postage_options=args.postage_stamp_selection
        if "all" in postage_options:
            logger.info("all")
            postage_options=["all"]
            if args.transients:
                postage_options.append("transients")
        else:
            if "good" in postage_options and "bad" in postage_options:
                if "transients" in postage_options:
                    postage_options=["all", "transients"]
                else:
                    postage_options=["all"]
            if "transients" in postage_options and not args.transients:
                logger.error("'transients' included in postage stamp options but '--transients' is not selected.")
                logger.error("Please confirm settings and launch again")
                exit(logger)
            else:
                for i in postage_options:
                    logger.info(i)
        
    if args.transients:
        if args.transient_max_separation==None:
            logger.error("Transients option selected but 'transient-max-separation' is not defined.")
            logger.error("Define the max-separation using '--transient-max-separation' and run again.")
            exit(logger)
        askap_sumss_snr_thresh=args.transients_askap_sumss_snr_thresh
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
    else:
        conv_source_find = False
    
    for image in args.images:
        logger.info("Beginning processing of {}".format(image))
        #First check that image exists and get abs path
        image = os.path.abspath(image)
        if not utils.checkfile(image):
            exit(logger)
        
        #Load in the image and read in the information ready.
        theimg = askapimage(image, readinfo=True)
        
        if theimg.freq == None:
            logger.info("Setting image frequency to {} Hz.".format(args.frequency))
            theimg.freq = args.frequency
        
        #Also get the rms
        theimg.get_rms_clipping()
        
        
        #Create output name and check if output already exists (as now there is easy access to the image name).
        output="{}".format(theimg.imagename.replace(".fits", "_results"))
        if args.output_tag!="":
            output=output.replace("_results", "_results_{}".format(args.output_tag))
        
        #Get full path
        full_output=os.path.abspath(output)
            
        if not utils.createdir(output, clobber=args.clobber):
            exit(logger)
        
        #Check if NVSS or SUMSS is required (or both)
        if theimg.centre.dec.degree <= -30.0:
            basecat="sumss"
            sumss=True
            if theimg.centre.dec.degree < -43:
                nvss=False
            else:
                nvss=True
        else:
            basecat="nvss"
            nvss=True
            sumss=False

        
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
            
        if sumss:
            if args.sumss_csv != None:
                sumss_source_cat=os.path.abspath(args.sumss_csv)
                if not utils.checkfile(sumss_source_cat):
                    exit(logger)
                fetch_sumss=False
            else:
                logger.info("No SUMSS catalog provided - will perform fetch of SUMSS sources.")
                fetch_sumss=True
            
        if nvss:
            if args.nvss_csv != None:
                nvss_source_cat=os.path.abspath(args.nvss_csv)
                if not utils.checkfile(nvss_source_cat):
                    exit(logger)
                fetch_nvss=False
            else:
                logger.info("No NVSS catalog provided - will perform fetch of NVSS sources.")
                fetch_nvss=True
                
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
            theimg = askapimage(weight_cropped_image, readinfo=True)
            theimg.original_name = original_theimg.imagename
            theimg.get_rms_clipping()
            
        if args.convolve:
            logger.info("Will covolve image to SUMSS resolution")
            if not args.weight_crop:
                original_theimg = theimg
            non_convolved = theimg.imagename
            if non_convolved_src_cat != None:
                logger.info("Loading provided pre-convolved source catalogue: {}".format(non_convolved_src_cat))
                subprocess.call(["cp", non_convolved_src_cat, "."])
                sf_sigmas=[99.,99.]
            else:
                logger.info("Running sourcefinder on pre-convolved image...")
                non_convolved_src_cat=theimg.imagename.replace(".fits", "_comp.csv")
                sf_sigmas=source_finding(theimg, sf, logger, options=sf_option_file)
                if not os.path.isfile(askap_no_conv_cat_file):
                    logger.error("Source finding failed.")
                    exit(logger)

            non_conv_askap_cat = pd.read_csv(non_convolved_src_cat, delimiter=",", engine="python")
            non_conv_askap_cat = Catalog(non_conv_askap_cat, "{}".format(theimg.imagename.replace(".fits", "_askap")), ra_col="ra", dec_col="dec", add_name_col=True)
                
            if args.convolved_image=="None":
                convolved_image = theimg.convolve2sumss()
                if convolved_image == "Failed":
                    logger.error("Something went wrong with CASA convolving. Check the logs.")
                    exit(logger)
            else:
                logger.info("Using already supplied convolved image: {}.".format(convolved_image))
                
            theimg = askapimage(convolved_image, readinfo=True)
            theimg.original_name = original_theimg.imagename
            theimg.non_convolved = non_convolved
            theimg.get_rms_clipping()
        
        else:
            non_conv_askap_cat = None
        
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
        askap_catalog=pd.read_csv(askap_cat_file, delimiter=",", engine="python")
        aegean_total = len(askap_catalog.index)
        logger.info("Total number of sources found by Aegean: {}".format(aegean_total))
        if not args.use_all_fits:
            askap_catalog=askap_catalog[askap_catalog["flags"]==0].reset_index(drop=True)
            aegean_good_fits=len(askap_catalog.index)
            logger.info("Total number of good fit sources found by Aegean: {} ({} sources removed).".format(aegean_good_fits, aegean_total - aegean_good_fits))
        theimg.total_askap_sources = len(askap_catalog.index)
        askap_catalog=Catalog(askap_catalog, "{}".format(theimg.imagename.replace(".fits", "_askap")), ra_col="ra", dec_col="dec", add_name_col=True)
            
        #Now check for SUMSS file or fetch if none
        if fetch_sumss:
            sumss_source_cat=theimg.imagename.replace(".fits", "_sumss_comp.csv")
            sumss_cat_df=theimg.get_sumss_catalogue(boundary_value=args.boundary_value)
            if len(sumss_cat_df.index)<1:
                logger.error("SUMSS catalog fetching seems to have failed.")
                exit(logger)
            #Write the filtered catalog to disk
            theimg.write_sumss_sources()
            logger.info("Fetching SUMSS catalogue complete.")
            logger.info("Written to disk as {}.".format(sumss_source_cat))
        else:
            logger.info("Using provided SUMSS catalog for {}:".format(theimg.imagename))
            subprocess.call(["cp", sumss_source_cat, "."])
            logger.info(sumss_source_cat)
            sumss_cat_df=pd.read_csv(sumss_source_cat, delimiter=",", engine="python")
        
        #Initialise the SUMSS catalogue object
        if fetch_sumss:
            raw_sumss=Catalog(theimg.raw_sumss_sources, "{}".format(theimg.imagename.replace(".fits", "_raw_sumss")))
        sumss_catalog=Catalog(sumss_cat_df, "{}".format(theimg.imagename.replace(".fits", "_sumss")), add_name_col=True)
        theimg.total_sumss_sources = len(sumss_cat_df.index)
        
        #Filter out extended sources if requested for diagnostics only
        if args.remove_extended:
            noext_sumss_catalog=Catalog(sumss_catalog.remove_extended(threshold=args.sumss_ext_thresh,sumss_psf=True),
                 "{}".format(theimg.imagename.replace(".fits", "_sumss_noext")))
            noext_askap_catalog=Catalog(askap_catalog.remove_extended(threshold=args.askap_ext_thresh,ellipse_a="a", ellipse_b="b", 
                 beam_a=theimg.bmaj*3600., beam_b=theimg.bmin*3600.),
                 "{}".format(theimg.imagename.replace(".fits", "_askap_noext")),ra_col="ra", dec_col="dec")

            sumss_touse = noext_sumss_catalog
            askap_touse = noext_askap_catalog
        else:
            sumss_touse = sumss_catalog
            askap_touse = askap_catalog
        
        
        #Produce two plots of the image with overlay of ASKAP and SUMSS sources
        theimg.plots={}
        if args.remove_extended:
            theimg.plots["sumss_overlay"]=theimg.create_overlay_plot(sumss_touse.df, overlay_cat_label="SUMSS Sources", overlay_cat_2=sumss_catalog.sumss_ext_cat, overlay_cat_label_2="SUMSS Extended Sources", sumss=True)
            theimg.plots["askap_overlay"]=theimg.create_overlay_plot(askap_touse.df, overlay_cat_label="ASKAP Extracted Sources", overlay_cat_2=askap_catalog.askap_ext_cat, overlay_cat_label_2="ASKAP Extended Extracted Sources")
        else:
            theimg.plots["sumss_overlay"]=theimg.create_overlay_plot(sumss_touse.df, overlay_cat_label="SUMSS Sources", sumss=True)
            theimg.plots["askap_overlay"]=theimg.create_overlay_plot(askap_touse.df, overlay_cat_label="ASKAP Extracted Sources")
            
        
        
        
        #Add the respective image to the ASKAP catalog for later
        askap_touse.add_single_val_col("image", theimg.image)
        askap_catalog.add_single_val_col("image", theimg.image)
        
        #Calculate distance from centre for each catalog for later
        askap_touse.add_distance_from_pos(theimg.centre)
        askap_catalog.add_distance_from_pos(theimg.centre)
        sumss_touse.add_distance_from_pos(theimg.centre)
        sumss_catalog.add_distance_from_pos(theimg.centre)
        
        #Add SUMSS S/N for ASKAP sources
        askap_touse.add_sumss_sn(flux_col="int_flux")
        askap_catalog.add_sumss_sn(flux_col="int_flux")
        sumss_touse.add_sumss_sn(flux_col="St", dec_col="_DEJ2000", flux_scaling=1.e-3)
        sumss_catalog.add_sumss_sn(flux_col="St", dec_col="_DEJ2000", flux_scaling=1.e-3)
        askap_touse._add_askap_sn()
        askap_catalog._add_askap_sn()
        
        #Annotation files
        if args.write_ann:
            if fetch_sumss:
                raw_sumss.write_ann(color="RED")
            sumss_catalog.write_ann()
            askap_catalog.write_ann(color="BLUE", ellipse_a="a", ellipse_b="b", ellipse_pa="pa")
        
        #Start crossmatching section
        
        #First we do a crossmatch for diagnostics only
        logger.info("Performing diagnostic cross-match")
        # basecat=args.crossmatch_base
        # logger.info("Base catalogue: {}".format(basecat))
        
        logger.debug("sumss_touse = {}".format(sumss_touse._crossmatch_catalog))
        logger.debug("askap_touse = {}".format(askap_touse._crossmatch_catalog))
        
        #Create new crossmatch object
        sumss_askap_crossmatch_diag=crossmatch(sumss_touse, askap_touse)
        sumss_askap_crossmatch_diag.perform_crossmatch(maxsep=args.diagnostic_max_separation)
        
        crossmatch_name=theimg.imagename.replace(".fits", "_sumss_askap_crossmatch_diagnostic.csv")
        
        #Calculate flux ratios and RA and Dec differences for plots
        logger.info("Calculating flux ratios and separations of crossmatches.")
        sumss_askap_crossmatch_diag.calculate_ratio("askap_int_flux", "sumss_St", "askap_sumss_int_flux_ratio", col2_scaling=1.e-3)
        sumss_askap_crossmatch_diag.calculate_ratio("sumss_St", "askap_int_flux", "sumss_askap_int_flux_ratio", col1_scaling=1.e-3)
        sumss_askap_crossmatch_diag.calculate_diff("askap_ra", "sumss__RAJ2000", "askap_sumss_ra_offset")
        sumss_askap_crossmatch_diag.calculate_diff("askap_dec", "sumss__DEJ2000", "askap_sumss_dec_offset")
        
        #Write out crossmatch df to file
        sumss_askap_crossmatch_diag.write_crossmatch(crossmatch_name)
        
        #Get acceptable seperation df for plots
        if args.diagnostic_max_separation!=None:
            plotting_df = sumss_askap_crossmatch_diag.crossmatch_df[sumss_askap_crossmatch_diag.crossmatch_df["d2d"]<= args.diagnostic_max_separation]
        else:
            plotting_df = sumss_askap_crossmatch_diag.crossmatch_df
            
        plotting_len=len(plotting_df.index)
        logger.debug("Length of plotting_df: {}".format(plotting_len))
        
        #Plots    
        logger.info("Producing plots.")
        #Ratio view plot
        theimg.plots["flux_ratio_image_view"]=plots.flux_ratio_image_view(plotting_df, title=theimg.imagename+" ASKAP / SUMSS flux ratio", base_filename=theimg.imagename.replace(".fits", ""))
        theimg.plots["position_offset"]=plots.position_offset(plotting_df, title=theimg.imagename+" SUMSS Position Offset", base_filename=theimg.imagename.replace(".fits", ""),
            bmaj=theimg.bmaj*3600., bmin=theimg.bmin*3600., pa=theimg.bpa)
        theimg.plots["source_counts"]=plots.source_counts(sumss_touse.df, askap_touse.df, sumss_askap_crossmatch_diag.crossmatch_df, args.diagnostic_max_separation, 
            title=theimg.imagename+" source counts", base_filename=theimg.imagename.replace(".fits", ""))
        theimg.plots["flux_ratios_from_centre"]=plots.flux_ratios_distance_from_centre(plotting_df, args.diagnostic_max_separation, 
            title=theimg.imagename+" ASKAP / SUMSS int. flux ratio vs. Distance from Image Centre", base_filename=theimg.imagename.replace(".fits", ""))
        theimg.plots["flux_ratios"]=plots.flux_ratios_askap_flux(plotting_df, args.diagnostic_max_separation, 
            title=theimg.imagename+" ASKAP / SUMSS int. flux ratio vs. ASKAP Int Flux", base_filename=theimg.imagename.replace(".fits", ""))
            
        #Now create crossmatch for the purpose of transient searching (and a crossmatch before convolving if required)
        logger.info("Performing transient cross-match")
        if non_conv_askap_cat != None:
            sumss_askap_preconv_crossmatch_transient=crossmatch(sumss_catalog, non_conv_askap_cat)
            sumss_askap_preconv_crossmatch_transient.perform_crossmatch(maxsep=args.transient_max_separation)
        else:
            sumss_askap_preconv_crossmatch_transient = None
        sumss_askap_crossmatch_transient=crossmatch(sumss_catalog, askap_catalog)
        sumss_askap_crossmatch_transient.perform_crossmatch(maxsep=args.transient_max_separation)
        
        crossmatch_name=theimg.imagename.replace(".fits", "_sumss_askap_crossmatch_transient.csv")
        logger.info("Calculating flux ratios and separations of crossmatches.")
        sumss_askap_crossmatch_transient.calculate_ratio("askap_int_flux", "sumss_St", "askap_sumss_int_flux_ratio", col2_scaling=1.e-3)
        sumss_askap_crossmatch_transient.calculate_ratio("sumss_St", "askap_int_flux", "sumss_askap_int_flux_ratio", col1_scaling=1.e-3)
        sumss_askap_crossmatch_transient.calculate_diff("askap_ra", "sumss__RAJ2000", "askap_sumss_ra_offset")
        sumss_askap_crossmatch_transient.calculate_diff("askap_dec", "sumss__DEJ2000", "askap_sumss_dec_offset")
        
        #Write out crossmatch df to file
        sumss_askap_crossmatch_transient.write_crossmatch(crossmatch_name)
        
        
        if args.transients:
            sumss_askap_crossmatch_transient.transient_search(max_separation=args.transient_max_separation, askap_sumss_snr_thresh=askap_sumss_snr_thresh, large_flux_thresh=large_flux_ratio_thresh,
                pre_conv_crossmatch=sumss_askap_preconv_crossmatch_transient)
            os.makedirs("transients/no-match")
            os.makedirs("transients/large-ratio")
            os.makedirs("transients/askap-notseen")
            subprocess.call("mv transients*.csv transients/", shell=True)
        
        if args.postage_stamps:
            logger.info("Starting postage stamp production.")
            sumss_askap_crossmatch_transient.produce_postage_stamps(args.sumss_mosaic_dir,postage_options, nprocs=1, max_separation=args.transient_max_separation, convolve=args.convolve, pre_convolve_image=theimg.non_convolved)
            os.makedirs("postage-stamps/good")
            os.makedirs("postage-stamps/bad")
            if args.transients:
                try:
                    subprocess.check_output("mv transient*NOMATCH*_sidebyside.jpg transients/no-match/", shell=True, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as e:
                    logger.warning("No transient 'NO MATCH' images to move.")
                try:
                    subprocess.check_output("mv transient*LARGERATIO*_sidebyside.jpg transients/large-ratio/", shell=True, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as e:
                    logger.warning("No transient 'LARGE RATIO' images to move.")
                try:
                    subprocess.check_output("mv transient_askapnotseen*_sidebyside.jpg transients/askap-notseen/", shell=True, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as e:
                    logger.warning("No transient 'ASKAP NOT SEEN' images to move.")
                
            subprocess.call("mv *GOOD*_sidebyside.jpg postage-stamps/good/", shell=True)
            subprocess.call("mv *BAD*_sidebyside.jpg postage-stamps/bad/", shell=True)
        
        
        # Database Entry
        # image_id=theimg.inject_db(datestamp=launchtime)
        #Processing settings
        image_id=theimg.inject_db(datestamp=launchtime, user=username, description=args.db_tag, db_engine=args.db_engine, db_username=args.db_username, db_host=args.db_host, 
            db_port=args.db_port, db_database=args.db_database)
        theimg.inject_processing_db(image_id, full_output, askap_cat_file, sumss_source_cat, args.askap_ext_thresh, 
            args.sumss_ext_thresh, args.transient_max_separation, sf_sigmas, db_engine=args.db_engine, db_username=args.db_username, db_host=args.db_host, 
            db_port=args.db_port, db_database=args.db_database)
        sumss_askap_crossmatch_transient.inject_good_db(image_id, db_engine=args.db_engine, db_username=args.db_username, db_host=args.db_host, 
            db_port=args.db_port, db_database=args.db_database, max_separation=args.transient_max_separation)
        if args.transients:
            sumss_askap_crossmatch_transient.inject_transients_db(image_id, db_engine=args.db_engine, db_username=args.db_username, db_host=args.db_host, 
            db_port=args.db_port, db_database=args.db_database)
        
        if args.website_media_dir!="none":
            media_dir=os.path.join(args.website_media_dir, str(image_id))
        else:
            media_dir=os.path.join("..", "static", "media", "{}".format(image_id))
        stamp_media_dir=os.path.join(media_dir, "stamps")
        # os.makedirs(media_dir)
        os.makedirs(stamp_media_dir)
        subprocess.call("cp postage-stamps/good/*.jpg {}/".format(stamp_media_dir), shell=True)
        subprocess.call("cp postage-stamps/bad/*.jpg {}/".format(stamp_media_dir), shell=True)
        if args.transients:
            subprocess.call("cp transients/askap-notseen/*.jpg {}/".format(stamp_media_dir), shell=True)
        subprocess.call("cp *.png {}/".format(media_dir), shell=True)
             
        os.chdir("..")
        
        logger.info("Analysis complete for {}.".format(theimg.imagename))
    
    
    logger.info("All images done. Log will be placed in {}.".format(output))
    subprocess.call(["cp", logname+".log", output])
    subprocess.call(["rm", logname+".log"])
        
        
        
        
        
        
        

    
    
        
        
    
    
    
    
