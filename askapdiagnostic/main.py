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
    
def source_finding(askapimg, logger, sf="aegean", options=None, save_diag_images=False):
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
            config = ConfigParser.SafeConfigParser()
            config.optionxform = str
            config.read([options])
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
            
        askapimg.find_sources(options=aegean_sf_options)
        if save_diag_images:
            aegean_sf_options["save"]=""
            askapimg.find_sources(options=aegean_sf_options)

def main():
    launchtime=datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    logger = logging.getLogger(__name__)

    try:
        import colorlog
        use_colorlog=True
    except ImportError:
        use_colorlog=False
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("images", type=str, nargs="+", help="Define the images to process")
    parser.add_argument("--output-tag", type=str, help="Add a tag to the output name.", default="")
    parser.add_argument("--log-level", type=str, choices=["WARNING", "INFO", "DEBUG"], help="Set the logging level.", default="INFO")
    parser.add_argument("--clobber", action="store_true", help="Overwrite output if already exists.")
    # parser.add_argument("--fetch-sumss", action="store_true", help="Fetch the SUMSS catalogue for the image area from Vizier.")
    parser.add_argument("--askap-csv", type=str, help="Manually define a aegean csv file containing the extracted sources to use for the ASKAP image.")
    parser.add_argument("--sumss-csv", type=str, help="Manually provide the SUMSS catalog csv.")
    parser.add_argument("--askap-csv-format", type=str, choices=["aegean"], help="Define which source finder provided the ASKAP catalog (currently only supports aegean).", default="aegean")
    parser.add_argument("--remove-extended", action="store_true", 
        help="Remove perceived extended sources from the catalogues. Uses the following arguments 'askap-ext-thresh' and 'sumss-ext-thresh' to set the threshold.")
    parser.add_argument("--askap-ext-thresh", type=float, 
        help="Define the maximum scaling threshold of the size of the ASKAP source compared to the PSF. Used to exclude extended sources. Only 1 axis has to exceed.", default=1.2)
    parser.add_argument("--sumss-ext-thresh", type=float, 
        help="Define the maximum scaling threshold of the size of the SUMSS source compared to the PSF. Use to exclude extended sources. Only 1 axis has to exceed.", default=1.2)
    #Option below provides framework for future source finder support, only aegaen for now.
    parser.add_argument("--use-all-fits", action="store_true", help="Use all the fits from Aegean ignoring all flags. Default only those with flag '0' are used.")
    parser.add_argument("--write-ann", action="store_true", help="Create kvis annotation files of the catalogues.")
    parser.add_argument("--boundary-value", type=str, choices=["nan", "zero"], help="Define whether the out-of-bounds value in the ASKAP FITS is 'nan' or 'zero'.", default="nan")
    parser.add_argument("--crossmatch-base", type=str, help="Define the base catalogue in the cross matching (currently not supported).", default="sumss")
    parser.add_argument("--max-separation", type=float, help="Maximum crossmatch distance (in arcsec) to be consdiered when creating plots.", default=None)
    parser.add_argument("--create-postage-stamps", action="store_true", help="Produce postage stamp plots of the cross matched sources within the max separation.")
    parser.add_argument("--sumss-mosaic-dir", type=str, help="Directory containing the SUMSS survey mosaic image files.", default=None)
    parser.add_argument("--aegean-settings-config", type=str, help="Select a config file containing the Aegean settings to be used (instead of defaults if none provided).", default=None)
    parser.add_argument("--transients", action="store_true", help="Perform a transient search analysis using the crossmatch data. Requires '--max-separation' to be defined.", default=None)
    # parser.add_argument("--dont-mask-sumss", action="store_true", help="Do not filter the SUMSS catalogue such that only sources that should be in the ASKAP image remain.")
    args = parser.parse_args()
    
    sf="aegean"  #For now hardcode the source finder to be aegean.
    
    logname="askapdiagnostic-{}".format(launchtime)
    utils.setup_logging(logname, args.log_level, use_colorlog)
    
    if args.create_postage_stamps:
        if args.sumss_mosaic_dir == None:
            logger.error("SUMSS mosaic directory has not been defined for generating the postage stamps!")
            logger.error("Define the SUMSS mosaic directory using '--sumss-mosaic-dir' and run again.")
            exit(logger)

    if args.transients:
        if args.max_separation==None:
            logger.error("Transients option selected but 'max-separation' is not defined.")
            logger.error("Define the max-separation using '--max-separation' and run again.")
            exit(logger)
        # if args.create_postage_stamps!=True:
        #     logger.warning("Transients option selected - turning on postage stamp production.")
        #     args.create_postage_stamps=True
    
    for image in args.images:
        logger.info("Beginning processing of {}".format(image))
        #First check that image exists and get abs path
        image = os.path.abspath(image)
        if not utils.checkfile(image):
            exit(logger)
        
        #Load in the image and read in the information ready.
        theimg = askapimage(image, readinfo=True)
        
        #Create output name and check if output already exists (as now there is easy access to the image name).
        output="{}".format(theimg.imagename.replace(".fits", "_results"))
        if args.output_tag!="":
            output=output.replace("_results", "_results_{}".format(args.output_tag))
            
        if not utils.createdir(output, clobber=args.clobber):
            exit(logger)
        
        #Check if catalogues have been provided
        if args.askap_csv != None:
            askap_cat_file = os.path.abspath(args.askap_csv)
            if not utils.checkfile(askap_cat_file):
                exit(logger)
            source_find=False
        else:
            logger.info("No ASKAP catalog provided - will perform source finding on the image using {}.".format(sf))
            source_find=True
            
        if args.sumss_csv != None:
            sumss_source_cat=os.path.abspath(args.sumss_csv)
            if not utils.checkfile(sumss_source_cat):
                exit(logger)
            fetch_sumss=False
        else:
            logger.info("No SUMSS catalog provided - will perform fetch of SUMSS sources.")
            fetch_sumss=True
            
        #If aegean settings provided need to check that the file is present
        if args.aegean_settings_config != None:
            aegean_settings_config=os.path.abspath(args.aegean_settings_config)
            if not utils.checkfile(aegean_settings_config):
                exit(logger)
            subprocess.call(["cp", aegean_settings_config, os.path.join(output, "aegean_settings_used.cfg")])
        else:
            aegean_settings_config = None
            
        #Now change directory
        os.chdir(output)
        
        #Perform source finding if needed
        if source_find:
            askap_cat_file=theimg.imagename.replace(".fits", "_comp.csv")
            source_finding(theimg, logger, options=aegean_settings_config)
            #Check that the source finding was successful
            if not os.path.isfile(askap_cat_file):
                logger.error("Source finding failed.")
                exit(logger)
        else:
            logger.info("Using provided catalog for {}:".format(theimg.imagename))
            logger.info(askap_cat_file)
            subprocess.call(["cp", askap_cat_file, "."])

        #Load the askap catalog into a new catalog object (but load it first to filter good sources)
        if sf=="aegean":
            askap_catalog=pd.read_csv(askap_cat_file, delimiter=",", engine="python")
            aegean_total = len(askap_catalog.index)
            logger.info("Total number of sources found by Aegean: {}".format(aegean_total))
            if not args.use_all_fits:
                askap_catalog=askap_catalog[askap_catalog["flags"]==0].reset_index(drop=True)
                aegean_good_fits=len(askap_catalog.index)
                logger.info("Total number of good fit sources found by Aegean: {} ({} sources removed).".format(aegean_good_fits, aegean_total - aegean_good_fits))
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
        
        #Filter extended sources if requested
        if args.remove_extended:
            noext_sumss_catalog=Catalog(sumss_catalog.remove_extended(threshold=args.sumss_ext_thresh,sumss_psf=True),
                 "{}".format(theimg.imagename.replace(".fits", "_sumss_noext")))
            noext_askap_catalog=Catalog(askap_catalog.remove_extended(threshold=args.askap_ext_thresh,ellipse_a="a", ellipse_b="b", 
                 beam_a=askap_catalog.df["psf_a"].iloc[0], beam_b=askap_catalog.df["psf_b"].iloc[0]),
                 "{}".format(theimg.imagename.replace(".fits", "_askap_noext")),ra_col="ra", dec_col="dec")

            sumss_touse = noext_sumss_catalog
            askap_touse = noext_askap_catalog
        else:
            sumss_touse = sumss_catalog
            askap_touse = askap_catalog
        
        #Add the respective image to the ASKAP catalog for later
        askap_touse.add_single_val_col("image", theimg.image)
        
        #Calculate distance from centre for each catalog for later
        askap_touse.add_distance_from_pos(theimg.centre)
        sumss_touse.add_distance_from_pos(theimg.centre)
        
        #Add SUMSS S/N for ASKAP sources
        askap_touse.add_sumss_sn(flux_col="int_flux")
        
        #Annotation files
        if args.write_ann:
            if fetch_sumss:
                raw_sumss.write_ann(color="RED")
            sumss_touse.write_ann()
            askap_touse.write_ann(color="BLUE", ellipse_a="a", ellipse_b="b", ellipse_pa="pa")
        
        #Start crossmatching
        logger.info("Performing cross-match")
        basecat=args.crossmatch_base
        logger.info("Base catalogue: {}".format(basecat))
        
        #Create new crossmatch object
        sumss_askap_crossmatch=crossmatch(sumss_touse, askap_touse)
        sumss_askap_crossmatch.perform_crossmatch()
        
        crossmatch_name=theimg.imagename.replace(".fits", "_sumss_askap_crossmatch.csv")
        
        #Calculate flux ratios and RA and Dec differences for plots
        logger.info("Calculating flux ratios and separations of crossmatches.")
        sumss_askap_crossmatch.calculate_ratio("askap_int_flux", "sumss_St", "askap_sumss_int_flux_ratio", col2_scaling=1.e-3)
        sumss_askap_crossmatch.calculate_ratio("sumss_St", "askap_int_flux", "sumss_askap_int_flux_ratio", col1_scaling=1.e-3)
        sumss_askap_crossmatch.calculate_diff("askap_ra", "sumss__RAJ2000", "askap_sumss_ra_offset")
        sumss_askap_crossmatch.calculate_diff("askap_dec", "sumss__DEJ2000", "askap_sumss_dec_offset")
        
        #Write out crossmatch df to file
        sumss_askap_crossmatch.write_crossmatch(crossmatch_name)
        
        #Get acceptable seperation df for plots
        if args.max_separation!=None:
            plotting_df = sumss_askap_crossmatch.crossmatch_df[sumss_askap_crossmatch.crossmatch_df["d2d"]<= args.max_separation]
        else:
            plotting_df = sumss_askap_crossmatch.crossmatch_df
            
        plotting_len=len(plotting_df.index)
        logger.debug("Length of plotting_df: {}".format(plotting_len))
        
        #Plots    
        logger.info("Producing plots.")
        plots.flux_ratio_image_view(plotting_df, title=theimg.imagename+" ASKAP / SUMSS flux ratio", base_filename=theimg.imagename.replace(".fits", ""))
        plots.position_offset(plotting_df, title=theimg.imagename+" SUMSS Position Offset", base_filename=theimg.imagename.replace(".fits", ""))
        plots.source_counts(sumss_touse.df, askap_touse.df, sumss_askap_crossmatch.crossmatch_df, args.max_separation, 
            title=theimg.imagename+" source counts", base_filename=theimg.imagename.replace(".fits", ""))
        plots.flux_ratios_distance_from_centre(plotting_df, args.max_separation, 
            title=theimg.imagename+" ASKAP / SUMSS int. flux ratio vs. Distance from Image Centre", base_filename=theimg.imagename.replace(".fits", ""))
        plots.flux_ratios_askap_flux(plotting_df, args.max_separation, 
            title=theimg.imagename+" ASKAP / SUMSS int. flux ratio vs. ASKAP Int Flux", base_filename=theimg.imagename.replace(".fits", ""))
        
        if args.transients:
            sumss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "data/sumss_images_info.csv"))
            sumss_askap_crossmatch.transient_search(max_separation=args.max_separation)
            os.mkdir("transients")
            subprocess.call("mv transients*.csv transients/", shell=True)
        
        if args.create_postage_stamps:
            logger.info("Beginning postage stamp production.")
            sumss_askap_crossmatch.produce_postage_stamps(args.sumss_mosaic_dir)# , max_separation=args.max_separation) for now hard lock all matches to be postage stamped.
            os.mkdir("postage-stamps")
            subprocess.call("mv *_sidebyside.png postage-stamps", shell=True)
            
        os.chdir("..")
        
        logger.info("Analysis complete for {}.".format(theimg.imagename))
        
    logger.info("All images done. Log will be placed in {}.".format(output))
    subprocess.call(["cp", logname+".log", output])
    subprocess.call(["rm", logname+".log"])
        
        
        
        
        
        
        

    
    
        
        
    
    
    
    
