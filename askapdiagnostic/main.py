#!/usr/bin/env python

import argparse
import logging
from askapdiagnostic.tools import utils
from askapdiagnostic.accessors.fitsimage import askapimage
from askapdiagnostic.tools.crossmatch import crossmatch
from askapdiagnostic.tools.catalog import catalog
from askapdiagnostic.plotting import plots
import datetime
import os
import pandas as pd
# from multiprocessing import Pool


#Steps are:
# 1. Read ASKAP image and get centre and size.
# 2. Perform source finding on image / or read already performed crossmatching.
# 3. Fetch SUMSS from Vizier / or read in SUMSS already performed.
# 4. Mask SUMSS to leave only sources that should be in the image.
# 5. Perform cross-matching using above outputs.
# 6. Produce plots and perform other analysis.

    
    
def source_finding(askapimg, sf="aegean", save_diag_images=False):
    if sf=="aegean":
        aegean_sf_options={
            "cores":1,
            "maxsummits":5,
            "seedclip":5,
            "floodclip":4,
            "nocov":""
        }
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
    
    parser = argparse.ArgumentParser()
    parser.add_argument("images", type=str, nargs="+", help="Define the images to process")
    parser.add_argument("--fetch-sumss", action="store_true", help="Fetch the SUMSS catalogue for the image area from Vizier.")
    parser.add_argument("--boundary-value", type=str, help="Define whether the out-of-bounds value in the FITS is 'nan' or 'zero'.", default="nan")
    parser.add_argument("--crossmatch-base", type=str, help="Define the base catalogue in the cross matching.", default="sumss")
    # parser.add_argument("--dont-mask-sumss", action="store_true", help="Do not filter the SUMSS catalogue such that only sources that should be in the ASKAP image remain.")
    args = parser.parse_args()
    
    
    utils.setup_logging("askapdiagnostic-{}".format(launchtime), use_colorlog)
    
    for image in args.images:
        #First load in the image and read in the information ready.
        theimg = askapimage(image, readinfo=True)
        
        #Check if source finding has already been done, if not do it
        askap_cat_file=theimg.imagename.replace(".fits", "_comp.csv")
        if not os.path.isfile(askap_cat_file):
            source_finding(theimg)
            #Check that the source finding was successful
            if not os.path.isfile(askap_cat_file):
                logger.error("Something has gone wrong")
        else:
            logger.info("ASKAP source catalogue for image {} found:".format(theimg.imagename))
            logger.info(askap_cat_file)
            logger.info("This catalogue will be used and no source extraction will be performed.")

        #Load the askap catalog into a new catalog object
        askap_catalog=catalog(pd.read_csv(askap_cat_file, delimiter=",", engine="python"),
         ra_col="ra", dec_col="dec")
            
        #Now check for SUMSS file or fetch if none
        sumss_source_cat=theimg.imagename.replace(".fits", "_sumss_comp.csv")
        if not os.path.isfile(sumss_source_cat):
            sumss_cat_df=theimg.get_sumss_catalogue(boundary_value=args.boundary_value, write_ann=False)
            #Write the filtered catalog to disk
            theimg.write_sumss_sources()
            logger.info("Fetching SUMSS catalogue complete.")
        else:
            logger.info("SUMSS catalogue already found.")
            sumss_cat_df=pd.read_csv(sumss_source_cat, delimiter=",", engine="python")
        
        #Initialise the SUMSS catalogue object
        sumss_catalog=catalog(sumss_cat_df)
        
        #Start crossmatching
        logger.info("Performing cross-match")
        basecat=args.crossmatch_base
        logger.info("Base catalogue: {}".format(basecat))
        
        #Create new crossmatch object
        sumss_askap_crossmatch=crossmatch(sumss_catalog.crossmatch_catalog, askap_catalog.crossmatch_catalog)
        sumss_askap_crossmatch.perform_crossmatch()
        
        #Get the fluxes from each catalog
        sumss_catalog_fluxes=sumss_catalog.get_col("Sp")
        #Need to convert SUMSS from mJy
        sumss_catalog_fluxes=sumss_catalog_fluxes/1.e3
        
        askap_catalog_fluxes=askap_catalog.get_col("peak_flux")
        
        
        plots.flux_ratio_image_view(sumss_catalog_fluxes, askap_catalog_fluxes, 20., 
            sumss_askap_crossmatch.idx, sumss_askap_crossmatch.d2d, askap_catalog.crossmatch_catalog)
        
        
        
        
        
        
        

    
    
        
        
    
    
    
    
