#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import logging
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.lines import Line2D
import aplpy
from astropy.coordinates import SkyCoord
from astropy import units as u
import multiprocessing as mp
from functools import partial
import os
import sys
import pkg_resources
import subprocess
import pandas as pd

logger = logging.getLogger(__name__)

def _sumss_rms(dec):
    if dec>-50:
        return 0.002
    else:
        return 0.0012

def _find_nearest_image(ra, dec, centres, images_data):
    askap_target=SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    seps = askap_target.separation(centres)
    min_index = np.argmin(seps.deg)
    # print d2d_sumss.deg[min_index]
    # print min_index
    image = images_data.iloc[min_index]["image"]
    # print image
    # print askap_target.to_string('hmsdms')
    return image

def _plotinitial(panels, key, figure, image, vmin=-99, vmax=-99, total_num_panels=2):
    panels[key]=aplpy.FITSFigure(image, figure=figure, subplot=(1,total_num_panels,key+1), slices=[0,0])
    if vmin == -99 and vmax == -99:
        panels[key].show_grayscale()
    elif vmin != -99 and vmax != -99:
        panels[key].show_grayscale(vmin=vmin, vmax=vmax)
    elif vmin == -99:
        panels[key].show_grayscale(vmax=vmax)
    else:
        panels[key].show_grayscale(vmin=vmin)
            
        # panels[key].show_contour(images[n-1], colors='red', levels=[3.*12e-6, 4.*12.e-6, 8.*12.e-6, 16*12.e-6])
    panels[key].set_theme('publication')
    return panels


def _copy_images_to_transient_folders(source_names, t_type, good_sources, bad_sources, softlink=False, move=False):
    if softlink:
        base_cmd="ln -s "
    elif move:
        base_cmd="mv "
    else:
        base_cmd="cp "
    for i in source_names:
        if i in good_sources:
            qual_tag="GOOD"
        else:
            qual_tag="BAD"
        thestamp="{}_{}_sidebyside.jpg".format(i, qual_tag)
        new_image_name="transient_{}_{}_{}_sidebyside.jpg".format(t_type, i, qual_tag)
        cmd=base_cmd+"{} {}".format(thestamp, new_image_name)
        try:
            subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            pass
            # logger.error("Copying of {} for transients failed.".format(thestamp))


def create_askap_postage_stamps(askap_df,crossmatch_df,askap_image, nprocs, sumss_mosaic_dir, nvss_mosaic_dir, cat_extractions, askap_extractions, cat2_extractions=[], radius=13./60., 
    convolve=False, askap_pre_convolve_image=None, askap_pre_convolve_catalog=None, dualmode=False, basecat="sumss", nvss_mosaic_data={}, nvss_centres=[]):
    #We need to load the SUMSS data
    #Load NVSS if needed
    # if dualmode or basecat=="sumss":
    try:
        logger.info("Loading SUMSS image data.")
        sumss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/sumss_images_info.csv"))
        sumss_centres = SkyCoord(ra=sumss_mosaic_data["center-ra"].values*u.deg, dec=sumss_mosaic_data["center-dec"].values*u.deg)
    except:
        logger.error("SUMSS mosaic data cannot be found!")
        return
        # load up the SUMSS image centre coords:
    
    if len(nvss_centres) == 0:
        try:
            logger.info("Loading NVSS image data.")
            nvss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/nvss_images_info.csv"))
            nvss_centres = SkyCoord(ra=nvss_mosaic_data["center-ra"].values*u.deg, dec=nvss_mosaic_data["center-dec"].values*u.deg)
        except:
            logger.error("NVSS mosaic data cannot be found!")
        
    matching_images=[]
    cat_fluxes=[]
    survey_used=[]
    
    for i,row in askap_df.iterrows():
        askap_target=SkyCoord(ra=row["ra"]*u.deg, dec=row["dec"]*u.deg)
        if askap_target.dec.degree <-30:
            search_centres=sumss_centres
            search_data=sumss_mosaic_data
            mosaic_dir = sumss_mosaic_dir
            survey_used.append("sumss")
        else:
            search_centres=nvss_centres
            search_data=nvss_mosaic_data
            mosaic_dir = nvss_mosaic_dir
            survey_used.append("nvss")
        seps = askap_target.separation(search_centres)
        min_index = np.argmin(seps.deg)
        # print d2d_sumss.deg[min_index]
        # print min_index
        image = search_data.iloc[min_index]["image"]
        # print image
        # print askap_target.to_string('hmsdms')
        matching_images.append(image.replace(".FITS", ""))
        #Open the image using aplpy as it has simple world2pixel function
        sumssfits = aplpy.FITSFigure(os.path.join(mosaic_dir, image), slices=[0,0])
        source_pixel_loc = sumssfits.world2pixel(askap_target.ra.deg, askap_target.dec.deg)
        source_pixel_loc_x = int(source_pixel_loc[0])
        source_pixel_loc_y = int(source_pixel_loc[1])
        temp_flux_values = []
        try:
            for i in [[0,0], [0, 1], [0,-1], [1,1], [1,0], [1,-1], [-1,-1], [-1, 0], [-1,1]]:
                temp_flux_values.append(sumssfits._data[source_pixel_loc_y+i[1], source_pixel_loc_x+i[0]])
        except:
            logger.error("Measuring flux failed for {}. Setting max flux to 0.".format(askap_df["name"]))
            temp_flux_values=[0,]
        peak_flux = np.max(temp_flux_values)
        cat_fluxes.append(peak_flux)
        sumssfits.close()
        
    if dualmode or basecat == "sumss":
        askap_df["sumss_Mosaic"]=matching_images
    else:
        askap_df["nvss_Mosaic"]=matching_images
    askap_df["cat_peak_flux"]=cat_fluxes
    askap_df["survey_used"]=survey_used
    
    if nprocs > 1:
        sumss_fits_mosaics = askap_df["sumss_Mosaic"].unique()
    
        sumss_looping_dict = _get_stamp_looping_parameters(askap_df, sumss_fits_mosaics, askap=True)
    
        postage_stamps_multi=partial(produce_comp_postage_stamps_multicore, askap_image=askap_image, params_dict=sumss_looping_dict, df=crossmatch_df, good_sources=[], bad_sources=[], 
            sumss_mosaic_dir=sumss_mosaic_dir, radius=radius, askap_only=True)
        askap_workers=mp.Pool(processes=nprocs)
        try:
            askap_workers.map(postage_stamps_multi, sumss_fits_mosaics)
        except KeyboardInterrupt:
            logger.warning("Caught KeyboardInterrupt, terminating jobs...")
            askap_workers.terminate()
            logger.info("Exiting...")
            askap_workers.close()
            sys.exit()
        
        askap_workers.close()
    else:
        produce_postage_stamps(askap_df, crossmatch_df, askap_image, sumss_mosaic_dir, nvss_mosaic_dir, [],[], cat_extractions, askap_extractions, cat2_extractions=cat2_extractions, radius=16./60., 
            max_separation=None, askap_only=True, convolve=convolve, askap_pre_convolve_image=askap_pre_convolve_image, askap_pre_convolve_catalog=askap_pre_convolve_catalog, 
            dualmode=dualmode, basecat=basecat)
    

def _get_stamp_looping_parameters(df, sumss_fits_mosaics, askap=False):
    sumss_looping_dict={}
    for sumss_img in sumss_fits_mosaics:
        sumss_looping_dict[sumss_img]=[]
        temp_df=df[df["sumss_Mosaic"]==sumss_img]
        for i,row in temp_df.iterrows():
            if askap==False:
                sumss_looping_dict[sumss_img].append([row["askap_ra"], row["askap_dec"], row["d2d"], row["askap_sumss_int_flux_ratio"], row["sumss_name"], row["askap_name"],
                    row["sumss__RAJ2000"], row["sumss__DEJ2000"]])
            else:
                sumss_looping_dict[sumss_img].append([row["ra"], row["dec"], "N/A", "N/A", "N/A", row["name"],
                    row["ra"], row["dec"]])
            if len(sumss_looping_dict[sumss_img])==5:
                break
    return sumss_looping_dict

def crossmatch_stamps(crossmatch, askap_image, postage_options, nprocs,sumss_mosaic_dir, nvss_mosaic_dir, radius=13./60., max_separation=15., convolve=False,
            askap_pre_convolve_image=None, askap_pre_convolve_catalog=None, dualmode=False, basecat="sumss", transients=False):
    # For reference, the transient df's are:
    # self.transients_no_matches_df=no_matches - crossmatch format
    # self.transients_large_ratios_df=large_ratios - crossmatch format
    # self.transients_not_matched_askap_should_see_df=not_matched_askap_sources_should_see - askap_cat -format
    
    #work out good and bad
    if transients:
        good_df=crossmatch.goodmatches_df_trans
        bad_df=crossmatch.badmatches_df_trans
    else:
        good_df=crossmatch.crossmatch_df[crossmatch.crossmatch_df["d2d"]<=max_separation].reset_index(drop=True)
        bad_df=crossmatch.crossmatch_df[crossmatch.crossmatch_df["d2d"]>max_separation].reset_index(drop=True)
    
    good_sources=good_df["master_name"].values
    bad_sources=bad_df["master_name"].values
    
    if "all" in postage_options:
        mode="all"
        filter_names=[]
        df = crossmatch.crossmatch_df
    elif postage_options==["transients",]:
        mode="transients"
        filter_names = crossmatch.transients_no_matches_df["master_name"].tolist() + crossmatch.transients_large_ratios_df["master_name"].tolist()
        df = crossmatch.crossmatch_df[crossmatch.crossmatch_df["master_name"].isin(filter_names)].reset_index(drop=True)
    
    
    
    cat_extractions=crossmatch.base_catalog.df
    if dualmode:
        cat2_extractions=crossmatch.merge_catalog.df 
    else:
        cat2_extractions = "None"
    askap_extractions=crossmatch.comp_catalog.df
    
    # We need to find the matching NVSS images.
    if dualmode or basecat=="nvss":
        try:
            logger.info("Loading NVSS image data.")
            nvss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/nvss_images_info.csv"))
            nvss_centres = SkyCoord(ra=nvss_mosaic_data["center-ra"].values*u.deg, dec=nvss_mosaic_data["center-dec"].values*u.deg)
        except:
            logger.error("NVSS mosaic data cannot be found!")
            return
        logger.info("Finding NVSS images")
        if dualmode:
            for i, row in df.iterrows():
                if row["survey_used"]=="nvss":
                    # print _find_nearest_image(row["master_ra"], row["master_dec"], nvss_centres, nvss_mosaic_data)
                    df.at[i, "sumss_Mosaic"] = _find_nearest_image(row["master_ra"], row["master_dec"], nvss_centres, nvss_mosaic_data)
                    # crossmatch.crossmatch_df.at[i, "sumss_Mosaic"] = _find_nearest_image(row["master_ra"], row["master_dec"], nvss_centres, nvss_mosaic_data)
        else:
            mosaics = []
            for i, row in df.iterrows():
                if row["survey_used"]=="nvss":
                    mosaics.append(_find_nearest_image(row["master_ra"], row["master_dec"], nvss_centres, nvss_mosaic_data))
            
            df["nvss_Mosaic"] = mosaics
            
    else:
        nvss_mosaic_data={}
        nvss_centres=[]
        
    
    if nprocs >1:
        #Get the unique SUMSS images.
        sumss_fits_mosaics = df["sumss_Mosaic"].unique()
        logger.info("{} unique SUMSS images needed for postage stamps".format(len(sumss_fits_mosaics)))
    
        sumss_looping_dict=_get_stamp_looping_parameters(df, sumss_fits_mosaics)
        #below is broken
        postage_stamps_multi=partial(produce_comp_postage_stamps_multicore, askap_image=askap_image, params_dict=sumss_looping_dict, df=df, 
            sumss_mosaic_dir=sumss_mosaic_dir, good_sources=good_sources, bad_sources=bad_sources, radius=radius)
        workers=mp.Pool(processes=nprocs)
        try:
            workers.map(postage_stamps_multi, sumss_fits_mosaics)
        except KeyboardInterrupt:
            logger.warning("Caught KeyboardInterrupt, terminating jobs...")
            workers.terminate()
            logger.info("Exiting...")
            workers.close()
            sys.exit()
    else:
        produce_postage_stamps(df, crossmatch.crossmatch_df, askap_image, 
            sumss_mosaic_dir, nvss_mosaic_dir, good_sources,bad_sources, cat_extractions, askap_extractions, cat2_extractions=cat2_extractions, radius=13./60., max_separation=None, askap_only=False, 
            convolve=convolve, askap_pre_convolve_image=askap_pre_convolve_image, askap_pre_convolve_catalog=askap_pre_convolve_catalog, dualmode=dualmode, 
            basecat=basecat, filter_names=filter_names)
        
    if "transients" in postage_options:
        _copy_images_to_transient_folders(crossmatch.transients_no_matches_df["master_name"].tolist(), "NOMATCH", good_sources, bad_sources)
        _copy_images_to_transient_folders(crossmatch.transients_large_ratios_df["master_name"].tolist(), "LARGERATIO", good_sources, bad_sources)
        # else:
        #     #all transients should have been done now above
        #     sumss_done=df["master_name"].tolist()
        #     df_required=crossmatch.crossmatch_df[~crossmatch.crossmatch_df.sumss_name.isin(sumss_done)].reset_index(drop=True)
        #     if nprocs >1:
        #         required_sumss_fits_mosaics = df_required["sumss_Mosaic"].unique()
        #         required_sumss_looping_dict=_get_stamp_looping_parameters(df_required, required_sumss_fits_mosaics)
        #         postage_stamps_multi=partial(produce_comp_postage_stamps_multicore, askap_image=askap_image, params_dict=required_sumss_looping_dict, df=df_required, sumss_mosaic_dir=sumss_mosaic_dir, radius=radius)
        #         try:
        #             workers.map(postage_stamps_multi, sumss_fits_mosaics)
        #         except KeyboardInterrupt:
        #             logger.warning("Caught KeyboardInterrupt, terminating jobs...")
        #             workers.terminate()
        #             logger.info("Exiting...")
        #             workers.close()
        #             sys.exit()
        #         #ASKAP only ones
        #         workers.close()
        #     else:
        #         produce_postage_stamps(df_required, crossmatch.crossmatch_df, askap_image, sumss_mosaic_dir, good_sources, bad_sources, radius=16./60., max_separation=None, askap_only=False, convolve=convolve, askap_pre_convolve_image=askap_pre_convolve_image,
        #             askap_pre_convolve_catalog=askap_pre_convolve_catalog, dualmode=, )
        
        create_askap_postage_stamps(crossmatch.transients_not_matched_askap_should_see_df, crossmatch.crossmatch_df, askap_image,nprocs, sumss_mosaic_dir, nvss_mosaic_dir, 
            cat_extractions, askap_extractions, cat2_extractions=cat2_extractions, convolve=convolve, askap_pre_convolve_image=askap_pre_convolve_image,
            askap_pre_convolve_catalog=askap_pre_convolve_catalog, dualmode=dualmode, basecat=basecat, nvss_mosaic_data=nvss_mosaic_data, nvss_centres=nvss_centres)
        
            
            
def produce_comp_postage_stamps_multicore(sumss_image, askap_image, params_dict, df, sumss_mosaic_dir, good_sources, bad_sources,radius=13./60., askap_only=False):
    try:
        sumss_image_path=os.path.join(sumss_mosaic_dir, sumss_image+".FITS")
        panels={}
        fig = plt.figure(figsize=(12, 6))
        key=1
        panels[key]=aplpy.FITSFigure(askap_image, figure=fig, subplot=(1,2,key+1))
        panels[key].show_grayscale()
            # panels[key].show_contour(images[n-1], colors='red', levels=[3.*12e-6, 4.*12.e-6, 8.*12.e-6, 16*12.e-6])
        panels[key].set_theme('publication')
    
        panels[key].show_ellipses(df["askap_ra"],df["askap_dec"],df["askap_b"]/3600., 
            df["askap_a"]/3600., angle=df["askap_pa"], layer="ASKAP Sources", color="#1f77b4")
        panels[key].show_ellipses(df["sumss__RAJ2000"],df["sumss__DEJ2000"],df["sumss_MinAxis"]/3600., 
            df["sumss_MajAxis"]/3600., angle=df["sumss_PA"], layer="SUMSS Sources", color="#d62728")
        
        key=0
    
        panels[key]=aplpy.FITSFigure(sumss_image_path, figure=fig, subplot=(1,2,key+1))
        panels[key].show_grayscale()
            # panels[key].show_contour(images[n-1], colors='red', levels=[3.*12e-6, 4.*12.e-6, 8.*12.e-6, 16*12.e-6])
        panels[key].set_theme('publication')
    
        panels[key].show_ellipses(df["askap_ra"],df["askap_dec"],df["askap_b"]/3600., 
            df["askap_a"]/3600., angle=df["askap_pa"], layer="ASKAP Sources", color="#1f77b4")
        panels[key].show_ellipses(df["sumss__RAJ2000"],df["sumss__DEJ2000"],df["sumss_MinAxis"]/3600., 
            df["sumss_MajAxis"]/3600., angle=df["sumss_PA"], layer="SUMSS Sources", color="#d62728")
        
        panels[key].axis_labels.hide()
        panels[key].tick_labels.hide()
    
        toplot=params_dict[sumss_image]
    
        for params in toplot:

            if askap_only:
                askap_ra=params[0]
                askap_dec=params[1]
                d2d=params[2]
                flux_ratio=params[3]
                sumss_name=""
                askap_name=params[5]
                sumss_ra=params[6]
                sumss_dec=params[7]
            else:
                askap_ra=params[0]
                askap_dec=params[1]
                d2d=params[2]
                flux_ratio=params[3]
                sumss_name=params[4]
                askap_name=params[5]
                sumss_ra=params[6]
                sumss_dec=params[7]
        
            panels[1].set_title("ASKAP {}".format(askap_name))
            panels[0].set_title("SUMSS {}".format(sumss_name))
            #Centre each image on the ASKAP coordinates for clarity
            recentre_ra=askap_ra
            recentre_dec=askap_dec
            for p in panels:
                panels[p].recenter(recentre_ra, recentre_dec, radius)
                panels[p].show_circles([recentre_ra], [recentre_dec], 120./3600., color='C1', label="ASKAP source", layer="ASKAP Source")
                if not askap_only:
                    panels[p].show_circles([sumss_ra], [sumss_dec],120./3600., color='C9', label="SUMSS source", layer="SUMSS Source")
        
            # sep_text=plt.text(0.02, 0.02, "Distance Separation = {:.2f} arcsec".format(d2d), transform=plt.gcf().transFigure)
            # ratio_text=plt.text(0.8, 0.02, "Int. Flux Ratio ASKAP/SUMSS = {:.2f}".format(flux_ratio), transform=plt.gcf().transFigure)
        
            #Figure name
            # plt.title(row["sumss_name"])
        
            if askap_only:
                figname = "transient_askapnotseen_ASKAP_{}_sidebyside.jpg".format(askap_name)
                custom_lines = [Line2D([0], [0], color='#1f77b4'),
                                Line2D([0], [0], color='#d62728'),    
                                Line2D([0], [0], color='C1')]    
                panels[0]._ax1.legend(custom_lines, ["ASKAP Sources", "SUMSS Sources", "Matched ASKAP"])
            else:
                if sumss_name in good_sources:
                    tag="GOOD"
                else:
                    tag="BAD"
                figname = "SUMSS_{}_{}_sidebyside.jpg".format(sumss_name, tag)
                custom_lines = [Line2D([0], [0], color='#1f77b4'),
                                Line2D([0], [0], color='#d62728'),    
                                Line2D([0], [0], color='C1'),    
                                Line2D([0], [0], color='C9')]
                panels[1]._ax1.legend(custom_lines, ["ASKAP Sources", "SUMSS Sources", "Matched ASKAP", "Matched SUMSS"])
                                

            fig.savefig(figname, bbox_inches="tight")

            logger.info("Saved figure {}.".format(figname))
        
            panels[1].remove_layer("ASKAP Source")
            panels[0].remove_layer("SUMSS Source")

            # sep_text.set_visible(False)
            # ratio_text.set_visible(False)
        
            panels[0].set_title("")
            panels[1].set_title("")
        plt.close(fig)
    except KeyboardInterrupt:
        logger.error("Keyboard Interrupt caught, will terminate")
        plt.close(fig)
        return
    
    #I think this will clear the SUMSS one, must be a more specific way.
    
    
def produce_postage_stamps(df, full_df, askap_fits, sumss_mosaic_dir, nvss_mosaic_dir, good_sources, bad_sources, cat_extractions, askap_extractions, cat2_extractions="None", 
    radius=13./60., max_separation=None, askap_only=False, convolve=False, askap_pre_convolve_image=None, askap_pre_convolve_catalog=None, dualmode=False, basecat="sumss", filter_names=[]):
        logger.info("Estimated time to completion = {:.2f} hours".format(len(df.index)*6./3600.))
        #Minimise fits opening so first get a list of all the unique SUMSS fits files to be used
        
        # Get the image paths here
        if dualmode or basecat=="sumss":
            fits_col = "sumss_Mosaic"
        else:
            fits_col = "nvss_Mosaic"
        
        for i, row in df.iterrows():
            # print sumss_mosaic_dir, nvss_mosaic_dir, row[fits_col]
            if row["survey_used"]=="sumss":
                df.at[i, fits_col] = os.path.join(sumss_mosaic_dir, row[fits_col]+".FITS")
            else:
                df.at[i, fits_col] = os.path.join(nvss_mosaic_dir, row[fits_col])
            
                    
        fits_mosaics = df[fits_col].unique()
        #For now support one ASKAP image at a time so can just take this from the first entry.
        
        #First initialise the figure and set up the askap panel with the image.
        if convolve:
            total_panels = 3
        else:
            total_panels = 2
        panels={}
        other={"sumss":"nvss", "nvss":"sumss"}
        fig = plt.figure(figsize=(18, 8))
        panels=_plotinitial(panels, 1, fig, askap_fits, total_num_panels=total_panels)#, vmin=-0.0025, vmax=0.007829)
        if convolve:
            panels=_plotinitial(panels, 2, fig, askap_pre_convolve_image, vmin=-0.0025, vmax=0.007829, total_num_panels=total_panels)
        
        #Load the ellipses onto the askap image
        for i in range(1,total_panels):
            panels[i].show_ellipses(askap_extractions["ra"],askap_extractions["dec"],askap_extractions["b"]/3600.*2.0, 
                askap_extractions["a"]/3600.*2.0, angle=askap_extractions["pa"], layer="ASKAP Sources", color="#1f77b4")
            panels[i].show_ellipses(cat_extractions["_RAJ2000"],cat_extractions["_DEJ2000"],cat_extractions["MinAxis"]/3600., 
                cat_extractions["MajAxis"]/3600., angle=cat_extractions["PA"], layer="{} Sources".format(basecat.upper()), color="#d62728")
            if dualmode:
                panels[i].show_ellipses(cat2_extractions["_RAJ2000"],cat2_extractions["_DEJ2000"],cat2_extractions["MinAxis"]/3600., 
                    cat2_extractions["MajAxis"]/3600., angle=cat2_extractions["PA"], layer="{} Sources".format(other[basecat].upper()), color="#F5B406")
            if convolve:
                panels[i].show_ellipses(askap_pre_convolve_catalog.df["ra"],askap_pre_convolve_catalog.df["dec"],askap_pre_convolve_catalog.df["b"]/3600.*2.0, 
                    askap_pre_convolve_catalog.df["a"]/3600.*2.0, angle=askap_pre_convolve_catalog.df["pa"], layer="ASKAP Non-Conv Sources", color="#11BA1C")
        
        bbox_dict=dict(boxstyle="round", ec="white", fc="white", alpha=0.7)
        
        #Now start the SUMSS loop
        for s_image in fits_mosaics:# [:1]:
            if s_image.split("/")[-1].startswith("J"):
                image_type="SUMSS"
            else:
                image_type="NVSS"
            #Get the actual fits file
            logger.debug("Image: {}".format(s_image))
            # s_image_path=os.path.join(sumss_mosaic_dir, s_image+".FITS")
            #Filter the dataframe such that only the SUMSS sources are present
            filtered_cross_matches=df[df[fits_col]==s_image].reset_index(drop=True)
            #Generate the base SUMSS panel
            if image_type=="SUMSS":
                vmin=-99.0
            else:
                vmin=-0.5e-2
                vmax=0.9e-2
            panels=_plotinitial(panels, 0, fig, s_image, total_num_panels=total_panels, vmin=vmin)
            # panels[1].set_title("SUMSS")
            #Add the sources
            panels[0].show_ellipses(askap_extractions["ra"],askap_extractions["dec"],askap_extractions["b"]/3600.*2.0, 
                askap_extractions["a"]/3600.*2.0, angle=askap_extractions["pa"], layer="ASKAP Sources", color="#1f77b4")
            panels[0].show_ellipses(cat_extractions["_RAJ2000"],cat_extractions["_DEJ2000"],cat_extractions["MinAxis"]/3600., 
                cat_extractions["MajAxis"]/3600., angle=cat_extractions["PA"], layer="{} Sources".format(basecat.upper()), color="#d62728")
            if dualmode:
                panels[0].show_ellipses(cat2_extractions["_RAJ2000"],cat2_extractions["_DEJ2000"],cat2_extractions["MinAxis"]/3600., 
                    cat2_extractions["MajAxis"]/3600., angle=cat2_extractions["PA"], layer="{} Sources".format(other[basecat].upper()), color="#F5B406")
            if convolve:
                panels[0].show_ellipses(askap_pre_convolve_catalog.df["ra"],askap_pre_convolve_catalog.df["dec"],askap_pre_convolve_catalog.df["b"]/3600.*2.0, 
                    askap_pre_convolve_catalog.df["a"]/3600.*2.0, angle=askap_pre_convolve_catalog.df["pa"], layer="ASKAP Non-Conv Sources", color="#11BA1C")
            for i in range(1,total_panels):
                panels[i].axis_labels.hide()
                panels[i].tick_labels.hide()
            #Now begin the main loop per source
            # debug_num=0
            for i, row in filtered_cross_matches.iterrows():
                if not askap_only:
                    panels[1].set_title("ASKAP "+row["askap_name"])
                    if convolve:
                        panels[2].set_title("ASKAP No Conv.")
                    panels[0].set_title(" "+row["master_name"])
                    recentre_ra=row["master_ra"]
                    recentre_dec=row["master_dec"]
                else:
                    panels[1].set_title("ASKAP "+row["name"])
                    panels[0].set_title(image_type)
                    if convolve:
                        panels[2].set_title("ASKAP No Conv.")
                    recentre_ra=row["ra"]
                    recentre_dec=row["dec"]
                #Centre each image on the ASKAP coordinates for clarity

                
                # try:
                for p in panels:
                    logger.debug("Fits Image: {}".format(s_image))
                    logger.debug("RA: {}".format(recentre_ra))
                    logger.debug("Dec: {}".format(recentre_dec))
                    logger.debug("Panel: {}".format(p))
                    panels[p].recenter(recentre_ra, recentre_dec, height=radius, width=radius)

                    if not askap_only:
                        #Edit here for transients
                        if row["master_name"] in bad_sources:
                            linestyle="--"
                        else:
                            linestyle="-"
                        panels[p].show_circles([row["askap_ra"]], [row["askap_dec"]], 120./3600., color='C1', linewidth=3, linestyle=linestyle, label="ASKAP source", layer="ASKAP Source")
                        panels[p].show_circles([recentre_ra], [recentre_dec], 120./3600., color='C9', linewidth=3, label="{} source".format(image_type), layer="Cat Source")
                    else:
                        panels[p].show_circles([recentre_ra], [recentre_dec], 120./3600., color='C1', linewidth=3, label="ASKAP position", layer="ASKAP Source")
                # except:
                    # logger.error("Can't zoom to region for source. Will skip.")
                    # logger.debug("SUMSS Image: {}".format(s_image))
                    # logger.debug("RA: {}".format(recentre_ra))
                    # logger.debug("Dec: {}".format(recentre_dec))
                    # continue
                
                if convolve:
                    pos_x_1 = 0.42
                    pos_y_1 = 0.68
                    pos_x_2 = 0.14
                    pos_y_2 = 0.7
                    pos_x_3 = 0.17
                    pos_y_3 = 0.15
                    pos_x_4 = 0.65
                    size=15
                else:
                    pos_x_1 = 0.57
                    pos_y_1 = 0.75
                    pos_x_2 = 0.14
                    pos_y_2 = 0.78
                    pos_x_3 = 0.02
                    pos_y_3 = 0.08
                    pos_x_4 = 0.8
                    size=12
                
                if not askap_only:
                    if row["master_name"] in good_sources:
                        sep_text=plt.text(pos_x_3, pos_y_3, "Distance Separation = {:.2f} arcsec".format(row["d2d"]), transform=plt.gcf().transFigure, size=size)
                        ratio_text=plt.text(pos_x_4, pos_y_3, "Int. Flux Ratio ASKAP/{} = {:.2f}".format(image_type, row["askap_cat_int_flux_ratio"]), transform=plt.gcf().transFigure, size=size)
                    if row["master_name"] in good_sources:
                        askap_snr_text=plt.text(pos_x_1, pos_y_1, "ASKAP Flux = {:.2f} mJy\nASKAP SNR = {:.2f}\n{} SNR = {:.2f}".format(row["askap_int_flux"]*1.e3, row["askap_askap_snr"],image_type, row["master_cat_snr"]), transform=plt.gcf().transFigure, bbox=bbox_dict)
                    if image_type=="SUMSS":
                        snr_col = "sumss_sumss_snr"
                    else:
                        snr_col = "nvss_nvss_snr"
                    sumss_snr_text=plt.text(pos_x_2, pos_y_2, "{0} Flux = {1:.2f} mJy\n{0} SNR ~ {2:.2f}".format(image_type, row["master_catflux"], row[snr_col]), transform=plt.gcf().transFigure, bbox=bbox_dict)
                else:
                    askap_snr_text=plt.text(pos_x_1, pos_y_1, "ASKAP Flux = {:.2f} mJy\nASKAP SNR = {:.2f}\n{} SNR = {}".format(row["int_flux"]*1.e3, row["askap_snr"], image_type, row["{}_snr".format(image_type.lower())]), transform=plt.gcf().transFigure, bbox=bbox_dict)
                    if image_type == "SUMSS":
                        cat_rms = _sumss_rms(row["dec"])
                    else:
                        cat_rms = 0.0005
                    sumss_snr_text=plt.text(pos_x_2, pos_y_2, "{0} Flux = {1:.2f} mJy\n{0} SNR ~ {2:.2f}".format(image_type, row["cat_peak_flux"]*1.e3, row["cat_peak_flux"]/cat_rms), transform=plt.gcf().transFigure, bbox=bbox_dict)
                
                #Figure name
                # plt.title(row["sumss_name"])
                if askap_only:
                    figname = "transient_askapnotseen_ASKAP_{}_sidebyside.jpg".format(row["name"])
                    custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                                    Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728')]  
                    if dualmode:
                        labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "{} Sources".format(other[basecat].upper()), "ASKAP Source Position"]
                        custom_lines+=[Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#F5B406'),]
                    else:
                        labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "ASKAP Source Position"]
                    custom_lines+=[Line2D([0], [0], color='C1'),]
                    if convolve:  
                        panels[2]._ax1.legend(custom_lines, labels)
                    else:
                        panels[1]._ax1.legend(custom_lines, labels)
                else:
                    if row["master_name"] in good_sources:
                        tag="GOOD"
                    else:
                        tag="BAD"
                    figname = "source_{}_{}_sidebyside.jpg".format(row["master_name"], tag)
                    if tag=="BAD":
                        askap_source_line=Line2D([0], [0], color='C1', linestyle=linestyle)
                    else:
                        askap_source_line=Line2D([0], [0], color='C1')
                    custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                                    Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),]   
                    if dualmode:
                        custom_lines += [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),]   
                                    
                    custom_lines+=[askap_source_line, Line2D([0], [0], color='C9')]
                                    
                    if convolve:
                        custom_lines.append(Line2D([0], [0], color="w", markeredgecolor='#11BA1C', marker="o", fillstyle="none", lw=2))
                    
                    labels=["ASKAP Sources", "{} Sources".format(basecat.upper())]
                    if dualmode:
                        labels+=["{} Source".format(other[basecat].upper())] 
                    
                    labels+=["Matched ASKAP", "Matched {}".format(image_type)]
                    if convolve:
                        panels[2]._ax1.legend(custom_lines,  labels+["Non-conv ASKAP Sources"])
                    else:
                        panels[1]._ax1.legend(custom_lines, labels)

                plt.savefig(figname, bbox_inches="tight")

                logger.info("Saved figure {}.".format(figname))
                
                for i in panels:                    
                    panels[i].remove_layer("ASKAP Source")
                    # panels[1].remove_layer("ASKAP Source")
                    if not askap_only:
                        panels[i].remove_layer("Cat Source")
                        # panels[1].remove_layer("SUMSS Source")

                if not askap_only and row["master_name"] in good_sources:
                    askap_snr_text.set_visible(False)
                elif askap_only:
                    askap_snr_text.set_visible(False)
                sumss_snr_text.set_visible(False)
                if not askap_only:
                    sep_text.set_visible(False)
                    ratio_text.set_visible(False)
                # debug_num+=1
                # if debug_num==1:
                #     debug_num=0
                #     break
            
            #I think this will clear the SUMSS one, must be a more specific way.    
            plt.gca().remove()

        plt.close()
    
