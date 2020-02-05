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
    # panels[key].show_colorbar()
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


def create_askap_postage_stamps(askap_df,crossmatch_df,askap_image, nprocs, cat_extractions, askap_extractions, cat2_extractions=[], radius=13./60., 
    convolve=False, askap_pre_convolve_image=None, askap_pre_convolve_catalog=None, dualmode=False, basecat="sumss"):
    #We need to load the SUMSS data
    #Load NVSS if needed
    # if dualmode or basecat=="sumss":
    # try:
    #     logger.info("Loading SUMSS image data.")
    #     sumss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/sumss_images_info.csv"))
    #     sumss_centres = SkyCoord(ra=sumss_mosaic_data["center-ra"].values*u.deg, dec=sumss_mosaic_data["center-dec"].values*u.deg)
    # except:
    #     logger.error("SUMSS mosaic data cannot be found!")
    #     return
    #     # load up the SUMSS image centre coords:
    #
    # if len(nvss_centres) == 0:
    #     try:
    #         logger.info("Loading NVSS image data.")
    #         nvss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/nvss_images_info.csv"))
    #         nvss_centres = SkyCoord(ra=nvss_mosaic_data["center-ra"].values*u.deg, dec=nvss_mosaic_data["center-dec"].values*u.deg)
    #     except:
    #         logger.error("NVSS mosaic data cannot be found!")
    #
    # matching_images=[]
    # cat_fluxes=[]
    # survey_used=[]
    
    # for i,row in askap_df.iterrows():
#         askap_target=SkyCoord(ra=row["ra"]*u.deg, dec=row["dec"]*u.deg)
#         if askap_target.dec.degree <-30:
#             search_centres=sumss_centres
#             search_data=sumss_mosaic_data
#             mosaic_dir = sumss_mosaic_dir
#             survey_used.append("sumss")
#         else:
#             search_centres=nvss_centres
#             search_data=nvss_mosaic_data
#             mosaic_dir = nvss_mosaic_dir
#             survey_used.append("nvss")
#         seps = askap_target.separation(search_centres)
#         min_index = np.argmin(seps.deg)
#         # print d2d_sumss.deg[min_index]
#         # print min_index
#         image = search_data.iloc[min_index]["image"]
#         # print image
#         # print askap_target.to_string('hmsdms')
#         matching_images.append(image.replace(".FITS", ""))
#         #Open the image using aplpy as it has simple world2pixel function
#         sumssfits = aplpy.FITSFigure(os.path.join(mosaic_dir, image), slices=[0,0])
#         source_pixel_loc = sumssfits.world2pixel(askap_target.ra.deg, askap_target.dec.deg)
#         source_pixel_loc_x = int(source_pixel_loc[0])
#         source_pixel_loc_y = int(source_pixel_loc[1])
#         temp_flux_values = []
#         try:
#             for i in [[0,0], [0, 1], [0,-1], [1,1], [1,0], [1,-1], [-1,-1], [-1, 0], [-1,1]]:
#                 temp_flux_values.append(sumssfits._data[source_pixel_loc_y+i[1], source_pixel_loc_x+i[0]])
#         except:
#             logger.error("Measuring flux failed for {}. Setting max flux to 0.".format(askap_df["name"]))
#             temp_flux_values=[0,]
#         peak_flux = np.max(temp_flux_values)
#         cat_fluxes.append(peak_flux)
#         sumssfits.close()
        
    # if dualmode or basecat == "sumss":
    #     askap_df["sumss_Mosaic"]=matching_images
    # else:
    #     askap_df["nvss_Mosaic"]=matching_images
    # askap_df["cat_peak_flux"]=cat_fluxes
    # askap_df["survey_used"]=survey_used
    
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
        produce_postage_stamps(askap_df, crossmatch_df, askap_image, [],[], cat_extractions, askap_extractions, cat2_extractions=cat2_extractions, radius=16./60., 
            max_separation=None, convolve=convolve, askap_pre_convolve_image=askap_pre_convolve_image, askap_pre_convolve_catalog=askap_pre_convolve_catalog, 
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

def _get_extracted_rectangle_size(panel, ra, dec, shift):
    pixel_x,pixel_y=panel.world2pixel([ra,], [dec,])
    #Currently shift by 5 so I'll hard fix this for now
    shifted_pixel_coords=[pixel_x[0]+shift, pixel_y[0]]
    new_ra1, new_dec1 = panel.pixel2world(shifted_pixel_coords[0], shifted_pixel_coords[1])
    ra_size = np.abs(ra-new_ra1)*2.
    shifted_pixel_coords=[pixel_x[0], pixel_y[0]+shift]
    new_ra2, new_dec2 = panel.pixel2world(shifted_pixel_coords[0], shifted_pixel_coords[1])
    dec_size = np.abs(dec-new_dec2)*2.
    return ra_size, dec_size
    

def crossmatch_stamps(crossmatch, askap_image, postage_options, nprocs, radius=13./60., max_separation=15., convolve=False,
            askap_pre_convolve_image=None, askap_pre_convolve_catalog=None, dualmode=False, basecat="sumss", transients=False, transient_min_ratio=2.0):

    #Here we just need to plot from the transients master df based on the ratio - good and bad no longer exists
    
    
    #work out good and bad
    if transients:
        #For now hardcode in limit
        # mean=crossmatch.transients_master_df["master_ratio"].mean()
        # std=crossmatch.transients_master_df["master_ratio"].std()
        # transient_min_ratio = mean + (1.0 *std)
        df_to_plot = crossmatch.transients_master_df[crossmatch.transients_master_df["master_ratio"]>=transient_min_ratio]
    else:
        df_to_plot = crossmatch.transients_master_df
    
    # good_sources=good_df["master_name"].values
    # bad_sources=bad_df["master_name"].values
    
    # if "all" in postage_options:
    #     mode="all"
    #     filter_names=[]
    #     df = crossmatch.crossmatch_df
    # elif postage_options==["transients",]:
    #     mode="transients"
    #     filter_names = crossmatch.transients_no_matches_df["master_name"].tolist() + crossmatch.transients_large_ratios_df["master_name"].tolist()
    #     # df = crossmatch.crossmatch_df[crossmatch.crossmatch_df["master_name"].isin(filter_names)].reset_index(drop=True)
    #     df = crossmatch.transients_no_matches_df.append(crossmatch.transients_large_ratios_df, ignore_index=True)
    
    
    
    cat_extractions=crossmatch.base_catalog.df
    if dualmode:
        cat2_extractions=crossmatch.merge_catalog.df 
    else:
        cat2_extractions = "None"
    askap_extractions=crossmatch.comp_catalog.df
    
    # We need to find the matching NVSS images. NOW ALREADY FOUND
    # if dualmode or basecat=="nvss":
    #     try:
    #         logger.info("Loading NVSS image data.")
    #         nvss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/nvss_images_info.csv"))
    #         nvss_centres = SkyCoord(ra=nvss_mosaic_data["center-ra"].values*u.deg, dec=nvss_mosaic_data["center-dec"].values*u.deg)
    #     except:
    #         logger.error("NVSS mosaic data cannot be found!")
    #         return
    #     logger.info("Finding NVSS images")
    #     if dualmode:
    #         for i, row in df.iterrows():
    #             if row["survey_used"]=="nvss":
    #                 # print _find_nearest_image(row["master_ra"], row["master_dec"], nvss_centres, nvss_mosaic_data)
    #                 df.at[i, "sumss_Mosaic"] = _find_nearest_image(row["master_ra"], row["master_dec"], nvss_centres, nvss_mosaic_data)
    #                 # crossmatch.crossmatch_df.at[i, "sumss_Mosaic"] = _find_nearest_image(row["master_ra"], row["master_dec"], nvss_centres, nvss_mosaic_data)
    #     else:
    #         mosaics = []
    #         for i, row in df.iterrows():
    #             if row["survey_used"]=="nvss":
    #                 mosaics.append(_find_nearest_image(row["master_ra"], row["master_dec"], nvss_centres, nvss_mosaic_data))
    #
    #         df["nvss_Mosaic"] = mosaics
    #
    # else:
    #     nvss_mosaic_data={}
    #     nvss_centres=[]
        
    
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
        produce_postage_stamps(df_to_plot, askap_image, 
            cat_extractions, askap_extractions, cat2_extractions=cat2_extractions, radius=13./60., max_separation=None, 
            convolve=convolve, askap_pre_convolve_image=askap_pre_convolve_image, askap_pre_convolve_catalog=askap_pre_convolve_catalog, dualmode=dualmode, 
            basecat=basecat)
        
    # if "transients" in postage_options:
        # _copy_images_to_transient_folders(crossmatch.transients_no_matches_df["master_name"].tolist(), "NOMATCH", good_sources, bad_sources)
        # _copy_images_to_transient_folders(crossmatch.transients_large_ratios_df["master_name"].tolist(), "LARGERATIO", good_sources, bad_sources)
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
        
        # create_askap_postage_stamps(crossmatch.transients_not_matched_askap_should_see_df, crossmatch.crossmatch_df, askap_image,nprocs,
        #     cat_extractions, askap_extractions, cat2_extractions=cat2_extractions, convolve=convolve, askap_pre_convolve_image=askap_pre_convolve_image,
        #     askap_pre_convolve_catalog=askap_pre_convolve_catalog, dualmode=dualmode, basecat=basecat)
        
            
            
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
    
    
def produce_postage_stamps(df, askap_fits, cat_extractions, askap_extractions, cat2_extractions="None", 
    radius=13./60., max_separation=None, convolve=False, askap_pre_convolve_image=None, askap_pre_convolve_catalog=None, dualmode=False, basecat="sumss"):
        # logger.info("Estimated time to completion = {:.2f} hours".format(len(df.index)*6./3600.))
        #Minimise fits opening so first get a list of all the unique SUMSS fits files to be used
        
        # Get the image paths here
        # if dualmode or basecat=="sumss":
        #     fits_col = "sumss_Mosaic"
        # else:
        #     fits_col = "nvss_Mosaic"
        #
        # for i, row in df.iterrows():
        #     # print sumss_mosaic_dir, nvss_mosaic_dir, row[fits_col]
        #     if row["survey_used"]=="sumss":
        #         df.at[i, fits_col] = os.path.join(sumss_mosaic_dir, row[fits_col]+".FITS")
        #     else:
        #         df.at[i, fits_col] = os.path.join(nvss_mosaic_dir, row[fits_col])
        # if askap_only:
        #     fits_col="catalog_Mosaic_path"
        # else:
        #     fits_col="master_catalog_Mosaic_path"
        fits_mosaics = df["master_catalog_mosaic_path"].unique()
        #For now support one ASKAP image at a time so can just take this from the first entry.
        
        #First initialise the figure and set up the askap panel with the image.
        if convolve:
            total_panels = 3
        else:
            total_panels = 2
        panels={}
        other={"sumss":"nvss", "nvss":"sumss"}
        fig = plt.figure(figsize=(18, 8))
        panels=_plotinitial(panels, 1, fig, askap_fits, total_num_panels=total_panels, vmin=-0.5e-02, vmax=1e-02)
        if convolve:
            panels=_plotinitial(panels, 2, fig, askap_pre_convolve_image, vmin=-2.e-3, vmax=2.5e-3, total_num_panels=total_panels)
        
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
        
        bbox_dict=dict(boxstyle="round", ec="white", fc="white", alpha=0.8)
        
        #Now start the Catalog FITS loop
        for s_image in fits_mosaics:# [:1]:
            if s_image.split("/")[-1] == "SKIP.FITS":
                logger.warning("Missing SUMSS FITS file: skipping these sources.")
                continue
            if s_image.split("/")[-1].startswith("J"):
                image_type="SUMSS"
            else:
                image_type="NVSS"
            #Get the actual fits file
            logger.debug("Image: {}".format(s_image))
            # s_image_path=os.path.join(sumss_mosaic_dir, s_image+".FITS")
            #Filter the dataframe such that only the SUMSS sources are present
            filtered_cross_matches=df[df["master_catalog_mosaic_path"]==s_image].reset_index(drop=True)
            # if not askap_only:
            if not pd.isna(filtered_cross_matches.iloc[0]["master_catalog_Mosaic_rms"]):
                image_rms = filtered_cross_matches.iloc[0]["master_catalog_Mosaic_rms"]
            else:
                image_rms = filtered_cross_matches.iloc[0]["catalog_Mosaic_rms"]
            # else:
                # image_rms = filtered_cross_matches.iloc[0]["catalog_Mosaic_rms"]
            #Generate the base SUMSS panel
            if image_type=="SUMSS":
                if int(s_image.split(".FITS")[0][-2:]) <= 48:
                    vmin=-0.7e-02
                    vmax=1.5e-02
                else:
                    vmin=-0.4e-02
                    vmax=0.9e-02
            else:
                vmin=-0.1e-2
                vmax=0.5e-2
            panels=_plotinitial(panels, 0, fig, s_image, total_num_panels=total_panels, vmin=vmin, vmax=vmax)
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
                skip=False
                if row["type"] == "goodmatch":
                    askap_title = "RACS "+row["askap_name"]
                    if pd.isna(row["aegean_convolved_int_flux"]):
                        askap_title_2 = "RACS "+row["askap_non_conv_name"]
                    else:
                        askap_title_2 = "RACS" + row["askap_name"]
                    sumss_title = "SUMSS "+row["master_name"]
                    recentre_ra=row["master_ra"]
                    recentre_dec=row["master_dec"]
                    
                elif row["type"] == "noaskapmatch":
                    askap_title = "RACS (Convolved)"
                    askap_title_2 = "RACS (Non-Convolved)"
                    sumss_title = "SUMSS "+row["master_name"]
                    recentre_ra=row["master_ra"]
                    recentre_dec=row["master_dec"]
                elif row["type"] == "nocatalogmatch":
                    askap_title = "RACS "+row["askap_name"]
                    askap_title_2 = "RACS "+row["askap_non_conv_name"]
                    sumss_title = "SUMSS ({})".format(row["master_catalog_mosaic"])
                    recentre_ra=row["askap_ra"]
                    recentre_dec=row["askap_dec"]
                
                panels[0].set_title(sumss_title)
                panels[1].set_title(askap_title)
                if convolve:
                    panels[2].set_title(askap_title_2)
               
                #Centre each image on the ASKAP coordinates for clarity
                if convolve:
                    pos_x_1 = 0.41
                    pos_y_1 = 0.68
                    pos_x_2 = 0.14
                    pos_y_2 = 0.27
                    pos_x_3 = 0.17
                    pos_y_3 = 0.15
                    pos_x_4 = 0.6
                    pos_x_5 = 0.68
                    pos_y_5 = 0.27
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
                
                snr_col = "sumss_sumss_snr"
                
                for p in panels:
                    logger.debug("Fits Image: {}".format(s_image))
                    logger.debug("RA: {}".format(recentre_ra))
                    logger.debug("Dec: {}".format(recentre_dec))
                    logger.debug("Panel: {}".format(p))
                    try:
                        panels[p].recenter(recentre_ra, recentre_dec, height=radius, width=radius)
                    except:
                        logger.error("Image out of range, can't recentre. Skipping Image.")
                        skip=True
                        break
                        

                    if row["type"] == "goodmatch":
                        panels[p].show_circles([recentre_ra], [recentre_dec], 120./3600., color='C9', linewidth=3, label="{} source".format(image_type), layer="Cat Source")
                        panels[p].show_circles([row["askap_ra"]], [row["askap_dec"]], 120./3600., color='C1', linewidth=3, label="ASKAP source", layer="ASKAP Source")
                        if p==0:
                            sep_text=plt.text(pos_x_3, pos_y_3, "Distance Separation = {:.2f} arcsec".format(row["d2d"]), transform=plt.gcf().transFigure, size=size)
                            ratio_text=plt.text(pos_x_4, pos_y_3, "Scaled Int. Flux Ratio = {:.2f} +/- {:.2f}".format(row["master_ratio"], row["master_ratio_err"]), transform=plt.gcf().transFigure, size=size)
                            askap_snr_text=plt.text(pos_x_1, pos_y_2, "ASKAP Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["askap_flux_to_use"]*1.e3, 
                                row["askap_flux_to_use_err"]*1.e3, row["askap_rms"]*1.e3, row["measured_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict)
                            sumss_snr_text=plt.text(pos_x_2, pos_y_2, "{0} Int. Flux = {1:.2f} +/- {2:.2f} mJy\n{0} RMS ~ {3:.2f} mJy\n{0} SNR = {4:.2f}".format(image_type, 
                                row["catalog_flux_to_use"]*1.e3,row["catalog_flux_to_use_err"]*1.e3,image_rms*1.e3, row[snr_col]), transform=plt.gcf().transFigure, bbox=bbox_dict)
                            if convolve:
                                askap_snr_text_preconv=plt.text(pos_x_5, pos_y_2, "ASKAP Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["askap_flux_to_use_2"]*1.e3, 
                                    row["askap_flux_to_use_2_err"]*1.e3, row["askap_rms_preconv"]*1.e3, row["measured_preconv_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict)
                                custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                                                Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),
                                                Line2D([0], [0], color="w", markeredgecolor='#11BA1C', marker="o", fillstyle="none", lw=2),
                                                Line2D([0], [0], color="w", markeredgecolor='C9', marker="o", fillstyle="none", lw=2),
                                                Line2D([0], [0], color='C1'),
                                                Line2D([0], [0], color='C9')
                                            ]
                                labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "Non-conv ASKAP sources", "Extracted Flux", "SUMSS Source Position", "ASKAP Source Position"]
                            else:
                                custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                                                Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),
                                                Line2D([0], [0], color="w", markeredgecolor='C9', marker="o", fillstyle="none", lw=2),
                                                Line2D([0], [0], color='C1'),
                                                Line2D([0], [0], color='C9')
                                            ]
                                labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "Extracted Flux", "SUMSS Source Position", "ASKAP Source Position"]
                    elif row["type"] == "noaskapmatch":
                        panels[p].show_circles([recentre_ra], [recentre_dec], 120./3600., color='C9', linewidth=3, label="{} source".format(image_type), layer="Cat Source")
                        if p==0:
                            ratio_text=plt.text(pos_x_4, pos_y_3, "Scaled Int. Flux Ratio = {:.2f} +/- {:.2f}".format(row["master_ratio"], row["master_ratio_err"]), transform=plt.gcf().transFigure, size=size)
                        if convolve:
                            if p==1:
                                panels[p].show_circles([recentre_ra], [recentre_dec], 22.5/3600., color='C1', linewidth=1, label="Extracted Flux", layer="Extracted Flux")
                            if p==2:
                                panels[p].show_circles([recentre_ra], [recentre_dec], 7./3600., color='C1', linewidth=1, label="Extracted Flux", layer="Extracted Flux_2")
                            if p==0:
                                askap_snr_text=plt.text(pos_x_1, pos_y_2, "ASKAP Measured Int. Flux = {:.2f} mJy +/- {:.2f} \nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["askap_flux_to_use"]*1.e3, 
                                    row["askap_flux_to_use_err"]*1.e3, row["askap_rms"]*1.e3, row["measured_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict)
                                askap_snr_text_preconv=plt.text(pos_x_5, pos_y_2, "ASKAP Measured Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["askap_flux_to_use_2"]*1.e3, 
                                    row["askap_flux_to_use_2_err"]*1.e3, row["askap_rms_preconv"]*1.e3, row["measured_preconv_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict) 
                                custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                                                Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),
                                                Line2D([0], [0], color="w", markeredgecolor='#11BA1C', marker="o", fillstyle="none", lw=2),
                                                Line2D([0], [0], color="w", markeredgecolor='C9', marker="o", fillstyle="none", lw=2),
                                                Line2D([0], [0], color='C1')
                                            ]
                                labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "Non-conv ASKAP sources", "Extracted Flux", "SUMSS Source Position"]
                        else:
                            if p==1:
                                panels[p].show_circles([recentre_ra], [recentre_dec], 7./3600., color='C1', linewidth=1, label="Extracted Flux", layer="Extracted Flux")
                            if p==0:
                                askap_snr_text=plt.text(pos_x_1, pos_y_2, "ASKAP Measured Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["askap_flux_to_use"]*1.e3, 
                                    row["askap_flux_to_use_err"]*1.e3, row["askap_rms"]*1.e3, row["measured_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict)
                                custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                                                Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),
                                                Line2D([0], [0], color="w", markeredgecolor='C9', marker="o", fillstyle="none", lw=2),
                                                Line2D([0], [0], color='C1')
                                            ]
                                labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "Extracted Flux", "SUMSS Source Position"]
                        if p==0:
                            sumss_snr_text=plt.text(pos_x_2, pos_y_2, "{0} Int. Flux = {1:.2f} +/- {2:.2f} mJy\n{0} RMS ~ {3:.2f} mJy\n{0} SNR = {4:.2f}".format(image_type, 
                                row["catalog_flux_to_use"]*1.e3, row["catalog_flux_to_use_err"]*1.e3,image_rms*1.e3, row[snr_col]), transform=plt.gcf().transFigure, bbox=bbox_dict)
                        #Rectangle on extracted flux on ASKAP images
                        # x_size,y_size=_get_extracted_rectangle_size(panels[1], recentre_ra, recentre_dec, 10)
                        # panels[1].show_rectangles([recentre_ra], [recentre_dec], x_size, y_size, color='C1', linewidth=1, label="Extracted Flux", layer="Extracted Flux")
                        # if convolve:
                        #     x_size,y_size=_get_extracted_rectangle_size(panels[2], recentre_ra, recentre_dec, 10)
                        #     panels[2].show_rectangles([recentre_ra], [recentre_dec], x_size, y_size, color='C1', linewidth=1, label="Extracted Flux", layer="Extracted Flux")
                        
                    elif row["type"] == "nocatalogmatch":
                        panels[p].show_circles([recentre_ra], [recentre_dec], 120./3600., color='C1', linewidth=3, label="ASKAP position", layer="ASKAP Source")
                        if p==0:
                            ratio_text=plt.text(pos_x_4, pos_y_3, "Scaled Int. Flux Ratio = {:.2f} +/- {:.2f}".format(row["master_ratio"], row["master_ratio_err"]), transform=plt.gcf().transFigure, size=size)
                            #Put a rectangle to show extracted peak flux measurement
                            # x_size,y_size=_get_extracted_rectangle_size(panels[0], recentre_ra, recentre_dec, 5)
                            panels[p].show_circles([recentre_ra], [recentre_dec], 22.5/3600., color='C9', linewidth=1, label="Extracted Flux", layer="Extracted Flux")
                            askap_snr_text=plt.text(pos_x_1, pos_y_2, "ASKAP Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["askap_flux_to_use"]*1.e3, 
                                row["askap_flux_to_use_err"]*1.e3,row["askap_rms"]*1.e3, row["measured_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict)
                            if convolve:
                                askap_snr_text_preconv=plt.text(pos_x_5, pos_y_2, "ASKAP Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["askap_flux_to_use_2"]*1.e3, 
                                    row["askap_flux_to_use_2_err"]*1.e3,row["askap_rms_preconv"]*1.e3, row["measured_preconv_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict) 
                                custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                                                Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),
                                                Line2D([0], [0], color="w", markeredgecolor='#11BA1C', marker="o", fillstyle="none", lw=2),
                                                Line2D([0], [0], color="w", markeredgecolor='C9', marker="o", fillstyle="none", lw=2),
                                                Line2D([0], [0], color='C9')
                                            ]
                                labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "Non-conv ASKAP sources", "Extracted Flux", "ASKAP Source Position"]
                            else:
                                custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                                                Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),
                                                Line2D([0], [0], color="w", markeredgecolor='C9', marker="o", fillstyle="none", lw=2),
                                                Line2D([0], [0], color='C9')
                                            ]
                                labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "Extracted Flux", "ASKAP Source Position"]
                            sumss_snr_text=plt.text(pos_x_2, pos_y_2, "Measured {0} Int. Flux = {1:.2f} +/- {2:.2f} mJy\n{0} RMS ~ {3:.2f} mJy\n{0} SNR = {4:.2f}".format(image_type, 
                                row["catalog_flux_to_use"]*1.e3,row["catalog_flux_to_use_err"]*1.e3,image_rms*1.e3, row["catalog_flux_to_use"]/image_rms), transform=plt.gcf().transFigure, bbox=bbox_dict)
                        # panels[0].show_rectangles([recentre_ra], [recentre_dec], x_size, y_size, color='C1', linewidth=1, label="Extracted Flux", layer="Extracted Flux")
                
                if skip:
                    continue
                # except:
                #     logger.error("Can't zoom to region for source. Will skip.")
                #     logger.debug("SUMSS Image: {}".format(s_image))
                #     logger.debug("RA: {}".format(recentre_ra))
                #     logger.debug("Dec: {}".format(recentre_dec))
                #     continue
                
                # if not askap_only:
                #     if row["master_name"] in good_sources:
                #
                #         if row["pipelinetag"]=="Match (non-convolved)":
                #
                #         else:
                #             askap_snr_text=plt.text(pos_x_1, pos_y_2, "ASKAP Flux = {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["askap_int_flux"]*1.e3,
                #                 row["askap_rms"]*1.e3, row["measured_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict)
                #         if convolve:
                #             if row["pipelinetag"]=="Match (non-convolved)":
                #
                #             else:
                #                 askap_snr_text_preconv=plt.text(pos_x_5, pos_y_2, "ASKAP Flux = {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["askap_non_conv_int_flux"]*1.e3,
                #                     row["askap_rms_preconv"]*1.e3, row["measured_preconv_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict)
                #     else:
                #         askap_snr_text=plt.text(pos_x_1, pos_y_2, "ASKAP Measured Peak Flux = {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["measured_askap_peak_flux"]*1.e3,
                #             row["askap_rms"]*1.e3, row["measured_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict)
                #         if convolve:
                #             askap_snr_text_preconv=plt.text(pos_x_5, pos_y_2, "ASKAP Measured Peak Flux = {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["measured_preconv_askap_peak_flux"]*1.e3,
                #                  row["askap_rms_preconv"]*1.e3, row["measured_preconv_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict)
                #     if image_type=="SUMSS":
                #         snr_col = "sumss_sumss_snr"
                #     else:
                #         snr_col = "nvss_nvss_snr"
                #     logger.debug("a {} b {} c {}".format(row["master_catflux"],image_rms*1.e3, row[snr_col]))
                #
                # else:
                #     askap_snr_text=plt.text(pos_x_1, pos_y_2, "ASKAP Flux = {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["int_flux"]*1.e3,
                #         row["rms"]*1.e3, row["measured_askap_local_rms"]*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict)
                #     if convolve:
                #         askap_snr_text_preconv=plt.text(pos_x_5, pos_y_2, "ASKAP Flux = {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row["non_conv_int_flux"]*1.e3,
                #         row["rms_preconv"]*1.e3,row["measured_preconv_askap_local_rms"]*1.e3, image_type, row["{}_snr".format(image_type.lower())]), transform=plt.gcf().transFigure, bbox=bbox_dict)
                #     # if image_type == "SUMSS":
                #     #     cat_rms = _sumss_rms(row["dec"])
                #     # else:
                #     #     cat_rms = 0.0005
                #     sumss_snr_text=plt.text(pos_x_2, pos_y_5, "Measured {0} Peak Flux = {1:.2f} mJy\n{0} RMS ~ {2:.2f} mJy".format(image_type, row["measured_catalog_peak_flux"]*1.e3, image_rms*1.e3), transform=plt.gcf().transFigure, bbox=bbox_dict)
                
                #Figure name
                # plt.title(row["sumss_name"])
                figname = "source_{}_postagestamps.jpg".format(row["master_name"])
                if convolve:
                    panels[2]._ax1.legend(custom_lines,  labels)
                else:
                    panels[1]._ax1.legend(custom_lines, labels)
                    # custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                    #                 Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728')]
                    # if dualmode:
                    #     labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "{} Sources".format(other[basecat].upper()), "ASKAP Source Position"]
                    #     custom_lines+=[Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#F5B406'),]
                    # else:
                    #     labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "ASKAP Source Position"]
                    # custom_lines+=[Line2D([0], [0], color='C1'),]
                    # custom_lines.append(Line2D([0], [0], color="w", markeredgecolor='C1', marker="s", fillstyle="none", lw=2))
                    # labels.append("Measured Peak Flux Area")
                    #
                    # if convolve:
                    #     panels[2]._ax1.legend(custom_lines, labels)
                    # else:
                    #     panels[1]._ax1.legend(custom_lines, labels)
                    #
                    # if tag=="BAD":
                    #     askap_source_line=Line2D([0], [0], color='C1', linestyle=linestyle)
                    # else:
                    #     askap_source_line=Line2D([0], [0], color='C1')
                    # custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                    #                 Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),]
                    # if dualmode:
                    #     custom_lines += [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),]
                    #
                    # if row["master_name"] in good_sources:
                    #     custom_lines+=[askap_source_line]
                    # custom_lines+=[Line2D([0], [0], color='C9')]
                    #
                    # if convolve:
                    #     custom_lines.append(Line2D([0], [0], color="w", markeredgecolor='#11BA1C', marker="o", fillstyle="none", lw=2))
                    #
                    # labels=["ASKAP Sources", "{} Sources".format(basecat.upper())]
                    # if dualmode:
                    #     labels+=["{} Source".format(other[basecat].upper())]
                    #
                    # if row["master_name"] in good_sources:
                    #     labels+=["Matched ASKAP"]
                    # labels+=["Matched {}".format(image_type)]
                    #
                    # if convolve:
                    #     labels+=["Non-conv ASKAP Sources"]
                    #
                    # if row["master_name"] not in good_sources:
                    #     custom_lines.append(Line2D([0], [0], color="w", markeredgecolor='C1', marker="s", fillstyle="none", lw=2))
                    #     labels.append("Measured Peak Flux Area")
                        


                plt.savefig(figname, bbox_inches="tight")

                logger.info("Saved figure {}.".format(figname))
                
                for i in panels:
                    if row["type"] == "goodmatch" or row["type"]=="nocatalogmatch":
                        panels[i].remove_layer("ASKAP Source")
                    if row["type"] == "goodmatch" or row["type"]=="noaskapmatch":
                        panels[i].remove_layer("Cat Source")
                        # panels[1].remove_layer("SUMSS Source")
                if row["type"]=="nocatalogmatch":
                    panels[0].remove_layer("Extracted Flux")
                if row["type"]=="noaskapmatch":
                    panels[1].remove_layer("Extracted Flux")
                    if convolve:
                        panels[2].remove_layer("Extracted Flux_2")
                # print "Delete texts"
                # print fig.texts
                # for txt in fig.texts:
                #     txt.set_visible(False)
                askap_snr_text.set_visible(False)
                if convolve:
                    askap_snr_text_preconv.set_visible(False)
                sumss_snr_text.set_visible(False)
                if row["type"] == "goodmatch":
                    sep_text.set_visible(False)
                ratio_text.set_visible(False)
                # debug_num+=1
                # if debug_num==1:
                #     debug_num=0
                #     break
            
            #I think this will clear the SUMSS one, must be a more specific way.    
            plt.gca().remove()

        plt.close()
    
