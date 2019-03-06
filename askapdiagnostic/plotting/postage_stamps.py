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

def _plotinitial(panels, key, figure, image, vmin=-99, vmax=-99):
    panels[key]=aplpy.FITSFigure(image, figure=figure, subplot=(1,2,key+1))
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
        thestamp="SUMSS_{}_{}_sidebyside.jpg".format(i, qual_tag)
        new_image_name="transient_{}_SUMSS_{}_{}_sidebyside.jpg".format(t_type, i, qual_tag)
        cmd=base_cmd+"{} {}".format(thestamp, new_image_name)
        try:
            subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            logger.error("Copying of {} for transients failed.".format(thestamp))


def create_askap_postage_stamps(askap_df,crossmatch_df,askap_image, nprocs, sumss_mosaic_dir, sumss_extractions, askap_extractions, radius=13./60.):
    try:
        logger.info("Loading SUMSS image data.")
        sumss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/sumss_images_info.csv"))
    except:
        logger.error("SUMSS mosaic data cannot be found!")
        return
    # load up the SUMSS image centre coords:
    sumss_centres = SkyCoord(ra=sumss_mosaic_data["center-ra"].values*u.deg, dec=sumss_mosaic_data["center-dec"].values*u.deg)
    
    matching_images=[]
    sumss_fluxes=[]
    
    for i,row in askap_df.iterrows():
        askap_target=SkyCoord(ra=row["ra"]*u.deg, dec=row["dec"]*u.deg)
        seps = askap_target.separation(sumss_centres)
        min_index = np.argmin(seps.deg)
        # print d2d_sumss.deg[min_index]
        # print min_index
        image = sumss_mosaic_data.iloc[min_index]["image"]
        # print image
        # print askap_target.to_string('hmsdms')
        matching_images.append(image.replace(".FITS", ""))
        #Open the image using aplpy as it has simple world2pixel function
        sumssfits = aplpy.FITSFigure(os.path.join(sumss_mosaic_dir, image))
        source_pixel_loc = sumssfits.world2pixel(askap_target.ra.deg, askap_target.dec.deg)
        source_pixel_loc_x = int(source_pixel_loc[0])
        source_pixel_loc_y = int(source_pixel_loc[1])
        temp_flux_values = []
        for i in [[0,0], [0, 1], [0,-1], [1,1], [1,0], [1,-1], [-1,-1], [-1, 0], [-1,1]]:
            temp_flux_values.append(sumssfits._data[source_pixel_loc_y+i[1], source_pixel_loc_x+i[0]])
        peak_flux = np.max(temp_flux_values)
        sumss_fluxes.append(peak_flux)
        sumssfits.close()
        
    askap_df["sumss_Mosaic"]=matching_images
    askap_df["sumss_peak_flux"]=sumss_fluxes
    
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
        produce_postage_stamps(askap_df, crossmatch_df, askap_image, sumss_mosaic_dir, [],[], sumss_extractions, askap_extractions, radius=13./60., max_separation=None, askap_only=True)
    

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

def crossmatch_stamps(crossmatch, askap_image, postage_options, nprocs,sumss_mosaic_dir, radius=13./60., max_separation=15.):
    # For reference, the transient df's are:
    # self.transients_no_matches_df=no_matches - crossmatch format
    # self.transients_large_ratios_df=large_ratios - crossmatch format
    # self.transients_not_matched_askap_should_see_df=not_matched_askap_sources_should_see - askap_cat -format
    
    #work out good and bad
    good_df=crossmatch.crossmatch_df[crossmatch.crossmatch_df["d2d"]<=max_separation].reset_index(drop=True)
    bad_df=crossmatch.crossmatch_df[crossmatch.crossmatch_df["d2d"]>max_separation].reset_index(drop=True)
    
    good_sources=good_df["sumss_name"].values
    bad_sources=bad_df["sumss_name"].values
    
    if "all" in postage_options:
        mode="all"
        df = crossmatch.crossmatch_df
    elif "good" in postage_options:
        mode="good"
        df = good_df
    else:
        mode="bad"
        df = bad_df

    sumss_extractions=crossmatch.base_catalog.df
    askap_extractions=crossmatch.comp_catalog.df

    if nprocs >1:
        #Get the unique SUMSS images.
        sumss_fits_mosaics = df["sumss_Mosaic"].unique()
        logger.info("{} unique SUMSS images needed for postage stamps".format(len(sumss_fits_mosaics)))
    
        sumss_looping_dict=_get_stamp_looping_parameters(df, sumss_fits_mosaics)
        #below is broken
        postage_stamps_multi=partial(produce_comp_postage_stamps_multicore, askap_image=askap_image, params_dict=sumss_looping_dict, df=df, sumss_mosaic_dir=sumss_mosaic_dir, good_sources=good_sources, 
            bad_sources=bad_sources, radius=radius)
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
            sumss_mosaic_dir,good_sources,bad_sources, sumss_extractions, askap_extractions, radius=13./60., max_separation=None, askap_only=False)
        
    if "transients" in postage_options:
        if mode=="all":
            _copy_images_to_transient_folders(crossmatch.transients_no_matches_df["sumss_name"].tolist(), "NOMATCH", good_sources, bad_sources)
            _copy_images_to_transient_folders(crossmatch.transients_large_ratios_df["sumss_name"].tolist(), "LARGERATIO", good_sources, bad_sources)
        else:
            sumss_done=df["sumss_name"].tolist()
            df_required=crossmatch.crossmatch_df[~crossmatch.crossmatch_df.sumss_name.isin(sumss_done)].reset_index(drop=True)
            if nprocs >1:
                required_sumss_fits_mosaics = df_required["sumss_Mosaic"].unique()
                required_sumss_looping_dict=_get_stamp_looping_parameters(df_required, required_sumss_fits_mosaics)
                postage_stamps_multi=partial(produce_comp_postage_stamps_multicore, askap_image=askap_image, params_dict=required_sumss_looping_dict, df=df_required, sumss_mosaic_dir=sumss_mosaic_dir, radius=radius)
                try:
                    workers.map(postage_stamps_multi, sumss_fits_mosaics)
                except KeyboardInterrupt:
                    logger.warning("Caught KeyboardInterrupt, terminating jobs...")
                    workers.terminate()
                    logger.info("Exiting...")
                    workers.close()
                    sys.exit()
                #ASKAP only ones
                workers.close()
            else:
                produce_postage_stamps(df_required, crossmatch.crossmatch_df, askap_image, sumss_mosaic_dir, good_sources, bad_sources, radius=13./60., max_separation=None, askap_only=False) 
        
        create_askap_postage_stamps(crossmatch.transients_not_matched_askap_should_see_df, crossmatch.crossmatch_df, askap_image,nprocs, sumss_mosaic_dir, sumss_extractions, askap_extractions)
        
            
            
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
    
    
def produce_postage_stamps(df, full_df, askap_fits, sumss_mosaic_dir, good_sources, bad_sources, sumss_extractions, askap_extractions, radius=13./60., max_separation=None, askap_only=False):
        logger.info("Estimated time to completion = {:.2f} hours".format(len(df.index)*6./3600.))
        #Minimise fits opening so first get a list of all the unique SUMSS fits files to be used
        sumss_fits_mosaics = df["sumss_Mosaic"].unique()
        #For now support one ASKAP image at a time so can just take this from the first entry.
        
        #First initialise the figure and set up the askap panel with the image.
        panels={}
        fig = plt.figure(figsize=(12, 6))
        panels=_plotinitial(panels, 1, fig, askap_fits, vmin=-0.0025, vmax=0.007829)
        
        #Load the ellipses onto the askap image
        panels[1].show_ellipses(askap_extractions["ra"],askap_extractions["dec"],askap_extractions["b"]/3600.*1.5, 
            askap_extractions["a"]/3600.*1.5, angle=askap_extractions["pa"], layer="ASKAP Sources", color="#1f77b4")
        panels[1].show_ellipses(sumss_extractions["_RAJ2000"],sumss_extractions["_DEJ2000"],sumss_extractions["MinAxis"]/3600., 
            sumss_extractions["MajAxis"]/3600., angle=sumss_extractions["PA"], layer="SUMSS Sources", color="#d62728")
        
        bbox_dict=dict(boxstyle="round", ec="white", fc="white", alpha=0.7)
        
        #Now start the SUMSS loop
        for s_image in sumss_fits_mosaics:# [:1]:
            #Get the actual fits file
            logger.debug("SUMSS image: {}".format(s_image))
            s_image_path=os.path.join(sumss_mosaic_dir, s_image+".FITS")
            #Filter the dataframe such that only the SUMSS sources are present
            filtered_cross_matches=df[df["sumss_Mosaic"]==s_image].reset_index(drop=True)
            #Generate the base SUMSS panel
            panels=_plotinitial(panels, 0, fig, s_image_path)
            # panels[1].set_title("SUMSS")
            #Add the sources
            panels[0].show_ellipses(askap_extractions["ra"],askap_extractions["dec"],askap_extractions["b"]/3600., 
                askap_extractions["a"]/3600., angle=askap_extractions["pa"], layer="ASKAP Sources", color="#1f77b4")
            panels[0].show_ellipses(sumss_extractions["_RAJ2000"],sumss_extractions["_DEJ2000"],sumss_extractions["MinAxis"]/3600., 
                sumss_extractions["MajAxis"]/3600., angle=sumss_extractions["PA"], layer="SUMSS Sources", color="#d62728")
            panels[1].axis_labels.hide()
            panels[1].tick_labels.hide()
            #Now begin the main loop per source
            # debug_num=0
            for i, row in filtered_cross_matches.iterrows():
                if not askap_only:
                    panels[1].set_title("ASKAP "+row["askap_name"])
                    panels[0].set_title("SUMSS "+row["sumss_name"])
                    recentre_ra=row["sumss__RAJ2000"]
                    recentre_dec=row["sumss__DEJ2000"]
                else:
                    panels[1].set_title("ASKAP "+row["name"])
                    panels[0].set_title("SUMSS")
                    recentre_ra=row["ra"]
                    recentre_dec=row["dec"]
                #Centre each image on the ASKAP coordinates for clarity

                
                try:
                    for p in panels:
                        panels[p].recenter(recentre_ra, recentre_dec, radius)

                        if not askap_only:
                            panels[p].show_circles([row["askap_ra"]], [row["askap_dec"]], 120./3600., color='C1', linewidth=3, label="ASKAP source", layer="ASKAP Source")
                            panels[p].show_circles([recentre_ra], [recentre_dec],120./3600., color='C9', linewidth=3, label="SUMSS source", layer="SUMSS Source")
                        else:
                            panels[p].show_circles([recentre_ra], [recentre_dec], 120./3600., color='C1', linewidth=3, label="ASKAP position", layer="ASKAP Source")
                except:
                    logger.error("Can't zoom to region for source. Will skip.")
                    continue

                
                if not askap_only:
                    sep_text=plt.text(0.02, 0.02, "Distance Separation = {:.2f} arcsec".format(row["d2d"]), transform=plt.gcf().transFigure)
                    ratio_text=plt.text(0.8, 0.02, "Int. Flux Ratio ASKAP/SUMSS = {:.2f}".format(row["askap_sumss_int_flux_ratio"]), transform=plt.gcf().transFigure)
                    askap_snr_text=plt.text(0.57, 0.75, "ASKAP Flux = {:.2f} mJy\nASKAP Sigma = {:.2f}\nSUMSS Sigma = {:.2f}".format(row["askap_int_flux"]*1.e3, row["askap_snr"], row["askap_sumss_snr"]), transform=plt.gcf().transFigure, bbox=bbox_dict)
                    sumss_snr_text=plt.text(0.14, 0.78, "SUMSS Flux = {:.2f} mJy\nSUMSS Sigma ~ {:.2f}".format(row["sumss_St"], row["sumss_sumss_snr"]), transform=plt.gcf().transFigure, bbox=bbox_dict)
                else:
                    askap_snr_text=plt.text(0.57, 0.75, "ASKAP Flux = {:.2f} mJy\nASKAP Sigma = {:.2f}\nSUMSS Sigma ~ {:.2f}".format(row["int_flux"]*1.e3, row["snr"], row["sumss_snr"]), transform=plt.gcf().transFigure, bbox=bbox_dict)
                    sumss_snr_text=plt.text(0.14, 0.78, "SUMSS Peak Flux = {:.2f} mJy\nSUMSS Sigma ~ {:.2f}".format(row["sumss_peak_flux"]*1.e3, row["sumss_peak_flux"]/_sumss_rms(row["dec"])), transform=plt.gcf().transFigure, bbox=bbox_dict)
                
                #Figure name
                # plt.title(row["sumss_name"])
                if askap_only:
                    figname = "transient_askapnotseen_ASKAP_{}_sidebyside.jpg".format(row["name"])
                    custom_lines = [Line2D([0], [0], color='#1f77b4'),
                                    Line2D([0], [0], color='#d62728'),    
                                    Line2D([0], [0], color='C1')]    
                    panels[1]._ax1.legend(custom_lines, ["ASKAP Sources", "SUMSS Sources", "ASKAP Source Position"])
                else:
                    if row["sumss_name"] in good_sources:
                        tag="GOOD"
                    else:
                        tag="BAD"
                    figname = "SUMSS_{}_{}_sidebyside.jpg".format(row["sumss_name"], tag)
                    custom_lines = [Line2D([0], [0], color='#1f77b4'),
                                    Line2D([0], [0], color='#d62728'),    
                                    Line2D([0], [0], color='C1'),    
                                    Line2D([0], [0], color='C9')]
                    panels[1]._ax1.legend(custom_lines, ["ASKAP Sources", "SUMSS Sources", "Matched ASKAP", "Matched SUMSS"])

                plt.savefig(figname, bbox_inches="tight")

                logger.info("Saved figure {}.".format(figname))
        
                panels[0].remove_layer("ASKAP Source")
                panels[1].remove_layer("ASKAP Source")
                if not askap_only:
                    panels[0].remove_layer("SUMSS Source")
                    panels[1].remove_layer("SUMSS Source")

                askap_snr_text.set_visible(False)
                sumss_snr_text.set_visible(False)
                if not askap_only:
                    sep_text.set_visible(False)
                    ratio_text.set_visible(False)
                # debug_num+=1
                # if debug_num==3:
                #     debug_num=0
                #     break
            
            #I think this will clear the SUMSS one, must be a more specific way.    
            plt.gca().remove()

        plt.close()
    
