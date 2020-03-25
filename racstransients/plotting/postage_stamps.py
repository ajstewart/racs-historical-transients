#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import logging
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.lines import Line2D
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D
from astropy import units as u
import multiprocessing as mp
from functools import partial
import os
import sys
import pkg_resources
import subprocess
import pandas as pd
from matplotlib import patches
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import Angle
from astropy.visualization import ZScaleInterval,ImageNormalize, LinearStretch, PercentileInterval
from astropy.wcs.utils import proj_plane_pixel_scales
from racstransients.tools.fitsimage import Askapimage
from astropy.visualization.wcsaxes import SphericalCircle

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


def filter_selavy_components(catalog, catalog_coords, src_coord, imsize):
    #Filter out selavy components outside field of image
    seps = src_coord.separation(catalog_coords)
    mask = seps <= imsize/1.4 #I think cutout2d angle means the width of the image, not a radius hence /2
    #drop the ones we don't need
    return catalog[mask].reset_index(drop=True)

def create_ellipses(df, wcs, color):
    #build ellipse collection for ASKAP image that do not change
    pix_scale = proj_plane_pixel_scales(wcs)
    sx = pix_scale[0]
    sy = pix_scale[1]
    degrees_per_pixel = np.sqrt(sx * sy)
    logger.debug(degrees_per_pixel)

    if "_RAJ2000" in df.columns:
        ra = "_RAJ2000"
        dec = "_DEJ2000"
        a = "MajAxis"
        b = "MinAxis"
        pa = "PA"
    else:
        ra = "ra"
        dec = "dec"
        a = "a"
        b = "b"
        pa = "pa"

    ww = df[a].astype(float)/3600.
    hh = df[b].astype(float)/3600.
    ww /= degrees_per_pixel
    hh /= degrees_per_pixel
    aa = df[pa].astype(float)
    x = df[ra].astype(float)
    y = df[dec].astype(float)
    coordinates = np.column_stack((x, y))
    coordinates = wcs.wcs_world2pix(coordinates, 0)

    patches = [Ellipse(coordinates[i], hh[i], ww[i], aa[i]) for i in range(len(x))]

    collection = PatchCollection(patches, facecolor="None", edgecolor=color, linewidth=1.5)

    return collection

def crossmatch_stamps(crossmatch, askap_data, askap_wcs, selection, nprocs, radius=13./60., contrast=0.2,
            convolve=False, askap_nonconv_img=None, askap_pre_convolve_catalog=None, dualmode=False,
            basecat="sumss"):

    #Here we just need to plot from the transients master df based on the ratio - good and bad no longer exists

    df_to_plot = crossmatch.transients_master_df

    cat_extractions=crossmatch.base_catalog.df
    if dualmode:
        cat2_extractions=crossmatch.merge_catalog.df
    else:
        cat2_extractions = "None"
    askap_extractions=crossmatch.comp_catalog.df

    #Get the unique mosaic files
    fits_mosaics = df_to_plot["master_catalog_Mosaic_path"].unique()

    if convolve:
        if askap_pre_convolve_catalog is not None:
            askap_pre_convolve_catalog = askap_pre_convolve_catalog.df

        if askap_nonconv_img is not None:
            askap_nonconv_data = askap_nonconv_img.data
            askap_nonconv_wcs = askap_nonconv_img.wcs

    else:
        askap_nonconv_data = None
        askap_nonconv_wcs = None

    for mos in fits_mosaics:
        #get the crossmatches
        sub_to_plot = df_to_plot[df_to_plot.master_catalog_Mosaic_path == mos]

        looping_dicts = [row.to_dict() for i,row in sub_to_plot.iterrows()]

        #open the mosaic
        mosimg = Askapimage(mos)
        mosimg.load_wcs()
        mosimg.load_fits_data()
        mosimg.load_position_dimensions()

        produce_multi = partial(produce_postage_stamps_new, askap_data=askap_data, askap_wcs=askap_wcs, mos_data=mosimg.data, mos_wcs=mosimg.wcs,
                            askap_extractions=askap_extractions, cat_extractions=cat_extractions, askap_pre_convolve_catalog=askap_pre_convolve_catalog,
                            radius=radius, contrast=contrast, convolve=convolve, askap_nonconv_data=askap_nonconv_data, askap_nonconv_wcs=askap_nonconv_wcs, basecat=basecat)

        workers=mp.Pool(processes=nprocs)
        try:
            workers.map(produce_multi, looping_dicts)
        except KeyboardInterrupt:
            logger.warning("Caught KeyboardInterrupt, terminating jobs...")
            workers.terminate()
            logger.info("Exiting...")
            workers.close()
            sys.exit()

        #close workers
        workers.close()


def produce_postage_stamps_new(row_dict, askap_data, askap_wcs, mos_data, mos_wcs, askap_extractions, cat_extractions, askap_pre_convolve_catalog, radius=13./60., contrast=0.2,
                            convolve=False, askap_nonconv_data=None, askap_nonconv_wcs=None, basecat="sumss"):
    #For now support one ASKAP image at a time so can just take this from the first entry.

    #First initialise the figure and set up the askap panel with the image.
    radius = Angle(radius * u.arcmin)
    if convolve:
        total_panels = 3
    else:
        total_panels = 2
    panels={}
    other={"sumss":"nvss", "nvss":"sumss"}
    fig = plt.figure(figsize=(18, 8))
    thissurvey = row_dict["survey_used"].upper()

    target = SkyCoord(row_dict["master_ra"]*u.degree, row_dict["master_dec"]*u.degree)
    askap_cutout = Cutout2D(askap_data, position=target, size=radius, wcs=askap_wcs, mode='partial')
    mos_cutout = Cutout2D(mos_data, position=target, size=radius, wcs=mos_wcs, mode='partial')
    if convolve:
        askap_nonconv_cutout = Cutout2D(askap_nonconv_data, position=target, size=radius, wcs=askap_nonconv_wcs, mode='partial')

    askap_norm = ImageNormalize(askap_cutout.data, interval=ZScaleInterval(contrast=contrast))
    mos_norm = ImageNormalize(mos_cutout.data, interval=ZScaleInterval(contrast=contrast))
    if convolve:
        askap_nonconv_norm = ImageNormalize(askap_nonconv_cutout.data, interval=ZScaleInterval(contrast=contrast))

    panels[1] = fig.add_subplot(1,total_panels,1, projection=mos_cutout.wcs)
    panels[2] = fig.add_subplot(1,total_panels,2, projection=askap_cutout.wcs)
    if convolve:
        panels[3] = fig.add_subplot(1,total_panels,3, projection=askap_nonconv_cutout.wcs)

    panels[1].imshow(mos_cutout.data,norm=mos_norm,cmap='gray_r')
    panels[2].imshow(askap_cutout.data,norm=askap_norm,cmap='gray_r')
    if convolve:
        panels[3].imshow(askap_nonconv_cutout.data,norm=askap_nonconv_norm,cmap='gray_r')

    asp = np.diff(panels[1].get_xlim())[0] / np.diff(panels[1].get_ylim())[0]

    asp /= np.abs(np.diff(panels[2].get_xlim())[0] / np.diff(panels[2].get_ylim())[0])

    panels[1].set_aspect(asp)

    # if convolve:
        # panels[3].set_aspect(asp)

    for i in panels:
        panels[i].set_autoscale_on(False)

    askap_collection = create_ellipses(askap_extractions, askap_cutout.wcs, "#1f77b4")

    if convolve:
        askap_nonconv_collection = create_ellipses(askap_pre_convolve_catalog, askap_cutout.wcs, "#11BA1C")
        askap_nonconv_collection_2 = create_ellipses(askap_pre_convolve_catalog, askap_cutout.wcs, "#11BA1C")
        askap_collection_2 = create_ellipses(askap_extractions, askap_cutout.wcs, "#1f77b4")
    else:
        askap_nonconv_collection = None
        askap_nonconv_collection_2 = None
        askap_collection_2 = None

    askap_catalog_collection = create_ellipses(cat_extractions, askap_cutout.wcs, "#d62728")

    catalog_askap_collection = create_ellipses(askap_extractions, mos_cutout.wcs, "#1f77b4")

    if convolve:
        catalog_askap_nonconv_collection = create_ellipses(askap_pre_convolve_catalog, mos_cutout.wcs, "#11BA1C")
        askap_catalog_collection_2 = create_ellipses(cat_extractions, askap_cutout.wcs, "#d62728")
    else:
        catalog_askap_nonconv_collection = None
        catalog_askap_nonconv_collection_2 = None
        askap_catalog_collection_2 = None

    catalog_collection = create_ellipses(cat_extractions, mos_cutout.wcs, "#d62728")

    #Ellipes loading
    panels[1].add_collection(catalog_collection, autolim=False)
    panels[1].add_collection(catalog_askap_collection, autolim=False)

    panels[2].add_collection(askap_collection, autolim=False)
    panels[2].add_collection(askap_catalog_collection, autolim=False)

    if convolve:
        panels[1].add_collection(catalog_askap_nonconv_collection, autolim=False)
        panels[2].add_collection(askap_nonconv_collection, autolim=False)
        panels[3].add_collection(askap_catalog_collection_2, autolim=False)
        panels[3].add_collection(askap_collection_2, autolim=False)
        panels[3].add_collection(askap_nonconv_collection_2, autolim=False)


    bbox_dict=dict(boxstyle="round", ec="white", fc="white", alpha=0.8)

    # if not pd.isna(filtered_cross_matches.iloc[0]["master_catalog_Mosaic_rms"]):
    #     image_rms = filtered_cross_matches.iloc[0]["master_catalog_Mosaic_rms"]
    # else:
    #     image_rms = filtered_cross_matches.iloc[0]["catalog_Mosaic_rms"]

    #Now begin the main loop per source
    # debug_num=0
    if row_dict["type"] == "match":
        askap_title = "RACS "+row_dict["askap_name"]
        if np.isnan(row_dict["aegean_convolved_int_flux"]) and convolve:
            askap_title_2 = "RACS "+row_dict["askap_non_conv_name"]
        else:
            askap_title_2 = "RACS" + row_dict["askap_name"]
        sumss_title = "{} {}".format(thissurvey, row_dict["master_name"])
        recentre_ra=row_dict["master_ra"]
        recentre_dec=row_dict["master_dec"]

    elif row_dict["type"] == "noaskapmatch":
        askap_title = "RACS (Convolved)"
        askap_title_2 = "RACS (Non-Convolved)"
        sumss_title = "{} {}".format(thissurvey, row_dict["master_name"])
        recentre_ra=row_dict["master_ra"]
        recentre_dec=row_dict["master_dec"]
    elif row_dict["type"] == "nocatalogmatch":
        askap_title = "RACS "+row_dict["askap_name"]
        if convolve:
            askap_title_2 = "RACS "+row_dict["askap_non_conv_name"]
        sumss_title = "{} ({})".format(thissurvey, row_dict["master_catalog_mosaic"])
        recentre_ra=row_dict["askap_ra"]
        recentre_dec=row_dict["askap_dec"]

    panels[1].set_title(sumss_title)
    panels[2].set_title(askap_title)
    if convolve:
        panels[3].set_title(askap_title_2)

    #Centre each image on the ASKAP coordinates for clarity
    if convolve:
        pos_x_1 = 0.41
        pos_y_1 = 0.68
        pos_x_2 = 0.14
        pos_y_2 = 0.27
        pos_x_3 = 0.07
        pos_y_3 = 0.12
        pos_x_4 = 0.65
        pos_x_5 = 0.68
        pos_y_5 = 0.27
        size=15
        boxsize = 10
    else:
        pos_x_1 = 0.57
        pos_y_1 = 0.75
        pos_x_2 = 0.14
        pos_y_2 = 0.14
        pos_x_3 = 0.02
        pos_y_3 = 0.02
        pos_x_4 = 0.8
        size = 18
        boxsize = 14


    for i in panels:
        lon = panels[i].coords[0]
        lat = panels[i].coords[1]
        if i==1:
            lon.set_axislabel('RA', size=size)
            lat.set_axislabel('Dec.', size=size)
            lon.set_ticklabel(size=12)
            lat.set_ticklabel(size=12)
        else:
            lon.set_ticklabel_visible(False)
            lat.set_ticklabel_visible(False)
            lon.set_axislabel('')
            lat.set_axislabel('')

    image_rms = row_dict["master_catalog_Mosaic_rms"]
    snr_col = "{0}_{0}_snr".format(thissurvey.lower())

    if row_dict["type"] == "match":
        for p in panels:
            r = SphericalCircle((target.ra, target.dec), 1. * u.arcmin, edgecolor='C9', facecolor='none', transform=panels[p].get_transform('fk5'), linewidth=3)
            panels[p].add_patch(r)
            askap_target = SkyCoord(row_dict["askap_ra"]*u.degree, row_dict["askap_dec"]*u.degree)
            r = SphericalCircle((askap_target.ra, askap_target.dec), 1. * u.arcmin, edgecolor='C1', facecolor='none', transform=panels[p].get_transform('fk5'), linewidth=3)
            panels[p].add_patch(r)
            if p==2:
                sep_text=plt.text(pos_x_3, pos_y_3, "Distance Separation = {:.2f} arcsec".format(row_dict["d2d"]), transform=fig.transFigure, size=size)
                ratio_text=plt.text(pos_x_4, pos_y_3, "Scaled Int. Flux Ratio = {:.2f} +/- {:.2f}".format(row_dict["master_ratio"], row_dict["master_ratio_err"]), transform=fig.transFigure, size=size)
                askap_snr_text=plt.text(pos_x_1, pos_y_2, "ASKAP Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row_dict["askap_flux_to_use"]*1.e3,
                    row_dict["askap_flux_to_use_err"]*1.e3, row_dict["askap_rms"]*1.e3, row_dict["measured_askap_local_rms"]*1.e3), transform=fig.transFigure, bbox=bbox_dict, size=boxsize)
                sumss_snr_text=plt.text(pos_x_2, pos_y_2, "{0} Int. Flux = {1:.2f} +/- {2:.2f} mJy\n{0} RMS ~ {3:.2f} mJy\n{0} SNR = {4:.2f}".format(thissurvey,
                    row_dict["catalog_flux_to_use"]*1.e3,row_dict["catalog_flux_to_use_err"]*1.e3,image_rms*1.e3, row_dict[snr_col]), transform=fig.transFigure, bbox=bbox_dict, size=boxsize)
                if convolve:
                    askap_snr_text_preconv=plt.text(pos_x_5, pos_y_2, "ASKAP Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row_dict["askap_flux_to_use_2"]*1.e3,
                        row_dict["askap_flux_to_use_2_err"]*1.e3, row_dict["askap_rms_preconv"]*1.e3, row_dict["measured_preconv_askap_local_rms"]*1.e3), transform=fig.transFigure, bbox=bbox_dict, size=boxsize)
                    custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4', linestyle="none"),
                                    Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728', linestyle="none"),
                                    Line2D([0], [0], color="w", markeredgecolor='#11BA1C', marker="o", fillstyle="none", lw=2, linestyle="none"),
                                    Line2D([0], [0], color="w", markeredgecolor='C9', marker="o", fillstyle="none", lw=2, linestyle="none"),
                                    Line2D([0], [0], color='C1'),
                                    Line2D([0], [0], color='C9')
                                ]
                    labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "Non-conv ASKAP sources", "Extracted Flux", "ASKAP Source Position", thissurvey+" Source Position"]
                else:
                    custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4', linestyle="none"),
                                    Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728', linestyle="none"),
                                    Line2D([0], [0], color="w", markeredgecolor='C9', marker="o", fillstyle="none", lw=2, linestyle="none"),
                                    Line2D([0], [0], color='C1'),
                                    Line2D([0], [0], color='C9')
                                ]
                    labels = ["ASKAP Sources", "{} Sources".format(basecat.upper()), "Extracted Flux", "ASKAP Source Position", thissurvey+" Source Position"]

    elif row_dict["type"] == "noaskapmatch":
        for p in panels:
            r = SphericalCircle((target.ra, target.dec), 1. * u.arcmin, edgecolor='C9', facecolor='none', transform=panels[p].get_transform('fk5'), linewidth=3)
            panels[p].add_patch(r)
            if p==1:
                ratio_text=plt.text(pos_x_4, pos_y_3, "Scaled Int. Flux Ratio = {:.2f} +/- {:.2f}".format(row_dict["master_ratio"], row_dict["master_ratio_err"]), transform=fig.transFigure, size=size)
            if convolve:
                if p==2:
                    r = SphericalCircle((target.ra, target.dec), 22.5 * u.arcsec, edgecolor='C1', facecolor='none', transform=panels[p].get_transform('fk5'), linewidth=1)
                    panels[p].add_patch(r)
                if p==3:
                    r = SphericalCircle((target.ra, target.dec), 7. * u.arcsec, edgecolor='C1', facecolor='none', transform=panels[p].get_transform('fk5'), linewidth=1)
                    panels[p].add_patch(r)
                if p==1:
                    askap_snr_text=plt.text(pos_x_1, pos_y_2, "ASKAP Measured Int. Flux = {:.2f} mJy +/- {:.2f} \nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row_dict["askap_flux_to_use"]*1.e3,
                        row_dict["askap_flux_to_use_err"]*1.e3, row_dict["askap_rms"]*1.e3, row_dict["measured_askap_local_rms"]*1.e3), transform=fig.transFigure, bbox=bbox_dict, size=boxsize)
                    askap_snr_text_preconv=plt.text(pos_x_5, pos_y_2, "ASKAP Measured Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row_dict["askap_flux_to_use_2"]*1.e3,
                        row_dict["askap_flux_to_use_2_err"]*1.e3, row_dict["askap_rms_preconv"]*1.e3, row_dict["measured_preconv_askap_local_rms"]*1.e3), transform=fig.transFigure, bbox=bbox_dict, size=boxsize)
                    custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                                    Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),
                                    Line2D([0], [0], color="w", markeredgecolor='#11BA1C', marker="o", fillstyle="none", lw=2),
                                    Line2D([0], [0], color="w", markeredgecolor='C9', marker="o", fillstyle="none", lw=2),
                                    Line2D([0], [0], color='C1')
                                ]
                    labels = ["ASKAP Sources", "{} Sources".format(thissurvey), "Non-conv ASKAP sources", "Extracted Flux", thissurvey+" Source Position"]
            else:
                if p==2:
                    r = SphericalCircle((target.ra, target.dec), 7. * u.arcsec, edgecolor='C1', facecolor='none', transform=panels[p].get_transform('fk5'), linewidth=1)
                    panels[p].add_patch(r)
                if p==1:
                    askap_snr_text=plt.text(pos_x_1, pos_y_2, "ASKAP Measured Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row_dict["askap_flux_to_use"]*1.e3,
                        row_dict["askap_flux_to_use_err"]*1.e3, row_dict["askap_rms"]*1.e3, row_dict["measured_askap_local_rms"]*1.e3), transform=fig.transFigure, bbox=bbox_dict, size=boxsize)
                    custom_lines = [Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#1f77b4'),
                                    Line2D([0], [0], color="w", marker="o", fillstyle="none", lw=2, markeredgecolor='#d62728'),
                                    Line2D([0], [0], color="w", markeredgecolor='C9', marker="o", fillstyle="none", lw=2),
                                    Line2D([0], [0], color='C1')
                                ]
                    labels = ["ASKAP Sources", "{} Sources".format(thissurvey), "Extracted Flux", thissurvey+" Source Position"]
            if p==1:
                sumss_snr_text=plt.text(pos_x_2, pos_y_2, "{0} Int. Flux = {1:.2f} +/- {2:.2f} mJy\n{0} RMS ~ {3:.2f} mJy\n{0} SNR = {4:.2f}".format(thissurvey,
                    row_dict["catalog_flux_to_use"]*1.e3, row_dict["catalog_flux_to_use_err"]*1.e3,image_rms*1.e3, row_dict[snr_col]), transform=fig.transFigure, bbox=bbox_dict, size=boxsize)
        #Rectangle on extracted flux on ASKAP images
        # x_size,y_size=_get_extracted_rectangle_size(panels[1], recentre_ra, recentre_dec, 10)
        # panels[1].show_rectangles([recentre_ra], [recentre_dec], x_size, y_size, color='C1', linewidth=1, label="Extracted Flux", layer="Extracted Flux")
        # if convolve:
        #     x_size,y_size=_get_extracted_rectangle_size(panels[2], recentre_ra, recentre_dec, 10)
        #     panels[2].show_rectangles([recentre_ra], [recentre_dec], x_size, y_size, color='C1', linewidth=1, label="Extracted Flux", layer="Extracted Flux")

    elif row_dict["type"] == "nocatalogmatch":
        for p in panels:
            r = SphericalCircle((target.ra, target.dec), 1. * u.arcmin, edgecolor='C1', facecolor='none', transform=panels[p].get_transform('fk5'), linewidth=3)
            panels[p].add_patch(r)
            if p==1:
                ratio_text=plt.text(pos_x_4, pos_y_3, "Scaled Int. Flux Ratio = {:.2f} +/- {:.2f}".format(row_dict["master_ratio"], row_dict["master_ratio_err"]), transform=fig.transFigure, size=size)
                #Put a rectangle to show extracted peak flux measurement
                # x_size,y_size=_get_extracted_rectangle_size(panels[0], recentre_ra, recentre_dec, 5)
                r = SphericalCircle((target.ra, target.dec), 22.5 * u.arcsec, edgecolor='C9', facecolor='none', transform=panels[p].get_transform('fk5'), linewidth=1)
                panels[p].add_patch(r)
                askap_snr_text=plt.text(pos_x_1, pos_y_2, "ASKAP Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row_dict["askap_flux_to_use"]*1.e3,
                    row_dict["askap_flux_to_use_err"]*1.e3,row_dict["askap_rms"]*1.e3, row_dict["measured_askap_local_rms"]*1.e3), transform=fig.transFigure, bbox=bbox_dict, size=boxsize)
                if convolve:
                    askap_snr_text_preconv=plt.text(pos_x_5, pos_y_2, "ASKAP Int. Flux = {:.2f} +/- {:.2f} mJy\nASKAP Image RMS ~ {:.2f} mJy\nASKAP Local RMS ~ {:.2f} mJy".format(row_dict["askap_flux_to_use_2"]*1.e3,
                        row_dict["askap_flux_to_use_2_err"]*1.e3,row_dict["askap_rms_preconv"]*1.e3, row_dict["measured_preconv_askap_local_rms"]*1.e3), transform=fig.transFigure, bbox=bbox_dict, size=boxsize)
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
                sumss_snr_text=plt.text(pos_x_2, pos_y_2, "Measured {0} Int. Flux = {1:.2f} +/- {2:.2f} mJy\n{0} RMS ~ {3:.2f} mJy\n{0} SNR = {4:.2f}".format(thissurvey,
                    row_dict["catalog_flux_to_use"]*1.e3,row_dict["catalog_flux_to_use_err"]*1.e3,image_rms*1.e3, row_dict["catalog_flux_to_use"]/image_rms), transform=fig.transFigure, bbox=bbox_dict, size=boxsize)
            # panels[0].show_rectangles([recentre_ra], [recentre_dec], x_size, y_size, color='C1', linewidth=1, label="Extracted Flux", layer="Extracted Flux")
    #         #
    #         # if skip:
    #         #     continue
            #
    figname = "source_{}_postagestamps.jpg".format(row_dict["master_name"])

    if convolve:
        panels[3].legend(custom_lines,  labels, fontsize=10, loc=1)
    else:
        panels[2].legend(custom_lines, labels, fontsize=14, loc=1)



    plt.savefig(figname, bbox_inches="tight")
    fig.clf()
    #plt.show()
    plt.close(fig)
    return
