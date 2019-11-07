#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import logging
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
import matplotlib.axes as maxes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import aplpy
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.visualization import ZScaleInterval,ImageNormalize, LinearStretch, PercentileInterval
import pandas as pd

logger = logging.getLogger(__name__)


def flux_ratio_image_view_astropy(df, fitsimage, title="Flux ratio plot", save=True, base_filename="image", basecat="sumss", ratio_col="askap_sumss_int_flux_ratio"):
    with fits.open(fitsimage) as hdu:
        header = hdu[0].header
           
    #Remove crazy ratios - likely errors.
    mask=[True if i < 5. else False for i in df[ratio_col]]
    wcs = WCS(header, naxis=2)
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111,projection=wcs)
    
    ratio_plot=ax.scatter(df["askap_ra"].values[mask], df["askap_dec"].values[mask], c=df[ratio_col][mask], cmap="Reds", marker="o", transform=ax.get_transform('world'))
    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_axislabel("Right Ascension (J2000)")
    lat.set_axislabel("Declination (J2000)")
    lon.set_major_formatter("hh:mm:ss")
    lat.set_major_formatter("dd:mm:ss")
    cb=plt.colorbar(ratio_plot, ax=ax)
    if "sumss" in ratio_col and "nvss" in ratio_col:
        cb.set_label("ASKAP / SUMSS & NVSS flux ratio")
    elif basecat=="sumss":
        cb.set_label("ASKAP / SUMSS flux ratio")
    else:
        cb.set_label("ASKAP / NVSS flux ratio")
    
    plt.title(title)
    
    filename="{}_flux_ratio_image_view.png".format(base_filename)
    
    if save:
        plt.savefig(filename, bbox_inches="tight")
        logger.info("Figure {} saved.".format(filename))
    plt.close()
    return filename    
    
def position_offset(df, title="Position offset plot", save=True, base_filename="image", bmaj=45., bmin=45., pa=0.0, basecat="sumss", 
        ra_offset_col="askap_sumss_ra_offset", dec_offset_col="askap_sumss_dec_offset", dualmode=False):
    # filter_df=df[df["d2d"]<=maxsep]
    df = df.sort_values(by=["askap_int_flux",], ascending=[True,])
    if dualmode:
        ra_offset_sumss=df[df["survey_used"]=="sumss"][ra_offset_col]*3600.
        dec_offset_sumss=df[df["survey_used"]=="sumss"][dec_offset_col]*3600.
        ra_offset_nvss=df[df["survey_used"]=="nvss"][ra_offset_col]*3600.
        dec_offset_nvss=df[df["survey_used"]=="nvss"][dec_offset_col]*3600.
    ra_offset=df[ra_offset_col]
    dec_offset=df[dec_offset_col]
    med_ra_offset=ra_offset.median()
    med_dec_offset=dec_offset.median()
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    beam=patches.Ellipse((0,0), bmaj, bmin, pa, 
         edgecolor="k", facecolor="gold", fill=False, alpha=0.8)
    # beam.set_linestyle="dashed"
    # beam.set_edgecolor="gold"
    # beam.set_facecolor="gold"
    # beam.fill=False
    ax.add_artist(beam)
    if dualmode:
        ax.scatter(ra_offset_sumss, dec_offset_sumss, label="SUMSS to ASKAP")
        ax.scatter(ra_offset_nvss, dec_offset_nvss, label="NVSS to ASKAP")
    else:
        theplot=ax.scatter(ra_offset, dec_offset, label="{} to ASKAP".format(basecat.upper()), c=np.log10(df["askap_int_flux"]*1.e3), cmap="plasma")
        cb=plt.colorbar(theplot, ax=ax)
        cb.set_label("log ASKAP Int. Flux (mJy)")
    plt.axhline(0, color="k", ls="--")
    plt.axhline(med_dec_offset, color='r', ls="--", label="Median Offset")
    plt.axvline(0, color="k", ls="--")
    plt.axvline(med_ra_offset, color='r', ls="--")
    plt.title(title)
    plt.xlabel("RA ASKAP - CATALOG (arcsec)")
    plt.ylabel("Dec ASKAP - CATALOG (arcsec)")
    plt.xlim([-20,20])
    plt.ylim([-20,20])
    # custom_lines=[Line2D([0], [0], color="r"),]
    # ax.legend(custom_lines, ['Median Offset'])
    ax.legend()
    ax.text(3,-19, "Median RA Offset = {:.2f} arcsec\nMedian Dec Offset = {:.2f} arcsec".format(med_ra_offset, med_dec_offset))
    filename="{}_position_offset_from_{}.png".format(base_filename, basecat)
    if save:
        plt.savefig(filename, bbox_inches="tight")
        logger.info("Figure {} saved.".format(filename))
    plt.close()
    return filename
    
    
def source_counts(askap_df, crossmatch_df, max_sep, sumss_cat=[], nvss_cat=[], title="Source counts plot", save=True, base_filename="image", basecat="sumss", dualmode=False):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    ax.grid(zorder=1)
    if dualmode or basecat=="sumss":
        num_sumss_sources=len(sumss_cat.df.index)
    else:
        num_sumss_sources=0
    if dualmode or basecat=="nvss":
        num_nvss_sources=len(nvss_cat.df.index)
    else:
        num_nvss_sources=0
    total_askap_sources = len(askap_df.index)
    
    #Calculate SUMSS sources should see
    if dualmode or basecat=="sumss":
        askap_in_sumss = askap_df[askap_df["dec"]<=-30.0].reset_index(drop=True)
        # askap_below_50 = askap_df[askap_df["dec"]<=-50.0].reset_index(drop=True)
        askap_seen_in_sumss = askap_in_sumss[askap_in_sumss["sumss_snr"]>=5.0].reset_index(drop=True)
        # askap_below_50 = askap_below_50[askap_below_50["peak_f"]>=5.0].reset_index(drop=True)
        num_askap_seen_in_sumss = len(askap_seen_in_sumss.index)
    else:
        num_askap_seen_in_sumss = 0
        
        
    if dualmode or basecat=="nvss":
        askap_in_nvss = askap_df[askap_df["dec"]>=-30.0].reset_index(drop=True)
        # askap_below_50 = askap_df[askap_df["dec"]<=-50.0].reset_index(drop=True)
        askap_seen_in_nvss = askap_in_nvss[askap_in_nvss["nvss_snr"]>=5.0].reset_index(drop=True)
        num_askap_seen_in_nvss = len(askap_seen_in_nvss.index)
    else:
        num_askap_seen_in_nvss = 0

    # total_askap_expected=num_askap_below_50+num_askap_above_50+num_askap_above_30
    if max_sep!=None:
        if dualmode or basecat=="sumss":
            num_sumss_sources_matched = len(crossmatch_df[(crossmatch_df["d2d"]<=max_sep) & (crossmatch_df["survey_used"]=="sumss")].index)
        else:
            num_sumss_sources_matched = 0
        if dualmode or basecat=="nvss":
            num_nvss_sources_matched = len(crossmatch_df[(crossmatch_df["d2d"]<=max_sep) & (crossmatch_df["survey_used"]=="nvss")].index)
        else:
            num_nvss_sources_matched = 0
    else:
        if dualmode or basecat=="sumss":
            num_sumss_sources_matched = len(crossmatch_df[crossmatch_df["survey_used"]=="sumss"].index)
        else:
            num_sumss_sources_matched = 0
        if dualmode or basecat=="nvss":
            num_nvss_sources_matched = len(crossmatch_df[crossmatch_df["survey_used"]=="nvss"].index)
        else:
            num_nvss_sources_matched = 0
    num_catalog_sources_matched = num_nvss_sources_matched + num_sumss_sources_matched
    logger.info("Number of ASKAP sources present expected to be seen:")
    logger.info("In SUMSS: {}".format(num_askap_seen_in_sumss))
    logger.info("In NVSS: {}".format(num_askap_seen_in_nvss))
    logger.info("Total: {}".format(num_askap_seen_in_nvss+num_askap_seen_in_sumss))
    logger.info("Catalog sources with match <= {} arcsec: {}".format(max_sep, num_catalog_sources_matched))
    ax.bar((1,2,3,4, 5, 6, 7),(total_askap_sources, num_askap_seen_in_sumss, num_sumss_sources, num_sumss_sources_matched, num_askap_seen_in_nvss,num_nvss_sources, num_nvss_sources_matched), zorder=10, color=['C0', 'C9', 'C9', 'C9', 'C1', 'C1', 'C1'])
    plt.xticks([1,2,3,4,5,6,7], ["Total ASKAP", "Expected ASKAP in SUMSS", "SUMSS Total", "SUMSS matched < {}\" (no extended)".format(max_sep),"Expected ASKAP in NVSS", "NVSS Total", "NVSS matched < {}\" (no extended)".format(max_sep)], rotation=45, ha="right")
    ax.set_ylabel("Number of Sources")
    for i,val in enumerate((total_askap_sources, num_askap_seen_in_sumss, num_sumss_sources, num_sumss_sources_matched, num_askap_seen_in_nvss,num_nvss_sources, num_nvss_sources_matched)):
        ax.text((i+1-0.1), val+20, "{}".format(val),zorder=10)
    filename="{}_source_counts.png".format(base_filename)
    if save:
        plt.savefig(filename, bbox_inches="tight")
        logger.info("Figure {} saved.".format(filename))
    plt.close()
    return filename
    
def flux_ratios_askap_flux(df, max_sep, title="Flux ratio", save=True, base_filename="image", basecat="sumss", ratio_col="askap_sumss_int_flux_ratio", dualmode=False):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    # filter_df=df[df["d2d"]<=maxsep].reset_index(drop=True)
    # ratios=filter_df["askap_int_flux"]/(filter_df["sumss_St"]/1.e3)
    mask=[True if i < 5. else False for i in df[ratio_col]]
    to_plot=df[mask]
    # ratios=ratios[mask]
    if dualmode:
        label="SUMSS/NVSS to ASKAP Crossmatch < {}\"".format(max_sep)
        ylabel="ASKAP / SUMSS & NVSS Int. Flux Ratio"
    elif basecat=="sumss":
        label="SUMSS to ASKAP Crossmatch < {}\"".format(max_sep)
        ylabel="ASKAP / SUMSS Int. Flux Ratio"
    else:
        label="NVSS to ASKAP Crossmatch < {}\"".format(max_sep)
        ylabel="ASKAP / NVSS Int. Flux Ratio"
    if dualmode:
        to_plot_sumss=to_plot[to_plot["survey_used"]=="sumss"]
        to_plot_nvss=to_plot[to_plot["survey_used"]=="nvss"]
        ax.scatter(to_plot_sumss["askap_int_flux"]*1000., to_plot_sumss[ratio_col], label="SUMSS to ASKAP Crossmatch < {}\"".format(max_sep))
        ax.scatter(to_plot_nvss["askap_int_flux"]*1000., to_plot_nvss[ratio_col], label="NVSS to ASKAP Crossmatch < {}\"".format(max_sep))
    else:
        ax.scatter(to_plot["askap_int_flux"]*1000., to_plot[ratio_col], label="{} to ASKAP Crossmatch < {}\"".format(basecat.upper(), max_sep))
    median_flux_ratio=to_plot[ratio_col].median()
    std_flux_ratio=to_plot[ratio_col].std()
    ax.set_xlabel("log ASKAP Flux (mJy)")
    ax.set_ylabel(ylabel)
    plt.title(title)
    ax.axhline(1., color="k", ls="--")
    ax.axhline(median_flux_ratio, label="Median Flux Ratio ({:.2f})".format(median_flux_ratio))
    # ax.axhline(median_flux_ratio+std_flux_ratio, label="Median +/- STD (STD = {:.2f})".format(std_flux_ratio), color="gold")
    # ax.axhline(median_flux_ratio-std_flux_ratio, color="gold")
    ax.fill_between([-1,10000],median_flux_ratio-std_flux_ratio,median_flux_ratio+std_flux_ratio, alpha=0.3, label="Median +/- std", color="#A9E5F4")
    ax.set_xlim([min(to_plot["askap_int_flux"]*1000.)-(min(to_plot["askap_int_flux"]*1000.)*0.1), max(to_plot["askap_int_flux"]*1000.)+(max(to_plot["askap_int_flux"]*1000.)*0.05) ])
    ax.set_xscale("log")
    plt.legend()
    filename="{}_flux_ratio_vs_askap_flux.png".format(base_filename)
    if save:
        plt.savefig(filename, bbox_inches="tight")
        logger.info("Figure {} saved.".format(filename))
    plt.close()
    return filename
    
    
def flux_ratios_distance_from_centre(df, max_sep, title="Flux ratio", save=True, base_filename="image", basecat="sumss", ratio_col="askap_sumss_int_flux_ratio", dualmode=False):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    df = df.sort_values(by=["askap_int_flux",], ascending=[True,])
    # filter_df=df[df["d2d"]<=maxsep].reset_index(drop=True)
    # ratios=filter_df["askap_int_flux"]/(filter_df["sumss_St"]/1.e3)
    mask=[True if i <= 5. else False for i in df[ratio_col]]
    to_plot=df[mask]
    # ratios=ratios[mask]
    if dualmode:
        to_plot_sumss=to_plot[to_plot["survey_used"]=="sumss"]
        to_plot_nvss=to_plot[to_plot["survey_used"]=="nvss"]
        ax.scatter(to_plot_sumss["askap_distance_from_centre"], to_plot_sumss[ratio_col], label="SUMSS to ASKAP Crossmatch < {}\"".format(max_sep))
        ax.scatter(to_plot_nvss["askap_distance_from_centre"], to_plot_nvss[ratio_col], label="NVSS to ASKAP Crossmatch < {}\"".format(max_sep))
    else:
        theplot = ax.scatter(to_plot["askap_distance_from_centre"], to_plot[ratio_col], label="{} to ASKAP Crossmatch < {}\"".format(basecat.upper(), max_sep),c=np.log10(to_plot["askap_int_flux"]*1.e3), cmap="plasma")
        cb=plt.colorbar(theplot, ax=ax)
        cb.set_label("log ASKAP Int. Flux (mJy)")
    median_flux_ratio=to_plot[ratio_col].median()
    std_flux_ratio=to_plot[ratio_col].std()
    ax.set_xlabel("ASKAP Position From Image Centre (deg)")
    if "sumss" in ratio_col and "nvss" in ratio_col:
        ax.set_ylabel("ASKAP / SUMSS & NVSS Int. Flux Ratio")
    elif basecat=="sumss":
        ax.set_ylabel("ASKAP / SUMSS Int. Flux Ratio")
    else:
        ax.set_ylabel("ASKAP / NVSS Int. Flux Ratio")    
    plt.title(title)
    ax.axhline(median_flux_ratio, label="Median Flux Ratio ({:.2f})".format(median_flux_ratio))
    ax.axhline(1., color="k", ls="--")
    # ax.axhline(median_flux_ratio+std_flux_ratio, label="Median +/- STD (STD = {:.2f})".format(std_flux_ratio), color="gold")
    # ax.axhline(median_flux_ratio-std_flux_ratio, color="gold")
    ax.fill_between([-1,6],median_flux_ratio-std_flux_ratio,median_flux_ratio+std_flux_ratio, alpha=0.2, label="Median +/- std", color="#A9E5F4")
    ax.set_xlim([-0.2, 5.2])
    plt.legend()
    filename="{}_flux_ratio_vs_distance_from_centre.png".format(base_filename)
    if save:
        plt.savefig(filename, bbox_inches="tight")
        logger.info("Figure {} saved.".format(filename))
    plt.close()
    return filename
    
  
def image_sources_overlay(image_data, image_wcs, imagename, overlay_cat, overlay_cat_label="sources", overlay_cat_2=pd.DataFrame(), overlay_cat_label_2="None", sumss=False, nvss=False):
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111,projection=image_wcs)
    
    palette = matplotlib.cm.get_cmap('gray_r')
    palette.set_bad('w', 1.0)
    palette.set_over('k', 1.0)
    palette.set_under('w', 1.0)
    
    # masked_array = np.ma.masked_where(np.isnan(image_data),image_data)
    img_norms = ImageNormalize(image_data, interval=PercentileInterval(99.9), stretch=LinearStretch(), clip=False)

        # panels[key].show_contour(images[n-1], colors='red', levels=[3.*12e-6, 4.*12.e-6, 8.*12.e-6, 16*12.e-6])
    im = ax.imshow(image_data, norm=img_norms, cmap=palette)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05, axes_class=maxes.Axes)
    
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Jy')

    ax.set_autoscale_on(False)
    
    if not sumss and not nvss:
        ra = "ra"
        dec = "dec"
        a = "a"
        b = "b"
        pa = "pa"
        plotname=imagename.replace(".fits", "_askap_source_overlay.png")
    else:
        ra = "_RAJ2000"
        dec = "_DEJ2000"
        a = "MajAxis"
        b = "MinAxis"
        pa = "PA"
        if sumss:
            plotname=imagename.replace(".fits", "_sumss_source_overlay.png")
        else:
            plotname=imagename.replace(".fits", "_nvss_source_overlay.png")
    ww = overlay_cat[a].astype(float)/3600.
    hh = overlay_cat[b].astype(float)/3600.
    aa = overlay_cat[pa].astype(float)
    x = overlay_cat[ra].astype(float)
    y = overlay_cat[dec].astype(float)
    if not sumss and not nvss:
        patches = [Ellipse((x[i], y[i]), ww[i]*1.1, hh[i]*1.1, 90.+(180.-aa[i])) for i in range(len(x))]
    else:
        patches = [Ellipse((x[i], y[i]), ww[i]*1.1, hh[i]*1.1, aa[i]) for i in range(len(x))]
    collection = PatchCollection(patches, facecolor="None", edgecolor="#1f77b4", linestyle="--", linewidth=2, transform=ax.get_transform('world'))
    ax.add_collection(collection, autolim=False)
    custom_lines = [Line2D([0], [0], color='#1f77b4'),]
    labels=[overlay_cat_label]
    if not overlay_cat_2.empty:
        ww = overlay_cat_2[a].astype(float)/3600.
        hh = overlay_cat_2[b].astype(float)/3600.
        aa = overlay_cat_2[pa].astype(float)
        x = overlay_cat_2[ra].astype(float)
        y = overlay_cat_2[dec].astype(float)
        if not sumss and not nvss:
            patches = [Ellipse((x[i], y[i]), ww[i]*1.1, hh[i]*1.1, 90.+(180.-aa[i])) for i in range(len(x))]
        else:
            patches = [Ellipse((x[i], y[i]), ww[i]*1.1, hh[i]*1.1, aa[i]) for i in range(len(x))]
        collection = PatchCollection(patches, facecolor="None", edgecolor='#d62728', linestyle="--", linewidth=2, transform=ax.get_transform('world'))
        ax.add_collection(collection, autolim=False)
        custom_lines+=[Line2D([0], [0], color='#d62728'),]
        labels+=[overlay_cat_label_2]
    
    ax.legend(custom_lines, labels)
    
    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_axislabel("Right Ascension (J2000)")
    lat.set_axislabel("Declination (J2000)")
    lon.set_major_formatter("hh:mm:ss")
    lat.set_major_formatter("dd:mm:ss")
    
    fig.savefig(plotname, bbox_inches="tight", dpi=300)
    
    plt.close()
    
    return plotname
    
            

