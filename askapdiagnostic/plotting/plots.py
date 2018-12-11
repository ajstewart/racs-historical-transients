#!/usr/bin/env python

import logging
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.lines import Line2D

logger = logging.getLogger(__name__)


def flux_ratio_image_view(df, maxsep, title="FLux ratio plot", save=True):
    filter_df=df[df["d2d"]<=maxsep].reset_index(drop=True)
    ratios=filter_df["askap_peak_flux"]/(filter_df["sumss_Sp"]/1.e3).values
    #sloppy check
    mask=[True if i < 5. else False for i in ratios]
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    ratio_plot=ax.scatter(filter_df["askap_ra"].values[mask], filter_df["askap_dec"].values[mask], c=ratios[mask], cmap="Blues", marker="o")

    cb=plt.colorbar(ratio_plot, ax=ax)
    cb.set_label("ASKAP / SUMSS flux ratio")
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")
    plt.gca().invert_xaxis()
    
    plt.title(title)
    
    if save:
        logger.info("Plot saved")
        plt.savefig("flux_ratio_location_pandas.png", bbox_inches="tight")
    plt.close()
    
    
def position_offset(df, maxsep, title="Position offset plot", save=True):
    filter_df=df[df["d2d"]<=maxsep]
    ra_offset=(filter_df["askap_ra"]-filter_df["sumss__RAJ2000"])*3600.
    dec_offset=(filter_df["askap_dec"]-filter_df["sumss__DEJ2000"])*3600.
    med_ra_offset=ra_offset.median()
    med_dec_offset=dec_offset.median()
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    beam=patches.Ellipse((0,0), filter_df.iloc[0]["askap_psf_a"],filter_df.iloc[0]["askap_psf_b"], filter_df.iloc[0]["askap_psf_pa"], 
            linestyle="--", edgecolor="gold", facecolor="gold", fill=False)
    # beam.set_linestyle="dashed"
    # beam.set_edgecolor="gold"
    # beam.set_facecolor="gold"
    # beam.fill=False
    ax.add_artist(beam)
    ax.scatter(ra_offset, dec_offset)
    plt.axhline(0, color="k", ls="--")
    plt.axhline(med_dec_offset, color='r', ls="--")
    plt.axvline(0, color="k", ls="--")
    plt.axvline(med_ra_offset, color='r', ls="--")
    plt.title(title)
    plt.xlabel("RA (arcsec)")
    plt.ylabel("Dec (arcsec)")
    plt.xlim([-20,20])
    plt.ylim([-20,20])
    custom_lines=[Line2D([0], [0], color="r"),]
    ax.legend(custom_lines, ['Median Offset'])
    ax.text(10,-19, "Median RA Offset = {:.2f}\nMedian Dec Offset = {:.2f}".format(med_ra_offset, med_dec_offset))
    if save:
        logger.info("Plot saved")
        plt.savefig("sumss-offset.png", bbox_inches="tight")
    plt.close()
    
    
def source_counts(sumss_df, askap_df, title="Source counts plot", save=True):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    ax.grid(zorder=1)
    num_sumss_sources=len(sumss_df.index)
    total_askap_sources = len(askap_df.index)
    askap_above_50 = askap_df[askap_df["dec"]>-50.0].reset_index(drop=True)
    askap_below_50 = askap_df[askap_df["dec"]<=-50.0].reset_index(drop=True)
    askap_above_50 = askap_above_50[askap_above_50["peak_flux"]>=0.010].reset_index(drop=True)
    askap_below_50 = askap_below_50[askap_below_50["peak_flux"]>=0.006].reset_index(drop=True)
    num_askap_above_50=len(askap_above_50.index)
    num_askap_below_50=len(askap_below_50.index)
    total_askap_expected=num_askap_below_50+num_askap_above_50
    logger.info("Number of ASKAP sources present expected to be seen:")
    logger.info("Dec > -50 deg (>= 10 mJy): {}".format(num_askap_above_50))
    logger.info("Dec <= -50 deg (>= 6 mJy): {}".format(num_askap_below_50))
    logger.info("Total: {}".format(total_askap_expected))
    ax.bar((1,2,3),(total_askap_sources, total_askap_expected, num_sumss_sources), zorder=10)
    plt.xticks([1,2,3], ["Total ASKAP", "Expected ASKAP in SUMSS", "SUMSS Total"])
    ax.set_ylabel("Number of Sources")
    for i,val in enumerate((total_askap_sources, total_askap_expected, num_sumss_sources)):
        ax.text((i+1-0.1), val+20, "{}".format(val),zorder=10)
    if save:
        logger.info("Plot saved")
        plt.savefig("askap-sumss-source-counts.png", bbox_inches="tight")
    plt.close()
    
    
def flux_ratios(df, maxsep, title="Flux ratio", save=True):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    filter_df=df[df["d2d"]<=maxsep].reset_index(drop=True)
    ratios=filter_df["askap_peak_flux"]/(filter_df["sumss_Sp"]/1.e3)
    mask=[True if i < 5. else False for i in ratios]
    ratios=ratios[mask]
    ax.scatter(range(1, len(ratios)+1), ratios)
    median_flux_ratio=ratios.median()
    std_flux_ratio=ratios.std()
    ax.set_xlabel("Source number")
    ax.set_ylabel("ASKAP / SUMSS Flux Ratio")
    plt.title(title)
    ax.axhline(median_flux_ratio, label="Median Flux Ratio ({:.2f})".format(median_flux_ratio), color="red")
    ax.axhline(median_flux_ratio+std_flux_ratio, label="Median +/- STD (STD = {:.2f})".format(std_flux_ratio), color="gold")
    ax.axhline(median_flux_ratio-std_flux_ratio, color="gold")
    plt.legend()
    if save:
        logger.info("Plot saved")
        plt.savefig("askap-sumss-flux-ratio.png", bbox_inches="tight")
    plt.close()
    
        
    
    


