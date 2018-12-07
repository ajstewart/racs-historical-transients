#!/usr/bin/env python

import logging
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def flux_ratio_image_view(basecat_fluxes, secondcat_fluxes, max_sep, idx, d2d, coords, save=True):
    flux_ratio=[]
    plot_ra=[]
    plot_dec=[]
    for i,val in enumerate(d2d):
        if val.arcsec <= max_sep:
            baseflux=basecat_fluxes[i]
            compflux=secondcat_fluxes[idx[i]]
            ratio=compflux/baseflux
            plot_ra.append(coords[idx[i]].ra.deg)
            plot_dec.append(coords[idx[i]].dec.deg)
            flux_ratio.append(ratio)

   
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)


    # min_ratio=min(flux_ratio)
    # max_ratio=max(flux_ratio)
    #
    # norm = matplotlib.colors.Normalize(vmin=min_ratio, vmax=max_ratio, clip=True)
    # mapper = cm.ScalarMappable(norm=norm)

    ratio_plot=ax.scatter(plot_ra, plot_dec, c=flux_ratio, cmap="jet")

    plt.colorbar(ratio_plot, ax=ax)

    plt.gca().invert_xaxis()
    if save:
        logger.info("Plot saved")
        plt.savefig("flux_ratio_location.pdf", bbox_inches="tight")


