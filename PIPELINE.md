# racs-historical-transients pipeline description

Below is a description of what the pipeline does to search for transients and explanations of the different tags and fluxes that the pipeline provides.

## Pipeline Transient Steps

1. ASKAP image provided is read in and if the weight trim or convolve options are selected, it will produce/open these images as well. If both are selected then the weight cropping is done first, followed by convolving. If convolution is to be done by the pipeline then CASA is used to do so. The target resolution is decided by the declination in SUMSS case or is 45 x 45 arcsec for NVSS.

2. Fetch the SUMSS or NVSS sources for the image area from Vizier (recommended), or use an already obtained list.

3. Read in source catalogues or run aegean on the ASKAP image, both normal and convolved (if selected). All the ASKAP fluxes are scaled to the SUMSS or NVSS frequency assuming a spectral index of 0.8.

4. If convolved is selected, the non-convolved catalogue is crossmatched to the convolved one, so that each convolved component has a non-convolved counterpart.

5. **Transient Crossmatching**: if convolved is selected then all the base comparison is done to the convolved image, with the non-convolved information looked at when needed. The crossmatch is performed between SUMSS or NVSS and the convolved ASKAP catalogue.

6. The crossmatch is analysed and split into three categories:
    - Good matches: there is a match to the SUMSS/NVSS component within the user defined distance limit.
    - No ASKAP match: SUMSS/NVSS sources for which the ASKAP crossmatch is outside the user defined distance limit.
    - No catalog match: An ASKAP source that has not been matched to a catalog (SUMSS/NVSS) source when it should have been seen (i.e. it's flux is high enough to be seen in the previous surveys).
    
7. Local RMS values are measured for each crossmatch source (300 pixel box, using astropy sigma clipping to 4 sigma).
    
8. No ASKAP matches are checked to see if there is actually a match in the non-convolved catalogue. If there is then the information from the non-convolved catalogue is subbed in. This is the only situation where this happens.

9. No ASKAP match and no catalog match sources are forcefully extracted using aegean in the ASKAP and catalog image respectively. These measured values are then used to perform the ratio calculation.

10. All crossmatches are analysed and the pipeline flags those which meet a defined criteria that could signify bad matches.

11. Flux ratios and variability metrics are calculated using the catlaog integrated flux and the scaled ASKAP integrated flux. The `Vs` metric takes into account the errors on the fluxes for which the user can input a percentage error for the ASKAP fluxes.

12. If selected cutout plots are made for all crossmatches.

13. Results injected into database ready for viewing.


## Pipeline Fluxes
When looking at the website there are a range of fluxes that are reported and it can be confusing as to what flux is used where. This attempts to explain:

### Flux Ratio and Variability Metrics
These are always calculated using the **scaled convolved ASKAP integrated flux** and the catalogue integrated flux.

There are two situations where this isn't true:

1. Where the convolution has convolved out an ASKAP source, so the crossmatch uses the non-convolved information directly, specifically the **scaled non-convolved ASKAP integrated flux** If this has happened then the `Using Non-convolved info` flag will be true in the flux table on the crossmatch detail page.

2. The measured Aegean extracted flux in one of the `no match` types is less than 3x the measured local rms, so the respective value is replaced with 3x the local rms value. If this is true then the `3x RMS used` is True.

### Query
Here the fluxes you query are the direct measured integrated fluxes for each.

### Cutout Plots
The plot values show the directly measured fluxes from the catalogues or Aegean. So sometimes this will show negative fluxes for forced extractions and odd rms values as these can also be from Aegean. The RMS image value comes from running sigma clipping on the entire image.

## Pipeline Flags
* **Good:** Considered a good match with no obvious issues.
* **Created Convolved Source:** There is no non-convolved match to the convolved ASKAP source within 2 beam widths. Likely a source that has been created by convolving.
* **Edge of ASKAP Image:** Source is within 24 pixels of an image edge (1 arcmin).
* **Extended Source:** The source is classed as extended if either of the source axes in either the catalog or ASKAP image is larger than 1.75 x the respective image beam axes.
* **Failed:** Source has failed to be measured. Usually caused by Aegean returning NaN values.
* **Edge of ASKAP Image:** Source is within 24 pixels of an image edge (1 arcmin).
* **Inflated convolved flux error:** This means that the convolved ASKAP flux has been artificially inflated such that the source is now above the SUMSS detection limit when it really shouldn't have been seen.
* **Large island source:** In either the convolved or non-convolved island catalogue the componenet is a member of an island with 3 or more components.
* **Likely artefact (bright source):** The crossmatch is less than 10 arcminutes away from a previously defined bright source in the SUMSS or NVSS data. Bright sources are determined by being more than two standarded deviations away from the median of the flux distribution of the entire catalogue.
* **Likely bright extended structure:** The source in question has a signal to noise > 50 in the ASKAP image.
* **Likely double/multi source:** There is a source within 2.5 x 45 arcsec of the crossmatched source in question (convolved catalogue if used).

## Other Flags
* **Flux convolved error:** For convolved sources that have a definite non-convolved counterpart, if the difference between these fluxes is more than 25% then this is flagged as true. It can cause crossmatches to have a higher (or lower) ratio than they actually should. This is checked when the non-convolved crossmatch is less than 2 (non-convolved) major axis beam witdhs away from the convovled source and the nearest non-convolved neighbour is more than 3 (non-convolved) major axis beamwidths away.

## Variability Metrics
The standard two epoch metrics have been used, as defined in [this paper](https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.4024R/abstract). `Vs` takes into account the user entered percentage error on the ASKAP fluxes.

## Query Arguments
Most of the options in the query are self-explantory, however these options perhaps need clarifying:

* **Distance to nearest ASKAP neighbour:** This distance takes the crossmatch coordinate and finds the nearest next ASKAP source from the non-convolved catalogue (i.e. not the source that has been crossmatched).
* **m:** The m metric is queried using the absolute value. Hence querying a minimum of 1.0 will return sources with an `m > 1` and `m < -1`.

## Aegean errors
The forced extraction with aegean sometimes does not return valid values. Namely the the errors on integrated fluxes are sometimes 0 and so to can the local rms be reported as 0. Also the integrated fluxes returned can be negative. In this cases the following solutions are used:

* SUMSS image extraction:
    - Int flux = 0: Use 3x the base sensitvity limit (3 x 1.2 mJy < -50, 3 x 2 mJy > -50).
    - Int flux error = 0: Use the base sensitvity limit (1.2 mJy < -50, 2 mJy > -50).
    - Local rms = 0: Use the base sensitvity limit (1.2 mJy < -50, 2 mJy > -50).
* ASKAP image extraction:
    - Int flux = 0: Use 3x the scaled-to-catalog measured local rms.
    - Int flux error = 0: Use the scaled-to-catalog measured local rms.
    - Local rms = 0: Use the scaled-to-catalog measured local rms.
