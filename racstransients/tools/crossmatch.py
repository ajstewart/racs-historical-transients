#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import logging
import aplpy
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.nddata.utils import Cutout2D
from racstransients.tools.fitsimage import Askapimage
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D
from racstransients.plotting import postage_stamps
import sqlalchemy
import psycopg2
import pkg_resources
import numpy as np
import pandas as pd
import subprocess
import multiprocessing


class crossmatch(object):
    """docstring for crossmatch"""
    def __init__(self, base_catalog, comp_catalog, base_catalogue_name="sumss", logger=None):
        self.logger = logger or logging.getLogger(__name__)
        super(crossmatch, self).__init__()
        self.base_catalog=base_catalog
        self.comp_catalog=comp_catalog
        self.performed=False 
        self.base_catalogue_name=base_catalogue_name
        
        #Save transient values for later
        self.transients_noaskapmatchtocatalog_total=0
        self.transients_noaskapmatchtocatalog_candidates=0
        self.transients_nocatalogmatchtoaskap_total=0
        self.transients_nocatalogmatchtoaskap_candidates=0
        self.transients_largeratio_total=0
        self.transients_goodmatches_total=0
        
    def perform_crossmatch(self):
        self.idx, self.d2d, self.d3d = self.base_catalog._crossmatch_catalog.match_to_catalog_sky(self.comp_catalog._crossmatch_catalog)
        self.matches = self.comp_catalog.df.loc[self.idx].reset_index(drop=True)
        newcolnames={}
        for c in self.matches.columns:
            newcolnames[c]="askap_{}".format(c)
        self.matches=self.matches.rename(columns=newcolnames)
        self.crossmatch_df = self.base_catalog.df.copy(deep=True)
        newcolnames={}
        for c in self.crossmatch_df.columns:
            newcolnames[c]="{}_{}".format(self.base_catalogue_name, c)
        self.crossmatch_df=self.crossmatch_df.rename(columns=newcolnames)
        self.crossmatch_df=self.crossmatch_df.join(self.matches)
        # for i in self.crossmatch_df.columns:
        #     print i
        self.crossmatch_df["d2d"]=self.d2d.arcsec
        self.crossmatch_df["{}_d2d".format(self.base_catalogue_name)]=self.d2d.arcsec
        self.crossmatch_df["survey_used"]=self.base_catalogue_name
        self.crossmatch_df["master_name"]=self.crossmatch_df["{}_name".format(self.base_catalogue_name)]
        self.crossmatch_df["master_ra"]=self.crossmatch_df["{}_{}".format(self.base_catalogue_name, self.base_catalog.ra_col)]
        self.crossmatch_df["master_dec"]=self.crossmatch_df["{}_{}".format(self.base_catalogue_name, self.base_catalog.dec_col)]
        self.crossmatch_df["master_iflux"]=self.crossmatch_df["askap_{}".format(self.comp_catalog.flux_col)]
        self.crossmatch_df["master_catflux"]=self.crossmatch_df["{}_{}".format(self.base_catalogue_name, self.base_catalog.flux_col)]
        self.crossmatch_df["master_scaled_askap_iflux"]=self.crossmatch_df["askap_askap_scaled_to_{}".format(self.base_catalogue_name)]
        if "askap_non_conv_askap_scaled_to_sumss" in self.crossmatch_df.columns or "askap_non_conv_askap_scaled_to_nvss" in self.crossmatch_df.columns:
            self.crossmatch_df["master_scaled_non_conv_askap_iflux"]=self.crossmatch_df["askap_non_conv_askap_scaled_to_{}".format(self.base_catalogue_name)]
            self.crossmatch_df["master_askap_non_conv_d2d"]=self.crossmatch_df["askap_non_conv_d2d"]
        else:
            self.crossmatch_df["master_scaled_non_conv_askap_iflux"]=0.0
            self.crossmatch_df["master_askap_non_conv_d2d"]=0.0
        self.crossmatch_df["master_scaled_catalog_iflux"]=self.crossmatch_df["{0}_{0}_scaled_to_askap".format(self.base_catalogue_name)]
        # self.crossmatch_df["master_scaled_catflux"]=self.crossmatch_df["{0}_{0}_scaled_to_askap".format(self.base_catalogue_name)]
        self.crossmatch_df["master_askap_snr"]=self.crossmatch_df["askap_askap_snr"]
        self.crossmatch_df["master_cat_snr"]=self.crossmatch_df["{0}_{0}_snr".format(self.base_catalogue_name)]
        # self.crossmatch_df["master_nvss_snr"]=self.crossmatch_df["nvss_nvss_snr".]
        self.crossmatch_df["master_askap_to_cat_scaled_snr"]=self.crossmatch_df["askap_{}_snr".format(self.base_catalogue_name)]
        self.crossmatch_df["master_cat_to_askap_scaled_snr"]=self.crossmatch_df["{0}_{0}_askap_snr".format(self.base_catalogue_name)]
        self.crossmatch_df["master_catalog_Mosaic"]=self.crossmatch_df["{0}_Mosaic".format(self.base_catalogue_name)]
        self.crossmatch_df["master_catalog_Mosaic_path"]=self.crossmatch_df["{0}_Mosaic_path".format(self.base_catalogue_name)]
        self.crossmatch_df["master_catalog_Mosaic_rms"]=self.crossmatch_df["{0}_Mosaic_rms".format(self.base_catalogue_name)]
        
        
        # self.goodmatches_df = self.crossmatch_df[self.crossmatch_df["d2d"]<=maxsep].reset_index(drop=True)
        # self.badmatches_df = self.crossmatch_df[self.crossmatch_df["d2d"]>maxsep].reset_index(drop=True)
        
        self.performed=True
        self.logger.info("Crossmatch complete.")
        
    def get_good_matches(self, maxsep):
        self.goodmatches_df = self.crossmatch_df[self.crossmatch_df["d2d"]<=maxsep].reset_index(drop=True)
        # return self.goodmatches_df
        
    def get_bad_matches(self, maxsep):
        self.badmatches_df = self.crossmatch_df[self.crossmatch_df["d2d"]>maxsep].reset_index(drop=True)
        # return self.bad_df
        
    def write_crossmatch(self, outname):
        if self.performed:
            self.crossmatch_df.to_csv(outname, index=False)
            self.logger.info("Written crossmatch dataframe to {}.".format(outname))
        else:
            self.logger.error("Need to perform a cross match first!")
        
    def calculate_ratio(self, col1, col2, output_col_name, col1_scaling=1., col2_scaling=1., dualmode=False, basecat="sumss"):
        self.crossmatch_df[output_col_name]=(self.crossmatch_df[col1]*col1_scaling)/(self.crossmatch_df[col2]*col2_scaling)
                
            # self.crossmatch_df[output_col_name]=self.crossmatch_df.fillna((self.crossmatch_df[col1]*col1_scaling)/(self.crossmatch_df[col2]*col2_scaling))
        
    def calculate_diff(self, col1, col2, output_col_name, col1_scaling=1., col2_scaling=1., dualmode=False, basecat="sumss"):
        self.crossmatch_df[output_col_name]=(self.crossmatch_df[col1]*col1_scaling)-(self.crossmatch_df[col2]*col2_scaling)
        
    def calculate_offsets(self, tag="catalog"):
        self.logger.info("Calculating offsets.")
        askap_coords = SkyCoord(ra=self.crossmatch_df["askap_ra"].values*u.deg, dec=self.crossmatch_df["askap_dec"].values*u.deg)
        catalog_coords = SkyCoord(ra=self.crossmatch_df["master_ra"].values*u.deg, dec=self.crossmatch_df["master_dec"].values*u.deg)
        dra, ddec = askap_coords.spherical_offsets_to(catalog_coords)
        self.crossmatch_df["askap_{}_ra_offset".format(tag)] = dra.arcsec
        self.crossmatch_df["askap_{}_dec_offset".format(tag)] = ddec.arcsec
            # self.crossmatch_df[output_col_name]= self.crossmatch_df.fillna((self.crossmatch_df[col1]*col1_scaling)-(self.crossmatch_df[col2]*col2_scaling))
                    
    def _plotinitial(self, panels, key, figure, image):
        panels[key]=aplpy.FITSFigure(image, figure=figure, subplot=(1,2,key+1))
        panels[key].show_grayscale()
            # panels[key].show_contour(images[n-1], colors='red', levels=[3.*12e-6, 4.*12.e-6, 8.*12.e-6, 16*12.e-6])
        panels[key].set_theme('publication')
        return panels
    
    def produce_postage_stamps(self, askap_data, askap_wcs, selection, nprocs, radius=13./60., convolve=False,
            askap_nonconv_image=None, askap_pre_convolve_catalog=None, dualmode=False, 
            basecat="sumss"):
        postage_stamps.crossmatch_stamps(self, askap_data, askap_wcs, selection, nprocs, radius=radius, convolve=convolve,
            askap_nonconv_img=askap_nonconv_image, askap_pre_convolve_catalog=askap_pre_convolve_catalog, dualmode=dualmode, 
            basecat=basecat)
        

    def merge(self, to_merge_crossmatch, basecat, max_sep=15., raw=False):
        if raw:
            self.crossmatch_df=self.crossmatch_df.append(to_merge_crossmatch.crossmatch_df)
        else:
            #Steps needed to be done:
            #1. Add to the crossmatch_df those ASKAP sources that do not have a SUMSS match.
            #This will also generate the columns we need to update the others.
            self.logger.info("Merging crossmatches...")
            askap_sources_in_base=self.crossmatch_df["askap_name"].tolist()
            to_add = to_merge_crossmatch.crossmatch_df[~to_merge_crossmatch.crossmatch_df["askap_name"].isin(askap_sources_in_base)].reset_index(drop=True)
            len_to_add=len(to_add.index)
            self.logger.info("Number of sources newly matched to add: {}.".format(len_to_add))
            if len_to_add>0:
                self.crossmatch_df=self.crossmatch_df.append(to_add)
            else:
                for col in to_merge_crossmatch.columns:
                    if col not in self.crossmatch_df.columns:
                        self.crossmatch_df[col]=np.nan
            to_review = to_merge_crossmatch.crossmatch_df[to_merge_crossmatch.crossmatch_df["askap_name"].isin(askap_sources_in_base)]
    
            if basecat=="sumss":
                other="nvss"
            else:
                other="sumss"
            other_d2d="{}_d2d".format(other)
            # print other_d2d
            self.crossmatch_df[other_d2d]=np.nan
            for i, row in to_review.iterrows():
                thisaskap_name = row["askap_name"]
                # Find the rows in the base with this askap source
                indexes = self.crossmatch_df.index[self.crossmatch_df["askap_name"]==thisaskap_name].tolist()
                if len(indexes) > 1:
                    self.logger.warning("Double match found in base: {}. Will not merge information to match.".format(row["askap_name"]))
                    continue
                else:
                    the_index = indexes[0]
                    for c in self.crossmatch_df.columns:
                        if other in c:
                            self.crossmatch_df.at[the_index, c]=row[c]
                            self.crossmatch_df.at[the_index, other_d2d]=row["d2d"]
                    
        self.merge_catalog=to_merge_crossmatch.base_catalog
                        
        self.logger.info("Finished merging.")
                        
            
    def _check_for_large_island(self, master_name, askap_name, island_name, pre_conv_crossmatch, askap_cat_islands_df=pd.DataFrame([]), 
        non_convolved_isl_cat_df=pd.DataFrame([]), threshold=3, askap_only=False, transient_sep=45.0):
        if len(askap_cat_islands_df.index) == 0:
            return False
        else:
            #First check the island in the crossmatched image islands
            n_components=askap_cat_islands_df[askap_cat_islands_df["island_id"]==island_name].iloc[0]["n_components"]
            self.logger.debug("First n_component {}".format(n_components))
            if n_components >= threshold:
                self.logger.info("Large island found for source {}".format(askap_name))
                return True
            elif non_convolved_isl_cat_df is not None:
                if len(non_convolved_isl_cat_df.index) == 0:
                    return False
            elif askap_only:
                #We need to find the nearest non-convovled source to the convolved askap source
                askap_source = self.comp_catalog.df[self.comp_catalog.df["name"]==askap_name]
                askap_ra = askap_source.iloc[0][self.comp_catalog.ra_col]
                askap_dec = askap_source.iloc[0][self.comp_catalog.dec_col]
                askap_target=SkyCoord(ra=askap_ra*u.deg, dec=askap_dec*u.deg)
                seps = askap_target.separation(pre_conv_crossmatch.comp_catalog._crossmatch_catalog)
                min_index = np.argmin(seps.deg)
                # print d2d_sumss.deg[min_index]
                # print min_index
                matching_source = pre_conv_crossmatch.comp_catalog.df.iloc[min_index]
                self.logger.debug("matching source: {}".format(matching_source))
                matching_island = matching_source["island"]
                n_components=non_convolved_isl_cat_df[non_convolved_isl_cat_df["island_id"]==matching_island].iloc[0]["n_components"]
                self.logger.debug("Second n_component {}".format(n_components))
                if n_components >= threshold:
                    self.logger.info("Large island found for source {}".format(askap_name))
                    return True
                else:
                    return False
            else:
                if non_convolved_isl_cat_df is not None:
                    #Now need to get the closest match to the master source in the pre_conv_crossmatch
                    non_conv_match=pre_conv_crossmatch.crossmatch_df[pre_conv_crossmatch.crossmatch_df["master_name"]==master_name].iloc[0]
                    pre_conv_island=non_conv_match["askap_island"]
                    if non_conv_match["d2d"] > transient_sep:
                        self.logger.debug("Closest preconvolved match for {} is not within acceptable range ({}). Returning False.".format(non_conv_match["master_name"], transient_sep))
                        return False
                    n_components=non_convolved_isl_cat_df[non_convolved_isl_cat_df["island_id"]==pre_conv_island].iloc[0]["n_components"]
                    self.logger.debug("Second n_component {}".format(n_components))
                    if n_components >= threshold:
                        self.logger.info("Large island found for source {}".format(askap_name))
                        return True
                    else:
                        return False
                else:
                    return False
        
    
    def _find_nearest_sources(self, ra, dec, source_name, sources_to_match, num_of_sources=6):
        target=SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        #Create a crossmatch list from the DF and remove the source itself
        sources_to_match=sources_to_match[sources_to_match["master_name"]!=source_name]
        sources_to_match_coords=SkyCoord(ra=sources_to_match["master_ra"].values*u.degree, dec=sources_to_match["master_dec"].values*u.degree)
        
        seps = target.separation(sources_to_match_coords)
        
        lowest=sorted(seps)[:num_of_sources]
        indexes=[]
        
        for i,val in enumerate(seps):
            if val in lowest:
                indexes.append(i)
        
        matched_sources = sources_to_match.iloc[indexes,:].reset_index(drop=True)
        
        return matched_sources
    
    def transient_search(self, max_separation=15.0, askap_snr_thresh=5.0, large_flux_thresh=3.0, pre_conv_crossmatch=None, image_beam_maj=45., image_beam_min=45., 
        image_beam_pa=0.0, dualmode=False, sumss=False, nvss=False, askap_img_wcs="None", askap_img_header={}, askap_img_data=[], preconv_askap_img_wcs="None", preconv_askap_img_header={},
        preconv_askap_img_data=[], askap_cat_islands_df=[], non_convolved_isl_cat_df=[], askap_image="None", preconv_askap_image="None", clean_for_sumss=False, max_sumss=0.0,
        clean_for_nvss=False, min_nvss=0.0):
        # Stage 1 - Define SUMSS sources with no match to ASKAP.
        # Stage 2 - Find those that have large flux ratios.
        # Stage 3 - Find ASKAP sources above the SUMSS threshold that have not been matched.
        if pre_conv_crossmatch is None:
            convolve = False
        else:
            convolve = True
        self.logger.info("Starting transient analysis")
        self.logger.info("Creating matched and non-matched lists")
        # First sort out the matched and non-matched
        # Does not care about dual mode or not.
        no_matches=self.crossmatch_df[self.crossmatch_df["d2d"]>max_separation].reset_index(drop=True)
        
        matches=self.crossmatch_df[self.crossmatch_df["d2d"]<=max_separation].reset_index(drop=True)
        #Check for duplicate matches
        dupmask=matches.duplicated(subset="askap_name")
        dup_matches=matches[dupmask].reset_index(drop=True)
        num_dup_matches=len(dup_matches.index)
        if num_dup_matches==0:
            self.logger.info("No double ASKAP source matches found.")
            moved_names=[]
        else:
            self.logger.warning("Matches to the same ASKAP source found!")
            self.logger.warning("Will keep the best match and move other match to 'no match' list.")
            indexes_to_move=[]
            unique_doubles=dup_matches["askap_name"].unique()
            for name in unique_doubles:
                temp=matches[matches["askap_name"]==name]
                minid=temp["d2d"].idxmin()
                allindex=temp.index.values.astype(int)
                for j in allindex:
                    if j!=minid:
                        indexes_to_move.append(j)
            rows=matches.loc[indexes_to_move]
            no_matches = no_matches.append(rows, ignore_index=True)
            matches.drop(rows.index, inplace=True)  
            #We need to know which ones were moved to get the right image - this is a bit of a hacky fix unfortunately
            moved_names = rows["master_name"].tolist()
        #Check the pre-convolved catalogue (if available) to see if a match is found there
        #Also does not care about dualmode.
        if pre_conv_crossmatch != None:
            using_pre_conv = True
            self.logger.info("Checking pre-convolved catalogue to check for matches that have been convolved out")
            preconv_df=pre_conv_crossmatch.crossmatch_df
            #Get the master_names of sources that have an askap match > max_sep
            no_matches_names = no_matches["master_name"].tolist()
            #We need to check if the new preconvolved match is not a direct match to a convolved source that is already matched.
            already_matched_askap_preconv_name = matches["askap_non_conv_name"]
            #Select only these sources from the preconv catalog
            preconv_df_no_matches = preconv_df[preconv_df["master_name"].isin(no_matches_names)].reset_index(drop=True)
            #Get which ones of these are within the acceptable max separation and not already matched
            preconv_df_have_match = preconv_df_no_matches[(preconv_df_no_matches["d2d"]<max_separation) & (~preconv_df_no_matches["askap_name"].isin(already_matched_askap_preconv_name))].reset_index(drop=True)
            #If there are some then we need to replace the info in the good matches with the non-convolved data
            if len(preconv_df_have_match.index)>0:
                #We need to force measure these sources in the ASKAP convolved image
                aegean_to_extract_df = preconv_df_have_match.filter(["master_name", "master_ra", "master_dec", "master_catflux"])
                aegean_to_extract_df.columns=["name", "ra", "dec", "peak_flux"]
                aegean_to_extract_df["peak_flux"]=aegean_to_extract_df["peak_flux"]*1.e-3
                #force the source size to be that of the ASKAP beam (in this case same as the catalogue)
                aegean_to_extract_df["a"]=askap_image.bmaj*3600.
                aegean_to_extract_df["b"]=askap_image.bmin*3600.
                aegean_to_extract_df["pa"]=askap_image.bpa
                aegean_results = self.force_extract_aegean(aegean_to_extract_df, askap_image.image)
                # preconv_df_have_match.join(aegean_to_extract_df, rsuffix="aegean_convolved")
                preconv_df_have_match = self._merge_forced_aegean(preconv_df_have_match, aegean_results, tag="aegean_convolved")
                self.logger.info("{} no match sources have a match in the pre-convolved image.".format(len(preconv_df_have_match.index)))
                preconv_matches_names = preconv_df_have_match["master_name"].tolist()
                #Now need to update the crossmatched ASKAP values to that of the non-convolved catalog
                #Need to move these ones that do have a match to good matches
                # rows_to_move = preconv_df[preconv_df["master_name"].isin(preconv_matches_names)]
                # self.logger.info("Moving {} sources to matched list.".format(len(rows_to_move.index)))
                rows_to_drop = no_matches[no_matches["master_name"].isin(preconv_matches_names)]
                self.logger.info("Dropping {} sources from non-matched list.".format(len(rows_to_drop.index)))
                no_matches.drop(rows_to_drop.index, inplace=True)
                matches = matches.append(preconv_df_have_match, ignore_index=True)
                #We need to know which ones were moved to get the right image - this is a bit of a hacky fix unfortunately
                # conv_moved_names = rows["sumss_name"].tolist()
                #Update the goodmatches for injection later
                self.goodmatches_df=matches
                #We also need to update the main crossmatch_df - this might be used for postage stamps later
                rows_to_drop = self.crossmatch_df[self.crossmatch_df["master_name"].isin(preconv_matches_names)]
                self.crossmatch_df.drop(rows_to_drop.index, inplace=True)
                self.crossmatch_df=self.crossmatch_df.append(preconv_df_have_match, ignore_index=True)
            else:
                self.logger.info("No further matches found in the convolved image.")
        else:
            using_pre_conv = False
        matches = matches.reset_index(drop=True)
        no_matches = no_matches.reset_index(drop=True)          
        self.goodmatches_df_trans=matches
        self.badmatches_df_trans=no_matches        
        # Stage 1
        # Consider all SUMSS sources with matches > max_sep to be 'not matched'
        self.logger.info("Finding sources with no ASKAP matches...")
        #For moved double matched sources the tag needs to be changed from BAD -> GOOD to load the correct image.
        
        no_match_postage_stamps=["source_{}_postagestamps.jpg".format(val) for i,val in no_matches["master_name"].iteritems()]

    
        no_matches["postage_stamp"]=no_match_postage_stamps
        
        #Force extract fluxes using Aegean from the ASKAP image
        aegean_to_extract_df = no_matches.filter(["master_name", "master_ra", "master_dec", "master_catflux"])
        aegean_to_extract_df.columns=["name", "ra", "dec", "peak_flux"]
        aegean_to_extract_df["peak_flux"]=aegean_to_extract_df["peak_flux"]*1.e-3
        #force the source size to be that of the ASKAP beam (in this case same as the catalogue)
        aegean_to_extract_df["a"]=image_beam_maj
        aegean_to_extract_df["b"]=image_beam_min
        aegean_to_extract_df["pa"]=image_beam_pa
        
        aegean_results = self.force_extract_aegean(aegean_to_extract_df, askap_image.image)
        no_matches = self._merge_forced_aegean(no_matches, aegean_results, tag="aegean_convolved")
        # no_matches.join(aegean_results, rsuffix="aegean")
        
        if pre_conv_crossmatch is not None:
            aegean_to_extract_df["a"]=preconv_askap_image.bmaj*3600.
            aegean_to_extract_df["b"]=preconv_askap_image.bmin*3600.
            aegean_to_extract_df["pa"]=preconv_askap_image.bpa
            preconv_aegean_results = self.force_extract_aegean(aegean_to_extract_df, preconv_askap_image.image)
            no_matches= self._merge_forced_aegean(no_matches, preconv_aegean_results, tag="aegean_preconvolved")
            # no_matches.join(aegean_results, rsuffix="aegean_preconv")
        
        self.logger.info("Measuring ASKAP local rms...")
        # no_matches.to_csv("testing_no_matches.csv")
        #measure the peak flux at the position in the ASKAP image(s)
        # no_matches=self.measure_askap_image_peak_fluxes(no_matches, askap_img_wcs, askap_img_data, using_pre_conv, preconv_askap_img_wcs, preconv_askap_img_data)
        #also measure local rms
        no_matches=self.measure_askap_image_local_rms(no_matches, askap_img_wcs, askap_img_header, askap_img_data, using_pre_conv, preconv_askap_img_wcs, preconv_askap_img_header, preconv_askap_img_data)
        matches=self.measure_askap_image_local_rms(matches, askap_img_wcs, askap_img_header, askap_img_data, using_pre_conv, preconv_askap_img_wcs, preconv_askap_img_header, preconv_askap_img_data)
        
        #Now check the sources to assign guesstimate type
        #First check bright source proximity
        # try:
        if sumss or dualmode:
            self.logger.info("Loading SUMSS bright sources")
            sumss_bright_sources_df=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/sumss_bright_sources.csv"), delimiter=";")
            # sumss_bright_tomatch = SkyCoord(ra=sumss_bright_sources_df["_RAJ2000"].values*u.deg, dec=sumss_bright_sources_df["_DEJ2000"].values*u.deg)
        if nvss or dualmode:
            self.logger.info("Loading NVSS bright sources")
            nvss_bright_sources_df=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/nvss_bright_sources.csv"), delimiter=";")
            # nvss_bright_tomatch = SkyCoord(ra=nvss_bright_sources_df["_RAJ2000"].values*u.deg, dec=nvss_bright_sources_df["_DEJ2000"].values*u.deg)
        
        if dualmode:
            sumss_bright_sources_df.append(nvss_bright_sources_df)
            bright_sources_tomatch=SkyCoord(ra=sumss_bright_sources_df["_RAJ2000"].values*u.deg, dec=sumss_bright_sources_df["_DEJ2000"].values*u.deg)
        elif sumss:
            bright_sources_tomatch=SkyCoord(ra=sumss_bright_sources_df["_RAJ2000"].values*u.deg, dec=sumss_bright_sources_df["_DEJ2000"].values*u.deg)
        else:
            bright_sources_tomatch=SkyCoord(ra=nvss_bright_sources_df["_RAJ2000"].values*u.deg, dec=nvss_bright_sources_df["_DEJ2000"].values*u.deg)
        
        pipeline_tags=[]
        nearest_sources=[]
    
        for i, row in no_matches.iterrows():
            targets=SkyCoord(ra=row["master_ra"]*u.deg, dec=row["master_dec"]*u.deg)
            seps = targets.separation(bright_sources_tomatch)
            if using_pre_conv==True:
                local_noise_col = "measured_preconv_askap_local_rms"
            else:
                local_noise_col = "measured_askap_local_rms"
            if np.min(seps.arcminute) <= 10.:
                pipeline_tags.append("Likely artefact (bright source)")
                #If it passes the bright source, then check if it is extremely extended (but this shouldn't be linked to bright source in case that fails.)
            #Check the non-convolved island if available
            
            elif self._check_for_large_island(row["master_name"], row["askap_name"], row["askap_island"], pre_conv_crossmatch, askap_cat_islands_df=askap_cat_islands_df, 
                non_convolved_isl_cat_df=non_convolved_isl_cat_df, threshold=3):
                pipeline_tags.append("Large island source")
            
            #check for possible really bright extended
            elif row["aegean_convolved_int_flux"]/row["askap_rms"] > 50.:
                pipeline_tags.append("Likely bright extended structure")
            
            #Check if the source would actually be seen in the local rms

            elif row["aegean_convolved_local_rms"]!=0.0 and row["master_scaled_catalog_iflux"]*1.e-3/row["aegean_convolved_local_rms"] < 5.0:
                pipeline_tags.append("Local RMS too high")
                
            elif row["survey_used"]=="sumss":
                if (row["sumss_MajAxis"] > (1.75 * row["sumss_telescope_bmaj"])) or (row["sumss_MinAxis"] > (1.75 * row["sumss_telescope_bmin"])):
                    pipeline_tags.append("Likely diffuse/extended")
                #Edge of image check - NEED TO TIDY UP THIS SECTION OF THE CODE REPITITION
                elif self._check_for_edge_case(row["master_ra"], row["master_dec"], askap_img_wcs, askap_img_data):
                    pipeline_tags.append("Edge of ASKAP image")
                else:
                    pipeline_tags.append("Candidate")
            else:
                if (row["nvss_MajAxis"] > (1.75 * 45.)) or (row["nvss_MinAxis"] > (1.75 * 45.)):
                    pipeline_tags.append("Likely diffuse/extended")
                elif self._check_for_edge_case(row["master_ra"], row["master_dec"], askap_img_wcs, askap_img_data):
                    pipeline_tags.append("Edge of ASKAP image")
                else:
                    pipeline_tags.append("Candidate")
            # print d2d_sumss.deg[min_index]
            # print min_index
            nearest_good_match_sources=self._find_nearest_sources(row["master_ra"], row["master_dec"], row["master_name"], self.goodmatches_df_trans)
            nearest_sources.append(",".join(nearest_good_match_sources["master_name"].tolist()))
        no_matches["pipelinetag"]=pipeline_tags
        no_matches["nearest_sources"]=nearest_sources
        # except:
            # self.logger.error("Bright source data cannot be found!")
            # self.logger.info("Skipping pipeline transient checks.")
            # no_matches["pipelinetag"]="N/A"
        
        self.transients_no_matches_df=no_matches
        no_matches.to_csv("transients_sources_no_askap_match.csv", index=False)
        self.logger.info("Written 'transients_sources_no_askap_match.csv'.")
        
        #New Stage 2 Check good matches
        self.logger.info("Checking good match sources...")
        
        goodmatch_postage_stamps=[]
        goodmatch_pipelinetags=[]
        nearest_sources=[]
        
        for i,row in self.goodmatches_df_trans.iterrows():
            # if row["d2d"] <= max_separation:
            goodmatch_postage_stamps.append("source_{}_postagestamps.jpg".format(row["master_name"]))
            if self._check_for_edge_case(row["master_ra"], row["master_dec"], askap_img_wcs, askap_img_data):
                goodmatch_pipelinetags.append("Edge of ASKAP image")
            #Check if either is extended
            elif (row["askap_a"] > (1.75 * row["askap_telescope_bmaj"])) or (row["askap_b"] > (1.75 * row["askap_telescope_bmin"])):
                goodmatch_pipelinetags.append("Extended source")
            elif (row["{}_MajAxis".format(row["survey_used"])] > (1.75 * row["{}_telescope_bmaj".format(row["survey_used"])])) or (row["{}_MinAxis".format(row["survey_used"])] > (1.75 * row["{}_telescope_bmin".format(row["survey_used"])])):
                goodmatch_pipelinetags.append("Extended source")
                
            elif pre_conv_crossmatch !=None:
                if row["master_name"] in preconv_matches_names:
                    if self._check_for_large_island(row["master_name"], row["askap_name"], row["askap_island"], pre_conv_crossmatch, askap_cat_islands_df=non_convolved_isl_cat_df, 
                         threshold=3, transient_sep=max_separation):
                        goodmatch_pipelinetags.append("Large island source")
                    else:
                        goodmatch_pipelinetags.append("Candidate")
                else:
                    if self._check_for_large_island(row["master_name"], row["askap_name"], row["askap_island"], pre_conv_crossmatch, askap_cat_islands_df=askap_cat_islands_df, 
                        non_convolved_isl_cat_df=non_convolved_isl_cat_df, threshold=3, transient_sep=max_separation):
                        goodmatch_pipelinetags.append("Large island source")
                    else:
                        goodmatch_pipelinetags.append("Candidate")
            elif self._check_for_large_island(row["master_name"], row["askap_name"], row["askap_island"], pre_conv_crossmatch, askap_cat_islands_df=askap_cat_islands_df, 
                    threshold=3, transient_sep=max_separation):
                goodmatch_pipelinetags.append("Large island source")
            else:
                goodmatch_pipelinetags.append("Candidate")
            # else:
                # large_ratios_postage_stamps.append("source_{}_BAD_sidebyside.jpg".format(row["master_name"]))
            #Also find the nearest sources
            nearest_good_match_sources=self._find_nearest_sources(row["master_ra"], row["master_dec"], row["master_name"], self.goodmatches_df_trans)
            nearest_sources.append(",".join(nearest_good_match_sources["master_name"].tolist()))
        self.goodmatches_df_trans["postage_stamp"] = goodmatch_postage_stamps
        self.goodmatches_df_trans["pipelinetag"] = goodmatch_pipelinetags
        self.goodmatches_df_trans["nearest_sources"] = nearest_sources
        
        # self.transients_large_ratios_df=large_ratios
        self.goodmatches_df_trans.to_csv("transients_good_matches.csv", index=False)
        self.logger.info("Written 'transients_good_matches.csv'.")
        
        # # Stage 2
        # # Here we want actual matches but with 'large' flux ratios
        # self.logger.info("Finding matches with large integrated flux ratios...")
        # median_match_flux_ratio=matches["askap_cat_int_flux_ratio"].median()
        # std_match_flux_ratio=matches["askap_cat_int_flux_ratio"].std()
        # #For now define sources with 'large' ratios as being more than median +/- std.
        # lower_limit=median_match_flux_ratio-(large_flux_thresh*std_match_flux_ratio)
        # upper_limit=median_match_flux_ratio+(large_flux_thresh*std_match_flux_ratio)
        # large_ratios=matches[(matches["askap_cat_int_flux_ratio"]<(lower_limit)) |
        #     (matches["askap_cat_int_flux_ratio"]>(upper_limit))].reset_index(drop=True)
        # large_ratios_postage_stamps=[]
        # large_ratios_pipelinetags=[]
        # nearest_sources=[]
        # for i,row in large_ratios.iterrows():
        #     # if row["d2d"] <= max_separation:
        #     large_ratios_postage_stamps.append("source_{}_GOOD_sidebyside.jpg".format(row["master_name"]))
        #     #Check if either is extended
        #     if (row["askap_a"] > (1.75 * row["askap_telescope_bmaj"])) or (row["askap_b"] > (1.75 * row["askap_telescope_bmin"])):
        #         large_ratios_pipelinetags.append("Extended source")
        #     elif (row["{}_MajAxis".format(row["survey_used"])] > (1.75 * row["{}_telescope_bmaj".format(row["survey_used"])])) or (row["{}_MinAxis".format(row["survey_used"])] > (1.75 * row["{}_telescope_bmin".format(row["survey_used"])])):
        #         large_ratios_pipelinetags.append("Extended source")
        #     elif pre_conv_crossmatch !=None:
        #         if row["master_name"] in preconv_matches_names:
        #             if self._check_for_large_island(row["master_name"], row["askap_name"], row["askap_island"], pre_conv_crossmatch, askap_cat_islands_df=non_convolved_isl_cat_df,
        #                  threshold=3, transient_sep=max_separation):
        #                 large_ratios_pipelinetags.append("Large island source")
        #             else:
        #                 large_ratios_pipelinetags.append("Match (non-convolved)")
        #         else:
        #             if lower_limit < (row["askap_non_conv_int_flux"]*1.e3/row["master_scaled_catalog_iflux"]) < upper_limit:
        #                 large_ratios_pipelinetags.append("Convolved flux error")
        #             elif self._check_for_large_island(row["master_name"], row["askap_name"], row["askap_island"], pre_conv_crossmatch, askap_cat_islands_df=askap_cat_islands_df,
        #                 non_convolved_isl_cat_df=non_convolved_isl_cat_df, threshold=3, transient_sep=max_separation):
        #                 large_ratios_pipelinetags.append("Large island source")
        #             else:
        #                 large_ratios_pipelinetags.append("Match (convolved)")
        #     elif self._check_for_large_island(row["master_name"], row["askap_name"], row["askap_island"], pre_conv_crossmatch, askap_cat_islands_df=askap_cat_islands_df,
        #             threshold=3, transient_sep=max_separation):
        #         large_ratios_pipelinetags.append("Large island source")
        #     else:
        #         large_ratios_pipelinetags.append("Match (convolved)")
        #     # else:
        #         # large_ratios_postage_stamps.append("source_{}_BAD_sidebyside.jpg".format(row["master_name"]))
        #     #Also find the nearest sources
        #     nearest_good_match_sources=self._find_nearest_sources(row["master_ra"], row["master_dec"], row["master_name"], self.goodmatches_df_trans)
        #     nearest_sources.append(",".join(nearest_good_match_sources["master_name"].tolist()))
        # large_ratios["postage_stamp"] = large_ratios_postage_stamps
        # large_ratios["pipelinetag"] = large_ratios_pipelinetags
        # large_ratios["nearest_sources"] = nearest_sources
        
        # self.transients_large_ratios_df=large_ratios
        # large_ratios.to_csv("transients_sources_matched_large_flux_ratio.csv", index=False)
        # self.logger.info("Written 'transients_sources_matched_large_flux_ratio.csv'.")
        
        # Stage 3 
        # - Obtain a list of ASKAP sources that have been matched well.
        # - Go through ASKAP sources - if not in list check SUMSS S/N.
        # - If S/N means it should be seen, it is a candidate (order by S/N).
        # - Fetch which SUMSS image the ASKAP source should be in.
        # - Extract SUMSS peak (?) flux at that position +/- 1 pixel.
        self.logger.info("Finding non-matched ASKAP sources that should have been seen...")
        matched_askap_sources=matches["askap_name"].tolist()
        not_matched_askap_sources=self.comp_catalog.df[~self.comp_catalog.df.name.isin(matched_askap_sources)].reset_index(drop=True)
        if sumss and clean_for_sumss:
            self.logger.info("Removing sources above the SUMSS boundary ({} deg)".format(max_sumss))
            not_matched_askap_sources=not_matched_askap_sources[not_matched_askap_sources["dec"]<=(max_sumss)].reset_index(drop=True)
        if nvss and clean_for_nvss:
            self.logger.info("Removing sources below the NVSS boundary ({} deg)".format(min_sumss))
            not_matched_askap_sources=not_matched_askap_sources[not_matched_askap_sources["dec"]>=(min_nvss)].reset_index(drop=True)
        # Now get those sources with a SNR ratio above the user defined threshold #Note March 4 - switch to 5.0, too many candidates at 4.5. March 13 - Added user option.
        mask=[]
        snrs=[]
        for i, row in not_matched_askap_sources.iterrows():
            if row["dec"]<-30.0:
                if row["sumss_snr"]>=askap_snr_thresh:
                    mask.append(True)
                    snrs.append(row["sumss_snr"])
                else:
                    mask.append(False)
            else:
                if nvss:
                    if row["nvss_snr"]>=askap_snr_thresh:
                        mask.append(True)
                        snrs.append(row["nvss_snr"])
                    else:
                        mask.append(False)
                else:
                    mask.append(True)
                    snrs.append(row["sumss_snr"])
                    
        original_cols = not_matched_askap_sources.columns
        if sumss:
            not_matched_askap_sources["survey_used"]="sumss"
            not_matched_askap_sources["survey"]="sumss"
        else:
            not_matched_askap_sources["survey_used"]="nvss"
            not_matched_askap_sources["survey"]="nvss"
        # not_matched_askap_sources["survey_used"]=["sumss" if d < -30 else "nvss" for d in not_matched_askap_sources["dec"].tolist()]
        # not_matched_askap_sources_should_see=not_matched_askap_sources_should_seeap_sources[not_matched_askap_sources["sumss_snr"]>=askap_snr_thresh].reset_index(drop=True)
        not_matched_askap_sources_should_see=not_matched_askap_sources[mask].reset_index(drop=True)
        not_matched_askap_sources_should_see["cat_snr"]=snrs
        # find_sumss_image_and_flux(not_matched_askap_sources_should_see)
        askapnotseen_postage_stamps=["source_{}_postagestamps.jpg".format(val) for i,val in not_matched_askap_sources_should_see["name"].iteritems()]
        not_matched_askap_sources_should_see["postage_stamp"]=askapnotseen_postage_stamps
        
        #Measure fluxes using aegean
        
        #For this we want to source find on every unique catalog image required
        
        unique_mosaics = not_matched_askap_sources_should_see["catalog_Mosaic_path"].unique()
        for i,mos in enumerate(unique_mosaics):
            if mos.split("/")[-1].startswith("J"):
                this_nvss = False
            else:
                this_nvss = True
            self.logger.info("Source finding from mosaic {} ({}/{})".format(mos.split("/")[-1], i+1, len(unique_mosaics)))
            if mos.split("/")[-1]=="SKIP.FITS":
                self.logger.warning("Missing SUMSS FITS mosaic. Skipping.")
                continue
            mosaic = Askapimage(mos)
            mosaic.load_fits_header()
            if not mosaic._beam_loaded:
                if this_nvss is False:
                    mosaic.calculate_sumss_beam()
                    mosaic.bmaj=mosaic.img_sumss_bmaj
                    mosaic.bmin=mosaic.img_sumss_bmin
                    mosaic.bpa=0.0
                else:
                    mosaic.calculate_nvss_beam()
                    mosaic.bmaj=mosaic.img_sumss_bmaj
                    mosaic.bmin=mosaic.img_sumss_bmin
                    mosaic.bpa=0.0
                
            aegean_to_extract_df = not_matched_askap_sources_should_see[not_matched_askap_sources_should_see["catalog_Mosaic_path"]==mos]
            aegean_to_extract_df = aegean_to_extract_df.filter(["ra", "dec", "peak_flux"])
            #force the source size to be that of the ASKAP beam (in this case same as the catalogue)
            aegean_to_extract_df["a"]=mosaic.bmaj
            aegean_to_extract_df["b"]=mosaic.bmin
            aegean_to_extract_df["pa"]=mosaic.bpa
            aegean_results = self.force_extract_aegean(aegean_to_extract_df, mos, nvss_beam=this_nvss)
            not_matched_askap_sources_should_see_subset = self._merge_forced_aegean(not_matched_askap_sources_should_see[not_matched_askap_sources_should_see["catalog_Mosaic_path"]==mos], 
                aegean_results, tag="aegean_catalog", ra_col="ra", dec_col="dec")
            if i == 0:
                new_not_matched_askap_sources_should_see=not_matched_askap_sources_should_see.copy()
                for c in not_matched_askap_sources_should_see_subset.columns:
                    if c not in not_matched_askap_sources_should_see.columns:
                        new_not_matched_askap_sources_should_see[c]=np.nan
                new_not_matched_askap_sources_should_see.update(not_matched_askap_sources_should_see_subset)
            else:
                new_not_matched_askap_sources_should_see.update(not_matched_askap_sources_should_see_subset)
        
               
        if len(not_matched_askap_sources_should_see.index)>0:
            tojoin = new_not_matched_askap_sources_should_see.filter([c for c in new_not_matched_askap_sources_should_see.columns if c not in not_matched_askap_sources_should_see.columns])
            not_matched_askap_sources_should_see = not_matched_askap_sources_should_see.join(tojoin)
        
        #Meaure fluxes in catalog image
        # not_matched_askap_sources_should_see = self.measure_catalog_image_peak_fluxes(not_matched_askap_sources_should_see)
        not_matched_askap_sources_should_see=self.measure_askap_image_local_rms(not_matched_askap_sources_should_see, askap_img_wcs, askap_img_header, askap_img_data, using_pre_conv,
            preconv_askap_img_wcs, preconv_askap_img_header, preconv_askap_img_data, askap_only=True)
        
        #check that no nan sources made it through from the aegean extraction (means that the SUMSS 'source' is out of range)
        not_matched_askap_sources_should_see = not_matched_askap_sources_should_see[~not_matched_askap_sources_should_see["aegean_catalog_int_flux"].isna()].reset_index(drop=True)
        
        #Perform checks for these 'transients'
        #1. Doubles
        #2. Extended/Diffuse
        #3. Edge - need to work out how to do
        pipeline_tags=[]
        nearest_sources=[]
        #Get the askap catalogue to match to
        askap_sources_tomatch = self.comp_catalog._crossmatch_catalog
        for i,row in not_matched_askap_sources_should_see.iterrows():
            askap_target=SkyCoord(ra=row["ra"]*u.deg, dec=row["dec"]*u.deg)
            seps = askap_target.separation(askap_sources_tomatch)
            if sorted(seps.arcsecond)[1] <= 2.5*45.:
                pipeline_tags.append("Likely double")
            elif (row["a"] > (1.75 * image_beam_maj)) or (row["b"] > (1.75 * image_beam_min)):
                pipeline_tags.append("Likely diffuse/extended")
            elif self._check_for_large_island(row["name"], row["name"], row["island"], pre_conv_crossmatch, askap_cat_islands_df=askap_cat_islands_df, 
                non_convolved_isl_cat_df=non_convolved_isl_cat_df, threshold=3, askap_only=True, transient_sep=max_separation):
                pipeline_tags.append("Large island source")
            elif self._check_for_edge_case(askap_target.ra.degree, askap_target.dec.degree, askap_img_wcs, askap_img_data):
                pipeline_tags.append("Edge of ASKAP image")
            elif (pre_conv_crossmatch != None) and (row["non_conv_d2d"] < 10.) and (row["non_conv_askap_scaled_to_{}".format(row["survey"])]/row["catalog_Mosaic_rms"] < 5.):
                pipeline_tags.append("Convolved flux error")
            else:
                pipeline_tags.append("Candidate")
                
            nearest_good_match_sources=self._find_nearest_sources(row["ra"], row["dec"], row["name"], self.goodmatches_df_trans)
            nearest_sources.append(",".join(nearest_good_match_sources["master_name"].tolist()))
            
        not_matched_askap_sources_should_see["pipelinetag"]=pipeline_tags
        not_matched_askap_sources_should_see["master_name"]=not_matched_askap_sources_should_see["name"]
        not_matched_askap_sources_should_see["nearest_sources"]=nearest_sources
        not_matched_askap_sources_should_see["master_ra"]=not_matched_askap_sources_should_see["ra"]
        not_matched_askap_sources_should_see["master_dec"]=not_matched_askap_sources_should_see["dec"]
        not_matched_askap_sources_should_see["master_catalog_Mosaic_path"]=not_matched_askap_sources_should_see["catalog_Mosaic_path"]
        not_matched_askap_sources_should_see["master_catalog_Mosaic_rms"]=not_matched_askap_sources_should_see["catalog_Mosaic_rms"]
        not_matched_askap_sources_should_see["master_catalog_Mosaic"]=not_matched_askap_sources_should_see["catalog_Mosaic"]
        if convolve:
            not_matched_askap_sources_should_see["master_scaled_non_conv_askap_iflux"]=[row["non_conv_askap_scaled_to_{}".format(row["survey"])] for i,row in not_matched_askap_sources_should_see.iterrows()]
        else:
            not_matched_askap_sources_should_see["master_scaled_non_conv_askap_iflux"]=[0.0 for i in range(len(not_matched_askap_sources_should_see.index))]
        not_matched_askap_sources_should_see["cat_snr"]=[row["sumss_snr"] if row["survey"]=="sumss" else row["nvss_snr"] for i,row in not_matched_askap_sources_should_see.iterrows()]
        not_matched_askap_sources_should_see["master_scaled_to_catalog"]=[row["askap_scaled_to_sumss"] if row["survey"]=="sumss" else row["askap_scaled_to_nvss"] for i,row in not_matched_askap_sources_should_see.iterrows()]
        
        
        self.transients_not_matched_askap_should_see_df=not_matched_askap_sources_should_see
        new_cols = ["askap_{}".format(i) if i in original_cols else i for i in self.transients_not_matched_askap_should_see_df.columns]
        self.transients_not_matched_askap_should_see_df.columns=new_cols
        self.transients_not_matched_askap_should_see_df["survey_used"]=self.transients_not_matched_askap_should_see_df["askap_survey_used"]
        not_matched_askap_sources_should_see.to_csv("transients_askap_sources_no_match_should_be_seen.csv", index=False)
        self.logger.info("Written 'transients_askap_sources_no_match_should_be_seen.csv'.")
        
        #Saving totals
        self.transients_noaskapmatchtocatalog_total=len(self.transients_no_matches_df.index)
        self.transients_noaskapmatchtocatalog_candidates=len(self.transients_no_matches_df[self.transients_no_matches_df["pipelinetag"]=="Candidate"].index)
        self.transients_nocatalogmatchtoaskap_total=len(self.transients_not_matched_askap_should_see_df.index)
        self.transients_nocatalogmatchtoaskap_candidates=len(self.transients_not_matched_askap_should_see_df[self.transients_not_matched_askap_should_see_df["pipelinetag"]=="Candidate"].index)
        # self.transients_largeratio_total=len(self.transients_large_ratios_df.index)
        # self.transients_largeratio_candidates=len(self.transients_large_ratios_df[self.transients_large_ratios_df["pipelinetag"].str.contains("Match ")])
        self.transients_goodmatches_total=len(self.goodmatches_df_trans.index)
        
        #Create overall transient table - merge goodmatches, bad matches and askap should be seen
        
        self.goodmatches_df_trans["type"]="goodmatch"
        self.transients_no_matches_df["type"]="noaskapmatch"
        self.transients_not_matched_askap_should_see_df["type"]="nocatalogmatch"
        self.transients_master_df = self.goodmatches_df_trans.append(self.transients_no_matches_df).reset_index(drop=True)
        self.transients_master_df = self.transients_master_df.append(self.transients_not_matched_askap_should_see_df).reset_index(drop=True)
        
        #Now go through and sort the master table
        self.sort_master_transient_table(convolve=convolve)
        
        self.transients_master_total=len(self.transients_master_df.index)
        candidate_mask = ((self.transients_master_df["pipelinetag"]=="Candidate") & (self.transients_master_df["master_ratio"]>=2.0))
        flagged_mask = ((self.transients_master_df["pipelinetag"]!="Candidate") & (self.transients_master_df["master_ratio"]>=2.0))
        self.transients_master_candidates_total = len(self.transients_master_df[candidate_mask].index)
        self.transients_master_flagged_total = len(self.transients_master_df[flagged_mask].index)
        
        #Check for null entries and mark as failed:
        self.validate_master_table()
        
        self.transients_master_df.to_csv("transients_master.csv")
        self.logger.info("Written master transient table as 'transients_master.csv'.")
    
    def sort_master_transient_table(self, convolve=True):
        self.logger.info("Processing transient table...")
        #Generate the master ratio column along with the error
        #1. for good types, we just take askap scaled and SUMSS int + errors.
        #1b. but for matches only to non-convolved we need to scale measured flux and error
        #2. For noaskapmatch take SUMSS and scaled measured ASKAP
        #3. For nocatalogmatch take scaled ASKAP and SUMSS measured data
        self.transients_master_df["master_ratio"]=0.0
        self.transients_master_df["master_ratio_err"]=0.0
        self.transients_master_df["aegean_scaled_int_flux"]=0.0
        self.transients_master_df["aegean_scaled_int_flux_err"]=0.0
        self.transients_master_df["ratio_askap_flux"]=0.0
        self.transients_master_df["ratio_askap_flux_err"]=0.0
        self.transients_master_df["ratio_catalog_flux"]=0.0
        self.transients_master_df["ratio_catalog_flux_err"]=0.0
        self.transients_master_df["master_catalog_mosaic"]="N/A"
        self.transients_master_df["master_catalog_mosaic_path"]="N/A"
        self.transients_master_df["askap_flux_to_use"]=0.0
        self.transients_master_df["askap_flux_to_use_err"]=0.0
        self.transients_master_df["scaled_askap_flux_to_use"]=0.0
        self.transients_master_df["scaled_askap_flux_to_use_err"]=0.0
        self.transients_master_df["askap_flux_to_use_2"]=0.0
        self.transients_master_df["askap_flux_to_use_2_err"]=0.0
        self.transients_master_df["scaled_askap_flux_to_use_2"]=0.0
        self.transients_master_df["scaled_askap_flux_to_use_2_err"]=0.0
        self.transients_master_df["catalog_flux_to_use"]=0.0
        self.transients_master_df["catalog_flux_to_use_err"]=0.0
        self.transients_master_df["distance_from_centre"]=0.0
        self.transients_master_df["aegean_rms_used"]="False"
        self.transients_master_df["inflated_convolved_flux"]="False"
        askap_freq = self.comp_catalog.frequency
        self.transients_master_df["aegean_convolved_int_flux_scaled"]=self.transients_master_df["aegean_convolved_int_flux"].apply(self._calculate_scaled_flux, args=(askap_freq, 843.e6))
        self.transients_master_df["aegean_convolved_err_int_flux_scaled"]=self.transients_master_df["aegean_convolved_err_int_flux"].apply(self._calculate_scaled_flux, args=(askap_freq, 843.e6))
        self.transients_master_df["aegean_convolved_local_rms_scaled"]=self.transients_master_df["aegean_convolved_local_rms"].apply(self._calculate_scaled_flux, args=(askap_freq, 843.e6))
        if convolve:
            self.transients_master_df["aegean_preconvolved_int_flux_scaled"]=self.transients_master_df["aegean_preconvolved_int_flux"].apply(self._calculate_scaled_flux, args=(askap_freq, 843.e6))
            self.transients_master_df["aegean_preconvolved_err_int_flux_scaled"]=self.transients_master_df["aegean_preconvolved_err_int_flux"].apply(self._calculate_scaled_flux, args=(askap_freq, 843.e6))
            self.transients_master_df["aegean_preconvolved_local_rms_scaled"]=self.transients_master_df["aegean_preconvolved_local_rms"].apply(self._calculate_scaled_flux, args=(askap_freq, 843.e6))
        else:
            self.transients_master_df["aegean_preconvolved_int_flux_scaled"]=[0.0 for i in range(len(self.transients_master_df.index))]
            self.transients_master_df["aegean_preconvolved_err_int_flux_scaled"]=[0.0 for i in range(len(self.transients_master_df.index))]
            self.transients_master_df["aegean_preconvolved_local_rms_scaled"]=[0.0 for i in range(len(self.transients_master_df.index))]
        self.transients_master_df["aegean_catalog_int_flux_scaled"]=self.transients_master_df["aegean_catalog_int_flux"].apply(self._calculate_scaled_flux, args=(843.e6, askap_freq))
        self.transients_master_df["aegean_catalog_err_int_flux_scaled"]=self.transients_master_df["aegean_catalog_err_int_flux"].apply(self._calculate_scaled_flux, args=(843.e6, askap_freq))
        self.transients_master_df["aegean_catalog_local_rms_scaled"]=self.transients_master_df["aegean_catalog_local_rms"].apply(self._calculate_scaled_flux, args=(843.e6, askap_freq))
        
        for i,row in self.transients_master_df.iterrows():
            if row["type"] == "goodmatch":
                if not pd.isna(row["aegean_convolved_int_flux_scaled"]):
                    if row["aegean_convolved_int_flux_scaled"] < 0.0 or row["aegean_convolved_int_flux_scaled"] < 3.*row["aegean_convolved_local_rms_scaled"]:
                        flux_to_use = 3.* row["aegean_convolved_local_rms_scaled"]
                        err_to_use = row["aegean_convolved_local_rms_scaled"]
                        used_rms = True
                    else:
                        flux_to_use = row["aegean_convolved_int_flux_scaled"]
                        err_to_use = row["aegean_convolved_err_int_flux_scaled"]
                        used_rms = False
                    askap_flux_to_use = row["aegean_convolved_int_flux"]
                    askap_flux_to_use_err = row["aegean_convolved_err_int_flux"]
                    askap_flux_to_use_2 = row["askap_int_flux"]
                    askap_flux_to_use_2_err = row["askap_err_int_flux"]
                    scaled_askap_flux_to_use = row["aegean_convolved_int_flux_scaled"]
                    scaled_askap_flux_to_use_err = row["aegean_convolved_err_int_flux_scaled"]
                    scaled_askap_flux_to_use_2 = row["askap_askap_scaled_to_{}".format(row.survey_used)]
                    scaled_askap_flux_to_use_2_err = 0.0
                               
                else:
                    flux_to_use = row["askap_askap_scaled_to_{}".format(row.survey_used)]
                    err_to_use = row["askap_askap_scaled_to_{}_err".format(row.survey_used)]
                    askap_flux_to_use = row["askap_int_flux"]
                    askap_flux_to_use_err = row["askap_err_int_flux"]
                    scaled_askap_flux_to_use = row["askap_askap_scaled_to_{}".format(row.survey_used)]
                    scaled_askap_flux_to_use_err = row["askap_askap_scaled_to_{}_err".format(row.survey_used)]
                    used_rms = False
                    if convolve:
                        askap_flux_to_use_2 = row["askap_non_conv_int_flux"]
                        askap_flux_to_use_2_err = row["askap_non_conv_err_int_flux"]
                        scaled_askap_flux_to_use_2 = row["askap_non_conv_askap_scaled_to_{}".format(row.survey_used)]
                        scaled_askap_flux_to_use_2_err = 0.0
                if row.survey_used == "sumss":
                    other_flux_to_use = row["sumss_St"]*1.e-3
                    other_error_to_use = row["sumss_e_St"]*1.e-3
                else:
                    other_flux_to_use = row["nvss_S1.4"]*1.e-3
                    other_error_to_use = row["nvss_e_S1.4"]*1.e-3
                catalog_flux_to_use = other_flux_to_use
                catalog_flux_to_use_err = other_error_to_use
                distance_from_centre = row["askap_distance_from_centre"]
                        
            elif row["type"] == "noaskapmatch":
                if row["aegean_convolved_int_flux_scaled"] < 0.0 or row["aegean_convolved_int_flux_scaled"] < 3.*row["aegean_convolved_local_rms_scaled"]:
                    flux_to_use = 3.* row["aegean_convolved_local_rms_scaled"]
                    err_to_use = row["aegean_convolved_local_rms_scaled"]
                    used_rms = True
                else:
                    flux_to_use = row["aegean_convolved_int_flux_scaled"]
                    err_to_use = row["aegean_convolved_err_int_flux_scaled"]
                    used_rms = False
                    
                askap_flux_to_use = row["aegean_convolved_int_flux"]
                askap_flux_to_use_err = row["aegean_convolved_err_int_flux"]
                if convolve:
                    askap_flux_to_use_2 = row["aegean_preconvolved_int_flux"]
                    askap_flux_to_use_2_err = row["aegean_preconvolved_err_int_flux"]
                
                scaled_askap_flux_to_use = row["aegean_convolved_int_flux_scaled"]
                scaled_askap_flux_to_use_err = row["aegean_convolved_err_int_flux_scaled"]
                if convolve:
                    scaled_askap_flux_to_use_2 = row["aegean_preconvolved_int_flux_scaled"]
                    scaled_askap_flux_to_use_2_err = row["aegean_preconvolved_err_int_flux_scaled"]
                
                if row.survey_used == "sumss":    
                    other_flux_to_use = row["sumss_St"]*1.e-3
                    other_error_to_use = row["sumss_e_St"]*1.e-3
                else:
                    other_flux_to_use = row["nvss_S1.4"]*1.e-3
                    other_error_to_use = row["nvss_e_S1.4"]*1.e-3
                catalog_flux_to_use = other_flux_to_use
                catalog_flux_to_use_err = other_error_to_use
                distance_from_centre = row["{}_distance_from_centre".format(row.survey_used)]
            
            elif row["type"] == "nocatalogmatch":
                if row["aegean_catalog_int_flux"] <= 0.0 or row["aegean_catalog_int_flux"] < 3.*row["aegean_catalog_local_rms"]:
                    other_flux_to_use = 3.* row["aegean_catalog_local_rms"]
                    other_error_to_use = row["aegean_catalog_local_rms"]
                    if other_flux_to_use ==  0.0:
                        if row["askap_dec"]<=-50:
                            other_flux_to_use = 3.*(6./5.)*1.e-3
                            other_error_to_use = (6./5.)*1.e-3                            
                        else:
                            other_flux_to_use = 3.*(10./5.)*1.e-3
                            other_error_to_use = (10./5.)*1.e-3
                    catalog_flux_to_use = other_flux_to_use
                    catalog_flux_to_use_err = other_error_to_use
                    used_rms = True
                else:
                    other_flux_to_use = row["aegean_catalog_int_flux"]
                    other_error_to_use = row["aegean_catalog_err_int_flux"]
                    catalog_flux_to_use = other_flux_to_use
                    catalog_flux_to_use_err = other_error_to_use
                    used_rms = False
                
                flux_to_use = row["askap_askap_scaled_to_{}".format(row.survey_used)]
                err_to_use = row["askap_askap_scaled_to_{}_err".format(row.survey_used)]
                askap_flux_to_use = row["askap_int_flux"]
                askap_flux_to_use_err = row["askap_err_int_flux"]
                if convolve:
                    askap_flux_to_use_2 = row["askap_non_conv_int_flux"]
                    askap_flux_to_use_2_err = row["askap_non_conv_err_int_flux"]
                scaled_askap_flux_to_use = row["askap_askap_scaled_to_{}".format(row.survey_used)]
                scaled_askap_flux_to_use_err = row["askap_askap_scaled_to_{}_err".format(row.survey_used)]
                if convolve:
                    scaled_askap_flux_to_use_2 = row["askap_non_conv_askap_scaled_to_{}".format(row.survey_used)]
                    scaled_askap_flux_to_use_2_err = 0.0
                distance_from_centre = row["askap_distance_from_centre"]
                
            
            self.logger.debug("type: {}, flux_to_use: {}, err_to_use: {}".format(row["type"], flux_to_use, err_to_use))
            self.logger.debug("type: {}, other_flux_to_use: {}, other_error_to_use: {}".format(row["type"], other_flux_to_use, other_error_to_use))
            
            #For candidate sources we need to check whether the non-convolved flux brings the ratio back down when convolved is used.
            #Apart from noaskapmatch
            convolve_check = ["goodmatch", "nocatalogmatch"]
            
            
            if flux_to_use > other_flux_to_use:
                ratio=flux_to_use/other_flux_to_use
                ratio_error = self._calculate_ratio_error(ratio, flux_to_use, err_to_use,
                    other_flux_to_use, other_error_to_use)
                if convolve and (row["type"] in convolve_check) and (row["pipelinetag"]=="Candidate") and (ratio >= 2.0):
                    try:
                        non_convolve_ratio = scaled_askap_flux_to_use_2 / other_flux_to_use
                    except:
                        non_convolve_ratio = 2.0
                    if non_convolve_ratio < 2.0:
                        inflated_var = "True"
                    else:
                        inflated_var = "False"
                else:
                    inflated_var = "False"
            else:
                ratio=other_flux_to_use/flux_to_use
                ratio_error = self._calculate_ratio_error(ratio, other_flux_to_use, other_error_to_use,
                    flux_to_use, err_to_use)
                if convolve and (row["type"] in convolve_check) and (row["pipelinetag"]=="Candidate") and (ratio >= 2.0):
                    try:
                        non_convolve_ratio =  other_flux_to_use / scaled_askap_flux_to_use_2
                    except:
                        non_convolve_ratio = 2.0
                    if non_convolve_ratio < 2.0:
                        inflated_var = "True"
                    else:
                        inflated_var = "False"
                else:
                    inflated_var = "False"
            
            self.logger.debug("ratio: {}, error: {}".format(ratio, ratio_error))
         
            self.transients_master_df.at[i, "inflated_convolved_flux"] = inflated_var
            self.transients_master_df.at[i, "master_ratio"] = ratio
            self.transients_master_df.at[i, "master_ratio_err"] = ratio_error
            self.transients_master_df.at[i, "ratio_askap_flux"] = flux_to_use
            self.transients_master_df.at[i, "ratio_askap_flux_err"] = err_to_use
            self.transients_master_df.at[i, "ratio_catalog_flux"] = other_flux_to_use
            self.transients_master_df.at[i, "ratio_catalog_flux_err"] = other_error_to_use
            self.transients_master_df.at[i, "askap_flux_to_use"] = askap_flux_to_use
            self.transients_master_df.at[i, "askap_flux_to_use_err"] = askap_flux_to_use_err
            self.transients_master_df.at[i, "scaled_askap_flux_to_use"] = scaled_askap_flux_to_use
            self.transients_master_df.at[i, "scaled_askap_flux_to_use_err"] = scaled_askap_flux_to_use_err
            if convolve:
                self.transients_master_df.at[i, "askap_flux_to_use_2"] = askap_flux_to_use_2
                self.transients_master_df.at[i, "askap_flux_to_use_2_err"] = askap_flux_to_use_2_err
                self.transients_master_df.at[i, "scaled_askap_flux_to_use_2"] = scaled_askap_flux_to_use_2
                self.transients_master_df.at[i, "scaled_askap_flux_to_use_2_err"] = scaled_askap_flux_to_use_2_err
            self.transients_master_df.at[i, "catalog_flux_to_use"] = catalog_flux_to_use
            self.transients_master_df.at[i, "catalog_flux_to_use_err"] = catalog_flux_to_use_err
            mosaic_col="master_catalog_Mosaic_path"
            self.transients_master_df.at[i, "master_catalog_mosaic"] = row[mosaic_col].split("/")[-1]
            self.transients_master_df.at[i, "master_catalog_mosaic_path"] = row[mosaic_col]
            self.transients_master_df.at[i, "distance_from_centre"] = distance_from_centre
            if used_rms:
                self.transients_master_df.at[i, "aegean_rms_used"] = "True"
            
        self.logger.info("Transient ratios calculated")
    
    def validate_master_table(self):
        #Here we want to check if any ratios are null meaning the source has failed to be extracted by Aegean
        null_indexes = self.transients_master_df[self.transients_master_df["master_ratio"].isnull()].index.tolist()
        if len(null_indexes)>0:
            self.logger.info("Some ratios failed, marking as 'Failed'.")
            for i in null_indexes:
                self.transients_master_df.at[i, "pipelinetag"]="Failed"
        else:
            self.logger.info("No ratio failures.")
            return
    
    def _calculate_ratio_error(self, ratio, x, dx, y, dy):
        return np.abs(ratio) * np.sqrt( ((dx/x)*(dx/x)) + ((dy/y)*(dy/y)))
    
    def _calculate_scaled_flux(self, flux, from_freq, to_freq, si=-0.8):
        if flux!=np.nan:
            return ((to_freq/from_freq)**(si))*flux
        else:
            return np.nan
         
    def _check_for_edge_case(self, ra, dec, img_wcs, img_data, num_pixels=50):
        #Works with ASKAP Image for now
        coord = SkyCoord(ra*u.degree, dec*u.degree)
        cutout = Cutout2D(img_data, coord, num_pixels*2., wcs=img_wcs)
        
        if np.isnan(cutout.data).any():
            return True
        
        pixels = img_wcs.wcs_world2pix(ra, dec, 1)
        
        size = img_data.shape
        
        if pixels[0] > size[0]-num_pixels:
            return True
        
        if pixels[1] > size[1]-num_pixels:
            return True
            
        return False
        
    def _get_peak_flux(self, ra, dec, img_wcs, img_data, num_pixels=10, sumss=False):
        #Works with ASKAP Image for now
        coord = SkyCoord(ra*u.degree, dec*u.degree)
        cutout = Cutout2D(img_data, coord, num_pixels*2., wcs=img_wcs)
        themax = np.max(cutout.data)
        if themax == np.nan:
            self.logger.error("Peak flux measuring returned nan - likely out of range.")
            self.logger.warning("Setting peak flux to 0.0")
            return 0.0
        return themax
        
    def _get_local_rms(self, ra, dec, img_wcs, img_data, num_pixels=300, sumss=False, max_iter=10, sigma=4):
        #Works with ASKAP Image for now
        coord = SkyCoord(ra*u.degree, dec*u.degree)
        cutout = Cutout2D(img_data, coord, num_pixels*2., wcs=img_wcs)
        cutout_clip_mean, cutout_clip_median, cutout_clip_std = sigma_clipped_stats(cutout.data, sigma=sigma, maxiters=max_iter)
        # return cutout_clip_std
        # pixels=img_wcs.wcs_world2pix(ra, dec, 1)
        # y,x = pixels
        # self.logger.debug("x:{} y:{}".format(x,y))
        # y = int(y)
        # x = int(x)
        # #Search 50 pixels around, see if in nan is in there
        # row_idx = np.array([range(x-num_pixels, x+num_pixels+1)])
        # col_idx = np.array([range(y-num_pixels, y+num_pixels+1)])
        # data_selection=img_data[0,0,row_idx[:, None], col_idx]
        # data_selection = np.nan_to_num(data_selection)
        # angle=angle
        # bscale=bscale
        # data_selection=data_selection.squeeze()
        # data_selection=data_selection*bscale
        # med, std, mask = self._Median_clip(data_selection, full_output=True, ftol=0.0, max_iter=max_iter, sigma=sigma)
        
        if cutout_clip_std == 0:
            self.logger.warning("Local RMS returned 0 or null. Setting to 1.0 mJy.")
            cutout_clip_std = 1.0e-3
        return cutout_clip_std
    
    def _Median_clip(self, arr, sigma=3, max_iter=3, ftol=0.01, xtol=0.05, full_output=False, axis=None):
        """Median_clip(arr, sigma, max_iter=3, ftol=0.01, xtol=0.05, full_output=False, axis=None)
        Return the median of an array after iteratively clipping the outliers.
        The median is calculated upon discarding elements that deviate more than
        sigma * standard deviation the median.

        arr: array to calculate the median from.
        sigma (3): the clipping threshold, in units of standard deviation.
        max_iter (3): the maximum number of iterations. A value of 0 will
            return the usual median.
        ftol (0.01): fraction tolerance limit for convergence. If the number
            of discarded elements changes by less than ftol, the iteration is
            stopped.
        xtol (0.05): absolute tolerance limit for convergence. If the number
            of discarded elements increases above xtol with respect to the
            initial number of elements, the iteration is stopped.
        full_output (False): If True, will also return the indices that were good.
        axis (None): Axis along which the calculation is to be done. NOT WORKING!!!

        >>> med = Median_clip(arr, sigma=3, max_iter=3)
        >>> med, std, inds_good = Median_clip(arr, sigma=3, max_iter=3, full_output=True)
        """
        arr = np.ma.masked_invalid(arr)
        med = np.median(arr, axis=axis)
        std = np.std(arr, axis=axis)
        ncount = arr.count(axis=axis)
        for niter in xrange(max_iter):
            ncount_old = arr.count(axis=axis)
            if axis is not None:
                condition = (arr < np.expand_dims(med-std*sigma, axis)) + (arr > np.expand_dims(med+std*sigma, axis))
            else:
                condition = (arr < med-std*sigma) + (arr > med+std*sigma)
            arr = np.ma.masked_where(condition, arr)
            ncount_new = arr.count(axis)
            if np.any(ncount-ncount_new > xtol*ncount):
                self.logger.debug("xtol reached {}; breaking at iteration {}".format(1-1.*ncount_new/ncount, niter+1))
                break
            if np.any(ncount_old-ncount_new < ftol*ncount_old):
                self.logger.debug("ftol reached {}; breaking at iteration {}".format(1-1.*ncount_new/ncount_old, niter+1) )
                break
            med = np.median(arr, axis=axis)
            std = np.std(arr, axis=axis)
        if full_output:
            if isinstance(arr.mask, np.bool_):
                mask = np.ones(arr.shape, dtype=bool)
            else:
                mask = ~arr.mask
            if axis is not None:
                med = med.data
                std = std.data
            return med, std, mask
        if axis is not None:
            med = med.data
        return med
    
    def measure_askap_image_peak_fluxes(self, no_matches, askap_img_wcs, askap_img_data, using_pre_conv, preconv_askap_img_wcs, preconv_askap_img_data):
        #This is ASKAP data so all the same image
        askap_peaks=[]
        preconv_askap_peaks=[]
        for i, row in no_matches.iterrows():
            askap_peak=self._get_peak_flux(row["master_ra"], row["master_dec"], askap_img_wcs, askap_img_data)
            if using_pre_conv:
                preconv_askap_peak=self._get_peak_flux(row["master_ra"], row["master_dec"], preconv_askap_img_wcs, preconv_askap_img_data)
            else:
                preconv_askap_peak=0.0
            askap_peaks.append(askap_peak)
            preconv_askap_peaks.append(preconv_askap_peak)
        
        no_matches["measured_askap_peak_flux"]=askap_peaks
        no_matches["measured_preconv_askap_peak_flux"]=preconv_askap_peaks
        
        return no_matches
        
    def measure_askap_image_local_rms(self, no_matches, askap_img_wcs, askap_img_header, askap_img_data, using_pre_conv, 
        preconv_askap_img_wcs, preconv_askap_img_header, preconv_askap_img_data, askap_only=False):
        #This is ASKAP data so all the same image
        if askap_only:
            ra_col="ra"
            dec_col="dec"
        else:
            ra_col="master_ra"
            dec_col="master_dec"
        askap_local_rms=[]
        preconv_askap_local_rms=[]
        for i, row in no_matches.iterrows():
            askap_local_rms_val=self._get_local_rms(row[ra_col], row[dec_col], askap_img_wcs, askap_img_data)
            if using_pre_conv:
                preconv_askap_local_rms_val=self._get_local_rms(row[ra_col], row[dec_col], preconv_askap_img_wcs, preconv_askap_img_data)
            else:
                preconv_askap_local_rms_val=0.0
            askap_local_rms.append(askap_local_rms_val)
            preconv_askap_local_rms.append(preconv_askap_local_rms_val)
        
        no_matches["measured_askap_local_rms"]=askap_local_rms
        no_matches["measured_preconv_askap_local_rms"]=preconv_askap_local_rms
        
        return no_matches
    
    def measure_catalog_image_peak_fluxes(self, no_matches, ra_col="ra", dec_col="dec"):
        #In this one we will have to load each fits image to do measure at the location
        loaded_fits={}
        peak_fluxes=[]
        for i,row in no_matches.iterrows():
            fits_image=row["catalog_Mosaic"]
            self.logger.debug(fits_image)
            if fits_image.startswith("J"):
                sumss_image=True
            else:
                sumss_image=False
            #Try using own library to load survey images
            if fits_image not in loaded_fits:
                loaded_fits[fits_image]=Askapimage(row["catalog_Mosaic_path"])
                loaded_fits[fits_image].load_wcs()
                loaded_fits[fits_image].load_fits_data()
            catalog_peak_flux = self._get_peak_flux(row[ra_col], row[dec_col], loaded_fits[fits_image].wcs, loaded_fits[fits_image].data, num_pixels=5, sumss=sumss_image)
            peak_fluxes.append(catalog_peak_flux)
        no_matches["measured_catalog_peak_flux"]=peak_fluxes
        
        return no_matches
    
    # def force_extract_pyse(self, df, image):
    #     #Step 1 - create pyse file
    #     num_of_sources = len(df.index)
    #     pyse_coords_file = "pyse_coordinates_to_measure_{}.json".format(image)
    #     with open(pyse_coords_file, "w") as f:
    #         f.write("[\n")
    #         for i,row in sources.iterrows():
    #             if i+1 == num_of_sources:
    #                 f.write("[{},{}]\n".format(row["RA"], row["Dec"]))
    #             else:
    #                 f.write("[{},{}],\n".format(row["RA"], row["Dec"]))
    #         f.write("]")
    #
    #     #Step 2 Run Pyse
    #
    #     command="pyse --force-beam --fixed-posns-file {} --csv --regions {} --ffbox 2.0".format(pyse_coords_file, image)
    #     subprocess.call(command, shell=True)
    #
    #     #Step 3 read output
    #     outputfile=image.replace(".fits", ".csv")
    #     results = pd.read_csv(outputfile, sep=", ")
    #     #Step 4 join to dataframe (they are in the same order)
    #     df.join(results, rsuffix="pyse")
    #     return df
    
    def force_extract_aegean(self, df, image, nvss_beam=False):
        #Step 1 - create aegean file
        num_of_sources = len(df.index)
        aegean_coords_file = "aegean_coordinates_to_measure_{}.csv".format(image.split("/")[-1].replace(".fits", ""))
        aegean_output_file = aegean_coords_file.replace(".csv", "_results.csv")
        
        #aegean needs ra, dec, int_flux, a, b, pa
        df.to_csv(aegean_coords_file, index=False)
            
        #Step 2 Run Aegean
            
        command="aegean --cores {} --priorized 1 --input {} --floodclip -1 --table {} {} --ratio 1".format(int(multiprocessing.cpu_count()/2), aegean_coords_file, aegean_output_file, image)
        if nvss_beam:
            command+=" --beam {0} {0} 0.0 --slice 0".format(45./3600.)
        subprocess.call(command, shell=True)
        
        #Step 3 read output
        #aegean adds a '_comp' to the output name for components
        aegean_output_file = aegean_output_file.replace(".csv", "_comp.csv")
        results = pd.read_csv(aegean_output_file)
        
        return results
        
    def _merge_forced_aegean(self, df, aegean_results, tag="aegean", ra_col="master_ra", dec_col="master_dec"):
        #AEGEAN SEEMS TO RANDOMISE THE ORDER THAT THE SOURCES ARE ENTERED INTO THE RESULTS FILE (PROBABLY THREADING)
        #Need to crossmatch the results
        aegean_results.columns=["{}_{}".format(tag, i) for i in aegean_results.columns]
        #If all sources are found then we can merge
        if len(aegean_results.index) == len(df.index):
            self.logger.info("All sources found in Aegean extraction, performing direct merge.")
            missing = False
        else:
            diff=len(df.index)-len(aegean_results.index)
            self.logger.warning("Not all sources found by aegean! {} are missing.".format(diff))
            missing = True
        df_match_catalog = SkyCoord(ra=df[ra_col].tolist()*u.degree, dec=df[dec_col].tolist()*u.degree)
        aegean_match_catalog = SkyCoord(ra=aegean_results["{}_ra".format(tag)].tolist()*u.degree, dec=aegean_results["{}_dec".format(tag)].tolist()*u.degree)
        idx, d2d, d3d = df_match_catalog.match_to_catalog_sky(aegean_match_catalog)
        if missing:
            self.logger.debug(d2d)
            indexes_missing=[df.index[i] for i,val in enumerate(d2d) if val.arcsec > 1.0]
            self.logger.debug(indexes_missing)
            for i in indexes_missing:
                temp = pd.Series([np.nan for j in aegean_results.columns], index=aegean_results.columns)
                temp.name = aegean_results.index[-1]+1
                aegean_results=aegean_results.append(temp,sort=False)
        else:
            indexes_missing = []
        #generate a dict to form the new index:
        indexes_map={}
        for i,val in enumerate(idx):
            if df.index[i] not in indexes_missing:
                indexes_map[val]=df.index[i]
        self.logger.debug(indexes_map)
        self.logger.debug(indexes_missing)
        new_index = [indexes_map[i] for i in sorted(indexes_map)] + indexes_missing
        aegean_results.index=new_index
        newdf = df.join(aegean_results)

        return newdf
                          
        
    def inject_transients_db(self, image_id, sumss, nvss, db_engine="postgresql", 
            db_username="postgres", db_password="postgres", db_host="localhost", db_port="5432", db_database="postgres", max_separation=45.0, dualmode=False, basecat="sumss", askap_image="image", askap_image_2="N/A"):        

        engine = sqlalchemy.create_engine('{}://{}:{}@{}:{}/{}'.format(db_engine, db_username, db_password, db_host, db_port, db_database))
        
        self.transients_master_df["image_id"]=image_id
        self.transients_master_df["usertag"]="N/A"
        self.transients_master_df["userreason"]="N/A"
        self.transients_master_df["checkedby"]="N/A"
        self.transients_master_df["askap_image"]=askap_image
        self.transients_master_df["askap_image_2"]=askap_image_2
        
        if nvss:
            survey = "nvss"
        else:
            survey = "sumss"
            
        if "askap_non_conv_d2d" not in self.transients_master_df.columns:
            self.transients_master_df["askap_non_conv_d2d"]=0.0
            
        if "askap_rms_preconv" not in self.transients_master_df.columns:
            self.transients_master_df["askap_rms_preconv"]=-1.0
            
        # db_col_names_order = ["image_id", #need to add
        #                 # "match_id", #need to add
        #                 "master_name", #done
        #                 "askap_name", #done
        #                 "catalog_name", #done
        #                 "ra", #done
        #                 "dec", #done
        #                 "catalog_iflux", #done
        #                 "catalog_iflux_e", #done
        #                 "catalog_rms", #done
        #                 # "catalog_local_rms", #done
        #                 "catalog_mosaic",
        #                 "askap_iflux", #need to add
        #                 "askap_iflux_e", #need to add
        #                 "askap_scale_flux", #need to add
        #                 "askap_scale_flux_e", #need to add
        #                 "askap_non_conv_flux",
        #                 "askap_non_conv_flux_e",
        #                 "askap_non_conv_scaled_flux",
        #                 "askap_non_conv_scaled_flux_e",
        #                 "askap_non_conv_d2d", #need to add
        #                 "d2d_askap_centre",   #need to add
        #                 "askap_rms",
        #                 "askap_rms_2",
        #                 "ratio",
        #                 "ratio_e",
        #                 "ratio_catalog_flux"
        #                 "ratio_catalog_flux_err"
        #                 "ratio_askap_flux"
        #                 "ratio_askap_flux_err"
        #                 "ploturl",  #done
        #                 "pipelinetag",  #done
        #                 "usertag", #need to add
        #                 "userreason", #need to add
        #                 "checkedby", #need to add
        #                 "survey",
        #                 "nearest_sources",
        #                 "type"] #done

        db_col_names_map = {"image_id":"image_id", #need to add
                        # "match_id", #need to add
                        "master_name":"master_name", #done
                        "askap_name":"askap_name",
                        "{}_name".format(survey):"catalog_name", #done
                        "master_ra":"ra", #done
                        "master_dec":"dec", #done
                        "catalog_flux_to_use":"catalog_iflux", #done
                        "catalog_flux_to_use_err":"catalog_iflux_e", #done
                        "master_catalog_Mosaic_rms":"catalog_rms", #done
                        # "catalog_local_rms", #done
                        "master_catalog_Mosaic":"catalog_mosaic",
                        "askap_flux_to_use":"askap_iflux", #need to add
                        "askap_flux_to_use_err":"askap_iflux_e", #need to add
                        "scaled_askap_flux_to_use":"askap_scale_flux", #need to add
                        "scaled_askap_flux_to_use_err":"askap_scale_flux_e", #need to add
                        "askap_flux_to_use_2":"askap_non_conv_flux",
                        "askap_flux_to_use_2_err":"askap_non_conv_flux_e",
                        "scaled_askap_flux_to_use_2":"askap_non_conv_scaled_flux",
                        "scaled_askap_flux_to_use_2_err":"askap_non_conv_scaled_flux_e",
                        "askap_non_conv_d2d":"askap_non_conv_d2d",
                        "distance_from_centre":"d2d_askap_centre",   #need to add
                        "askap_rms":"askap_rms",
                        "askap_rms_preconv":"askap_rms_2",
                        "master_ratio":"ratio",
                        "master_ratio_err":"ratio_e",
                        "ratio_catalog_flux":"ratio_catalog_flux",
                        "ratio_catalog_flux_err":"ratio_catalog_flux_err",
                        "ratio_askap_flux":"ratio_askap_flux",
                        "ratio_askap_flux_err":"ratio_askap_flux_err",
                        "postage_stamp":"ploturl",  #done
                        "pipelinetag":"pipelinetag",  #done
                        "usertag":"usertag", #need to add
                        "userreason":"userreason", #need to add
                        "checkedby":"checkedby", #need to add
                        "survey_used":"survey",
                        "nearest_sources":"nearest_sources",
                        "type":"transient_type",
                        "aegean_rms_used":"aegean_rms_used",
                        "askap_image":"askap_image",
                        "askap_image_2":"askap_image_2",
                        "inflated_convolved_flux":"inflated_convolved_flux"} #done
                        
        filter_list = [i for i in db_col_names_map]
        
        to_write = self.transients_master_df.filter(filter_list, axis=1)
        
        to_write.columns=[db_col_names_map[i] for i in to_write.columns]
        
        to_write["ploturl"]=[os.path.join("media/{}/stamps/{}".format(image_id, i)) for i in to_write["ploturl"].values]
        
        to_write = to_write.sort_values(by=["pipelinetag","ratio"], ascending=[True, False])
        new_type_vals = {"goodmatch":"Good match", "nocatalogmatch":"No catalog match", "noaskapmatch":"No askap match"}
        to_write["transient_type"]=[new_type_vals[i] for i in to_write["transient_type"].values]
        values = {"catalog_name":"N/A", "askap_name":"N/A", "askap_non_conv_d2d":0.0,
                    "catalog_iflux":-1.0, #done
                    "catalog_iflux_e":-1.0, #done
                    "askap_iflux":-1.0, #need to add
                    "askap_iflux_e":-1.0, #need to add
                    "askap_scale_flux":-1.0, #need to add
                    "askap_scale_flux_e":-1.0, #need to add
                    "askap_non_conv_flux":-1.0,
                    "askap_non_conv_flux_e":-1.0,
                    "askap_non_conv_scaled_flux":-1.0,
                    "askap_non_conv_scaled_flux_e":-1.0,
                    "ratio_catalog_flux":-1.0,
                    "ratio_catalog_flux_err":-1.0,
                    "ratio_askap_flux":-1.0,
                    "ratio_askap_flux_err":-1.0,
                    "ratio":-1.0,
                    "ratio_e":-1.0
                    }
        to_write.fillna(value=values, inplace=True)
        to_write.to_sql("images_transients", engine, if_exists="append", index=False)
        