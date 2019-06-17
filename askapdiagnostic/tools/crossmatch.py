#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import logging
import aplpy
from astropy.coordinates import SkyCoord
from astropy import units as u
from askapdiagnostic.tools.fitsimage import askapimage
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D
from askapdiagnostic.plotting import postage_stamps
import sqlalchemy
import psycopg2
import pkg_resources
import numpy as np
import pandas as pd

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
        self.matches=self.matches.rename(str, columns=newcolnames)
        self.crossmatch_df = self.base_catalog.df.copy(deep=True)
        newcolnames={}
        for c in self.crossmatch_df.columns:
            newcolnames[c]="{}_{}".format(self.base_catalogue_name, c)
        self.crossmatch_df=self.crossmatch_df.rename(str, columns=newcolnames)
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
            # self.crossmatch_df[output_col_name]= self.crossmatch_df.fillna((self.crossmatch_df[col1]*col1_scaling)-(self.crossmatch_df[col2]*col2_scaling))
                    
    def _plotinitial(self, panels, key, figure, image):
        panels[key]=aplpy.FITSFigure(image, figure=figure, subplot=(1,2,key+1))
        panels[key].show_grayscale()
            # panels[key].show_contour(images[n-1], colors='red', levels=[3.*12e-6, 4.*12.e-6, 8.*12.e-6, 16*12.e-6])
        panels[key].set_theme('publication')
        return panels
    
    def produce_postage_stamps(self, postage_options, radius=15./60., nprocs=1, max_separation=15., convolve=False, pre_convolve_image=None, preconolve_catalog=None,
        dualmode=False, basecat="sumss", transients=False):
        askap_fits = self.crossmatch_df["askap_image"].iloc[0]
        postage_stamps.crossmatch_stamps(self, askap_fits, postage_options, nprocs, radius=13./60., max_separation=max_separation, convolve=convolve,
            askap_pre_convolve_image=pre_convolve_image, askap_pre_convolve_catalog=preconolve_catalog, dualmode=dualmode, basecat=basecat, transients=transients)
        

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
            elif len(non_convolved_isl_cat_df.index) == 0:
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
        dualmode=False, sumss=False, nvss=False, askap_img_wcs="None", askap_img_header={}, askap_img_data=[], preconv_askap_img_wcs="None", preconv_askap_img_header={},
        preconv_askap_img_data=[], askap_cat_islands_df=[], non_convolved_isl_cat_df=[]):
        # Stage 1 - Define SUMSS sources with no match to ASKAP.
        # Stage 2 - Find those that have large flux ratios.
        # Stage 3 - Find ASKAP sources above the SUMSS threshold that have not been matched.
        
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
            #Select only these sources from the preconv catalog
            preconv_df_no_matches = preconv_df[preconv_df["master_name"].isin(no_matches_names)].reset_index(drop=True)
            #Get which ones of these are within the acceptable max separation
            preconv_df_have_match = preconv_df_no_matches[preconv_df_no_matches["d2d"]<max_separation].reset_index(drop=True)
            #If there are some then we need to replace the info in the good matches with the non-convolved data
            if len(preconv_df_have_match.index)>0:
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
                   
        self.goodmatches_df_trans=matches
        self.badmatches_df_trans=no_matches        
        # Stage 1
        # Consider all SUMSS sources with matches > max_sep to be 'not matched'
        self.logger.info("Finding sources with no ASKAP matches...")
        #For moved double matched sources the tag needs to be changed from BAD -> GOOD to load the correct image.
        if len(moved_names)==0:
            no_match_postage_stamps=["source_{}_BAD_sidebyside.jpg".format(val) for i,val in no_matches["sumss_name"].iteritems()]
        else:
            no_match_postage_stamps=[]
            for i,val in no_matches["master_name"].iteritems():
                # if val not in moved_names:
                no_match_postage_stamps.append("source_{}_BAD_sidebyside.jpg".format(val))
                # else:
                    # no_match_postage_stamps.append("source_{}_GOOD_sidebyside.jpg".format(val))
    
        no_matches["postage_stamp"]=no_match_postage_stamps
        
        #measure the peak flux at the position in the ASKAP image(s)
        no_matches=self.measure_askap_image_peak_fluxes(no_matches, askap_img_wcs, askap_img_data, using_pre_conv, preconv_askap_img_wcs, preconv_askap_img_data)
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
            elif row["measured_askap_peak_flux"]/row["askap_rms"] > 50.:
                pipeline_tags.append("Likely bright extended structure")
            
            #Check if the source would actually be seen in the local rms
            #Due to the excessive RMS in the convolved case, we use non-convovled if available
            elif row["master_scaled_catalog_iflux"]*1.e-3/row[local_noise_col] < 5.0:
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
        
        # Stage 2
        # Here we want actual matches but with 'large' flux ratios
        self.logger.info("Finding matches with large integrated flux ratios...")
        median_match_flux_ratio=matches["askap_cat_int_flux_ratio"].median()
        std_match_flux_ratio=matches["askap_cat_int_flux_ratio"].std()
        #For now define sources with 'large' ratios as being more than median +/- std.
        lower_limit=median_match_flux_ratio-(large_flux_thresh*std_match_flux_ratio)
        upper_limit=median_match_flux_ratio+(large_flux_thresh*std_match_flux_ratio)
        large_ratios=matches[(matches["askap_cat_int_flux_ratio"]<(lower_limit)) | 
            (matches["askap_cat_int_flux_ratio"]>(upper_limit))].reset_index(drop=True)
        large_ratios_postage_stamps=[]
        large_ratios_pipelinetags=[]
        nearest_sources=[]
        for i,row in large_ratios.iterrows():
            # if row["d2d"] <= max_separation:
            large_ratios_postage_stamps.append("source_{}_GOOD_sidebyside.jpg".format(row["master_name"]))
            #Check if either is extended
            if (row["askap_a"] > (1.75 * row["askap_telescope_bmaj"])) or (row["askap_b"] > (1.75 * row["askap_telescope_bmin"])):
                large_ratios_pipelinetags.append("Extended source")
            elif (row["{}_MajAxis".format(row["survey_used"])] > (1.75 * row["{}_telescope_bmaj".format(row["survey_used"])])) or (row["{}_MinAxis".format(row["survey_used"])] > (1.75 * row["{}_telescope_bmin".format(row["survey_used"])])):
                large_ratios_pipelinetags.append("Extended source")
            elif pre_conv_crossmatch !=None:
                if row["master_name"] in preconv_matches_names:
                    if self._check_for_large_island(row["master_name"], row["askap_name"], row["askap_island"], pre_conv_crossmatch, askap_cat_islands_df=non_convolved_isl_cat_df, 
                         threshold=3, transient_sep=max_separation):
                        large_ratios_pipelinetags.append("Large island source")
                    else:
                        large_ratios_pipelinetags.append("Match (non-convolved)")
                else:
                    if lower_limit < (row["askap_non_conv_int_flux"]*1.e3/row["master_scaled_catalog_iflux"]) < upper_limit:
                        large_ratios_pipelinetags.append("Convolved flux error")
                    elif self._check_for_large_island(row["master_name"], row["askap_name"], row["askap_island"], pre_conv_crossmatch, askap_cat_islands_df=askap_cat_islands_df, 
                        non_convolved_isl_cat_df=non_convolved_isl_cat_df, threshold=3, transient_sep=max_separation):
                        large_ratios_pipelinetags.append("Large island source")
                    else:
                        large_ratios_pipelinetags.append("Match (convolved)")
            elif self._check_for_large_island(row["master_name"], row["askap_name"], row["askap_island"], pre_conv_crossmatch, askap_cat_islands_df=askap_cat_islands_df, 
                    threshold=3, transient_sep=max_separation):
                large_ratios_pipelinetags.append("Large island source")
            else:
                large_ratios_pipelinetags.append("Match (convolved)")
            # else:
                # large_ratios_postage_stamps.append("source_{}_BAD_sidebyside.jpg".format(row["master_name"]))
            #Also find the nearest sources
            nearest_good_match_sources=self._find_nearest_sources(row["master_ra"], row["master_dec"], row["master_name"], self.goodmatches_df_trans)
            nearest_sources.append(",".join(nearest_good_match_sources["master_name"].tolist()))
        large_ratios["postage_stamp"] = large_ratios_postage_stamps
        large_ratios["pipelinetag"] = large_ratios_pipelinetags
        large_ratios["nearest_sources"] = nearest_sources
        
        self.transients_large_ratios_df=large_ratios
        large_ratios.to_csv("transients_sources_matched_large_flux_ratio.csv", index=False)
        self.logger.info("Written 'transients_sources_matched_large_flux_ratio.csv'.")
        
        # Stage 3 
        # - Obtain a list of ASKAP sources that have been matched well.
        # - Go through ASKAP sources - if not in list check SUMSS S/N.
        # - If S/N means it should be seen, it is a candidate (order by S/N).
        # - Fetch watch SUMSS image the ASKAP source should be in.
        # - Extract SUMSS peak (?) flux at that position +/- 1 pixel.
        self.logger.info("Finding non-matched ASKAP sources that should have been seen...")
        matched_askap_sources=matches["askap_name"].tolist()
        not_matched_askap_sources=self.comp_catalog.df[~self.comp_catalog.df.name.isin(matched_askap_sources)].reset_index(drop=True)
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
                if row["nvss_snr"]>=askap_snr_thresh:
                    mask.append(True)
                    snrs.append(row["nvss_snr"])
                else:
                    mask.append(False)
        not_matched_askap_sources["survey"]=["sumss" if d < -30 else "nvss" for d in not_matched_askap_sources["dec"].tolist()]
        # not_matched_askap_sources_should_see=not_matched_askap_sources_should_seeap_sources[not_matched_askap_sources["sumss_snr"]>=askap_snr_thresh].reset_index(drop=True)
        not_matched_askap_sources_should_see=not_matched_askap_sources[mask].reset_index(drop=True)
        not_matched_askap_sources_should_see["cat_snr"]=snrs
        # find_sumss_image_and_flux(not_matched_askap_sources_should_see)
        askapnotseen_postage_stamps=["transient_askapnotseen_ASKAP_{}_sidebyside.jpg".format(val) for i,val in not_matched_askap_sources_should_see["name"].iteritems()]
        not_matched_askap_sources_should_see["postage_stamp"]=askapnotseen_postage_stamps
        
        #Meaure fluxes in catalog image
        not_matched_askap_sources_should_see = self.measure_catalog_image_peak_fluxes(not_matched_askap_sources_should_see)
        not_matched_askap_sources_should_see=self.measure_askap_image_local_rms(not_matched_askap_sources_should_see, askap_img_wcs, askap_img_header, askap_img_data, using_pre_conv, 
            preconv_askap_img_wcs, preconv_askap_img_header, preconv_askap_img_data, askap_only=True)
        
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
        not_matched_askap_sources_should_see["nearest_sources"]=nearest_sources
        not_matched_askap_sources_should_see["master_scaled_non_conv_askap_iflux"]=[row["non_conv_askap_scaled_to_{}".format(row["survey"])] for i,row in not_matched_askap_sources_should_see.iterrows()]
        not_matched_askap_sources_should_see["cat_snr"]=[row["sumss_snr"] if row["survey"]=="sumss" else row["nvss_snr"] for i,row in not_matched_askap_sources_should_see.iterrows()]
        not_matched_askap_sources_should_see["master_scaled_to_catalog"]=[row["askap_scaled_to_sumss"] if row["survey"]=="sumss" else row["askap_scaled_to_nvss"] for i,row in not_matched_askap_sources_should_see.iterrows()]
        
        
        self.transients_not_matched_askap_should_see_df=not_matched_askap_sources_should_see
        not_matched_askap_sources_should_see.to_csv("transients_askap_sources_no_match_should_be_seen.csv", index=False)
        self.logger.info("Written 'transients_askap_sources_no_match_should_be_seen.csv'.")
        
        #Saving totals
        self.transients_noaskapmatchtocatalog_total=len(self.transients_no_matches_df.index)
        self.transients_noaskapmatchtocatalog_candidates=len(self.transients_no_matches_df[self.transients_no_matches_df["pipelinetag"]=="Candidate"].index)
        self.transients_nocatalogmatchtoaskap_total=len(self.transients_not_matched_askap_should_see_df.index)
        self.transients_nocatalogmatchtoaskap_candidates=len(self.transients_not_matched_askap_should_see_df[self.transients_not_matched_askap_should_see_df["pipelinetag"]=="Candidate"].index)
        self.transients_largeratio_total=len(self.transients_large_ratios_df.index)
        self.transients_largeratio_candidates=len(self.transients_large_ratios_df[self.transients_large_ratios_df["pipelinetag"].str.contains("Match ")])
        self.transients_goodmatches_total=len(self.goodmatches_df_trans.index)
        
    def _check_for_edge_case(self, ra, dec, img_wcs, img_data, num_pixels=50):
        #Works with ASKAP Image for now
        pixels=img_wcs.wcs_world2pix(ra, dec, 1)
        y,x = pixels
        y = int(y)
        x = int(x)
        #Search 10 pixels around, see if in nan is in there
        row_idx = np.array([range(x-num_pixels, x+num_pixels+1)])
        col_idx = np.array([range(y-num_pixels, y+num_pixels+1)])
        data_selection=img_data[0,0,row_idx[:, None], col_idx]
        return np.isnan(data_selection).any()
        
    def _get_peak_flux(self, ra, dec, img_wcs, img_data, num_pixels=10, sumss=False):
        #Works with ASKAP Image for now
        self.logger.debug(ra)
        self.logger.debug(dec)
        pixels=img_wcs.wcs_world2pix(ra, dec, 1)
        y,x = pixels
        self.logger.debug("x:{} y:{}".format(x,y))
        y = int(y)
        x = int(x)
        #Search 10 pixels around, see if in nan is in there
        row_idx = np.array([range(x-num_pixels, x+num_pixels+1)])
        col_idx = np.array([range(y-num_pixels, y+num_pixels+1)])
        try:
            if sumss:
                data_selection=img_data[row_idx[:, None], col_idx]
            else:
                data_selection=img_data[0,0,row_idx[:, None], col_idx]
        except:
            self.logger.error("Peak flux measuring failed - likely out of range.")
            self.logger.warning("Setting peak flux to 0.0")
            data_selection=[0.0]
        return np.max(data_selection)
        
    def _get_local_rms(self, ra, dec, img_wcs, img_data, angle, bscale, num_pixels=300, sumss=False, max_iter=10, sigma=4):
        #Works with ASKAP Image for now
        pixels=img_wcs.wcs_world2pix(ra, dec, 1)
        y,x = pixels
        self.logger.debug("x:{} y:{}".format(x,y))
        y = int(y)
        x = int(x)
        #Search 50 pixels around, see if in nan is in there
        row_idx = np.array([range(x-num_pixels, x+num_pixels+1)])
        col_idx = np.array([range(y-num_pixels, y+num_pixels+1)])
        data_selection=img_data[0,0,row_idx[:, None], col_idx]
        data_selection = np.nan_to_num(data_selection)
        angle=angle
        bscale=bscale
        data_selection=data_selection.squeeze()
        data_selection=data_selection*bscale
        med, std, mask = self._Median_clip(data_selection, full_output=True, ftol=0.0, max_iter=max_iter, sigma=sigma)
        
        if std == 0:
            self.logger.warning("Local RMS returned 0 or null. Setting to 1.0 mJy.")
            std = 1.0e-3
        return std
    
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
            askap_local_rms_val=self._get_local_rms(row[ra_col], row[dec_col], askap_img_wcs, askap_img_data, askap_img_header["crval1"], askap_img_header["bscale"])
            if using_pre_conv:
                preconv_askap_local_rms_val=self._get_local_rms(row[ra_col], row[dec_col], preconv_askap_img_wcs, preconv_askap_img_data, preconv_askap_img_header["crval1"], preconv_askap_img_header["bscale"])
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
                loaded_fits[fits_image]=askapimage(row["catalog_Mosaic_path"])
                loaded_fits[fits_image].load_wcs()
                loaded_fits[fits_image].load_fits_data()
            catalog_peak_flux = self._get_peak_flux(row[ra_col], row[dec_col], loaded_fits[fits_image].wcs, loaded_fits[fits_image].data, num_pixels=5, sumss=sumss_image)
            peak_fluxes.append(catalog_peak_flux)
        no_matches["measured_catalog_peak_flux"]=peak_fluxes
        
        return no_matches
    
    def inject_transients_db(self, image_id, db_engine="postgresql", 
            db_username="postgres", db_host="localhost", db_port="5432", db_database="postgres", dualmode=False, basecat="sumss"):
        engine = sqlalchemy.create_engine('{}://{}@{}:{}/{}'.format(db_engine, db_username, db_host, db_port, db_database))
        
        #This is just a huge mess - if anyone happens to be here reading this and trying to understand it I can only apologise.
        
        #First do sumss no match
        # match_id=1
        # result = engine.execute("SELECT match_id FROM images_sumssnomatch")
        # try:
        #     match_id+=int(result.fetchall()[-1][-1])
        # except:
        #     match_id=1

        #Database columns now standardised across good matches and transients, the full list of columns is below
        
        full_list=["image_id", "master_name", "sumss_name", "nvss_name", 
                    "master_ra", "master_dec","sumss_St", "sumss_e_St", "nvss_S1.4", "nvss_e_S1.4", 
                    "master_catflux", "master_catalog_Mosaic_rms", "askap_int_flux", "askap_err_int_flux", 
                    "master_scaled_askap_iflux", "askap_non_conv_int_flux", "master_scaled_non_conv_askap_iflux", "master_askap_non_conv_d2d", "measured_askap_peak_flux", "measured_preconv_askap_peak_flux", 
                    "measured_askap_local_rms", "measured_preconv_askap_local_rms", "measured_catalog_peak_flux","askap_cat_int_flux_ratio", "master_scaled_catalog_iflux","askap_rms", 
                    "askap_rms_preconv", "sumss_sumss_snr", "nvss_nvss_snr", "master_cat_snr", "master_cat_to_askap_scaled_snr", "askap_askap_snr",
                    "master_askap_to_cat_scaled_snr", "postage_stamp", "pipelinetag", "usertag", "userreason", "checkedby", "survey_used", "nearest_sources"]
            
        db_col_names = ["image_id", #need to add
                        # "match_id", #need to add
                        "master_name", #done
                        "sumss_name", #done
                        "nvss_name", #done
                        "ra", #done
                        "dec", #done
                        "sumss_iflux", #done
                        "sumss_iflux_e", #done
                        "nvss_iflux", #done
                        "nvss_iflux_e", #done
                        "catalog_iflux", #done
                        "catalog_rms", #done
                        "askap_iflux", #need to add
                        "askap_iflux_e", #need to add
                        "askap_scale_flux", #need to add
                        "askap_non_conv_flux",
                        "askap_non_conv_scaled_flux",
                        "askap_non_conv_d2d",
                        "measured_askap_flux", #done
                        "measured_askap_flux_2", #done
                        "measured_askap_local_rms",
                        "measured_askap_local_rms_2",
                        "measured_catalog_flux", #need to add
                        "askap_cat_ratio",  #need to add 
                        "catalog_scale_flux", #done
                        "askap_rms",    #done
                        "askap_rms_2",    #done
                        "sumss_snr",    #done
                        "nvss_snr",    #done
                        "cat_snr",    #done
                        "scaled_cat_snr",    #done
                        "askap_snr",    #need to add
                        "scaled_askap_snr",    #need to add
                        "ploturl",  #done
                        "pipelinetag",  #done
                        "usertag", #need to add
                        "userreason", #need to add
                        "checkedby", #need to add
                        "survey",
                        "nearest_sources"] #done

        if dualmode:
            filter_list=["master_name", "sumss_name", "nvss_name", 
                    "master_ra", "master_dec","sumss_St", "sumss_e_St", "nvss_S1.4", "nvss_e_S1.4", 
                    "master_catflux","master_catalog_Mosaic_rms", "measured_askap_peak_flux", "measured_preconv_askap_peak_flux", 
                    "measured_askap_local_rms", "measured_preconv_askap_local_rms", "master_scaled_catalog_iflux","askap_rms", 
                    "askap_rms_preconv", "sumss_sumss_snr", "nvss_nvss_snr", "master_cat_snr", "master_cat_to_askap_scaled_snr",
                    "postage_stamp", "pipelinetag", "survey_used", "nearest_sources"]
                    
            sortby = "master_cat_snr"
        elif basecat == "sumss":
            filter_list=["master_name", "sumss_name", 
                    "master_ra", "master_dec","sumss_St", "sumss_e_St", 
                    "master_catflux","master_catalog_Mosaic_rms", "measured_askap_peak_flux", "measured_preconv_askap_peak_flux", 
                    "measured_askap_local_rms", "measured_preconv_askap_local_rms", "master_scaled_catalog_iflux","askap_rms", 
                    "askap_rms_preconv", "sumss_sumss_snr", "master_cat_snr", "master_cat_to_askap_scaled_snr",
                    "postage_stamp", "pipelinetag", "survey_used", "nearest_sources"]
            sortby="sumss_sumss_snr"
        else:
            filter_list=["master_name", "nvss_name", 
                    "master_ra", "master_dec", "nvss_S1.4", "nvss_e_S1.4", 
                    "master_catflux","master_catalog_Mosaic_rms", "measured_askap_peak_flux", "measured_preconv_askap_peak_flux", 
                    "measured_askap_local_rms", "measured_preconv_askap_local_rms","master_scaled_catalog_iflux","askap_rms", 
                    "askap_rms_preconv", "nvss_nvss_snr", "master_cat_snr", "master_cat_to_askap_scaled_snr",
                    "postage_stamp", "pipelinetag", "survey_used", "nearest_sources"]
            sortby="nvss_nvss_snr"
        db_df=self.transients_no_matches_df.filter(filter_list, axis=1).sort_values(by=["pipelinetag", sortby], ascending=[False, False])
        db_df["image_id"]=image_id
        # db_df["match_id"]=[i for i in range(match_id, match_id+len(db_df.index))]
        # db_df=db_df[["image_id",]+filter_list]
        db_df["postage_stamp"]=[os.path.join("media/{}/stamps/{}".format(image_id, i)) for i in db_df["postage_stamp"].values]
        # db_df["pipelinetag"]="N/A"
        db_df["usertag"]="N/A"
        db_df["userreason"]="N/A"
        db_df["checkedby"]="N/A"
        #As this is No ASKAP match to Catalog, we add the following null values
        db_df["askap_int_flux"] = 0.0
        db_df["askap_err_int_flux"] = 0.0
        db_df["askap_non_conv_int_flux"] = 0.0
        db_df["master_scaled_non_conv_askap_iflux"]=0.0
        db_df["master_askap_non_conv_d2d"]=0.0
        db_df["master_scaled_askap_iflux"] = 0.0
        db_df["measured_catalog_peak_flux"]= 0.0 #need to add
        db_df["askap_cat_int_flux_ratio"] = 0.0  #need to add
        db_df["askap_askap_snr"] = 0.0
        db_df["master_askap_to_cat_scaled_snr"] = 0.0
        if not dualmode:
            if basecat=="sumss":
                db_df["nvss_name"]="N/A"
                db_df["nvss_S1.4"]=0.0
                db_df["nvss_e_S1.4"]=0.0
                db_df["nvss_nvss_snr"]=0.0
            else:
                db_df["sumss_name"]="N/A"
                db_df["sumss_St"]=0.0
                db_df["sumss_e_St"]=0.0
                db_df["sumss_sumss_snr"]=0.0
        db_df=db_df[full_list]
        db_df.columns=db_col_names
        values = {'sumss_iflux': 0, 'sumss_iflux_e': 0, 'sumss_name': "N/A", 'sumss_snr': 0, 'nvss_iflux': 0, 'nvss_iflux_e': 0, 'nvss_name': "N/A", 'nvss_snr': 0,
            "nvss_askap_snr":0, "askap_nvss_snr":0, "sumss_askap_snr":0, "askap_sumss_snr":0, "measured_askap_flux":0.0, "measured_askap_flux_2":0.0, 'askap_non_conv_flux':0.0,
            "askap_non_conv_scaled_flux":0.0, "askap_non_conv_d2d":0.0}
        db_df.fillna(value=values, inplace=True)
        db_df.to_sql("images_sumssnomatch", engine, if_exists="append", index=False)
        
        #LargeRatio
        match_id = 1
        # result = engine.execute("SELECT match_id FROM images_largeratio")
        # try:
        #     match_id+=int(result.fetchall()[-1][-1])
        # except:
        #     match_id=1
        
        # if dualmode:
        #     filter_list=["sumss_name", "nvss_name", "sumss__RAJ2000", "sumss__DEJ2000","sumss_St", "sumss_e_St", "nvss_S1.4", "nvss_e_S1.4", "askap_int_flux", "askap_err_int_flux",
        #     "master_scaled_iflux", "sumss_sumss_snr", "nvss_nvss_snr", "postage_stamp", "pipelinetag"]
        # elif basecat == "sumss":
        #     filter_list=["sumss_name", "sumss__RAJ2000", "sumss__DEJ2000","sumss_St", "sumss_e_St", "askap_int_flux", "askap_err_int_flux","master_scaled_iflux",
        #     "sumss_sumss_snr", "master_scaled_iflux", "postage_stamp", "pipelinetag"]
        # else:
        #     filter_list=["nvss_name", "nvss__RAJ2000", "nvss__DEJ2000", "nvss_S1.4", "nvss_e_S1.4", "askap_int_flux", "askap_err_int_flux", "master_scaled_iflux",
        #                  "nvss_nvss_snr", "postage_stamp", "pipelinetag"]  
        
        if dualmode:
            filter_list=["image_id", "master_name", "sumss_name", "nvss_name", 
                            "master_ra", "master_dec","sumss_St", "sumss_e_St", "nvss_S1.4", "nvss_e_S1.4", 
                            "master_catflux","master_catalog_Mosaic_rms", "askap_int_flux", "askap_err_int_flux", 
                            "master_scaled_askap_iflux","askap_non_conv_int_flux", "master_scaled_non_conv_askap_iflux", "master_askap_non_conv_d2d","measured_askap_peak_flux", "measured_preconv_askap_peak_flux", 
                            "measured_askap_local_rms", "measured_preconv_askap_local_rms", "measured_catalog_peak_flux","askap_cat_int_flux_ratio", "master_scaled_catalog_iflux","askap_rms", 
                            "askap_rms_preconv", "sumss_sumss_snr", "nvss_nvss_snr", "master_cat_snr", "master_cat_to_askap_scaled_snr", "askap_askap_snr",
                            "master_askap_to_cat_scaled_snr", "postage_stamp", "pipelinetag", "usertag", "userreason", "checkedby", "survey_used", "nearest_sources"]
        elif basecat == "sumss":
            filter_list=["image_id", "master_name", "sumss_name", 
                    "master_ra", "master_dec","sumss_St", "sumss_e_St", 
                    "master_catflux","master_catalog_Mosaic_rms", "askap_int_flux", "askap_err_int_flux", 
                    "master_scaled_askap_iflux", "askap_non_conv_int_flux", "master_scaled_non_conv_askap_iflux", "master_askap_non_conv_d2d","measured_askap_peak_flux", "measured_preconv_askap_peak_flux", 
                    "measured_askap_local_rms", "measured_preconv_askap_local_rms", "measured_catalog_peak_flux","askap_cat_int_flux_ratio", "master_scaled_catalog_iflux","askap_rms", 
                    "askap_rms_preconv", "sumss_sumss_snr", "master_cat_snr", "master_cat_to_askap_scaled_snr", "askap_askap_snr",
                    "master_askap_to_cat_scaled_snr", "postage_stamp", "pipelinetag", "survey_used", "nearest_sources"]
        else:
            filter_list=["image_id", "master_name", "nvss_name", 
                                "master_ra", "master_dec", "nvss_S1.4", "nvss_e_S1.4", 
                                "master_catflux","master_catalog_Mosaic_rms", "askap_int_flux", "askap_err_int_flux", 
                                "master_scaled_askap_iflux","askap_non_conv_int_flux", "master_scaled_non_conv_askap_iflux", "master_askap_non_conv_d2d","measured_askap_peak_flux", "measured_preconv_askap_peak_flux", 
                                "measured_askap_local_rms", "measured_preconv_askap_local_rms", "measured_catalog_peak_flux","askap_cat_int_flux_ratio", "master_scaled_catalog_iflux","askap_rms", 
                                "askap_rms_preconv", "nvss_nvss_snr", "master_cat_snr", "master_cat_to_askap_scaled_snr", "askap_askap_snr",
                                "master_askap_to_cat_scaled_snr", "postage_stamp", "pipelinetag", "survey_used", "nearest_sources"] 
           
        db_df=self.transients_large_ratios_df.filter(filter_list, axis=1).sort_values(by=["pipelinetag", "askap_cat_int_flux_ratio"], ascending=False)
        db_df["image_id"]=image_id
        db_df["match_id"]=[i for i in range(match_id, match_id+len(db_df.index))]
        db_df["postage_stamp"]=[os.path.join("media/{}/stamps/{}".format(image_id, i)) for i in db_df["postage_stamp"].values]
        # pipeline_tags = ["Convolved match" if "_GOOD_" in i else "Non-convolved match" for i in db_df["postage_stamp"].values]
        # db_df["pipelinetag"]=pipeline_tags
        # db_df=db_df[["image_id", "match_id"]+filter_list]
        #Large ratio won't have any transient measurements
        db_df["measured_askap_peak_flux"]=0.0
        db_df["measured_preconv_askap_peak_flux"]=0.0
        db_df["measured_catalog_peak_flux"]=0.0
        # db_df["measured_askap_local_rms"]=0.0
        # db_df["measured_preconv_askap_local_rms"]=0.0
        db_df["usertag"]="N/A"
        db_df["userreason"]="N/A"
        db_df["checkedby"]="N/A"
        if not dualmode:
            if basecat=="sumss":
                db_df["nvss_name"]="N/A"
                db_df["nvss_S1.4"]=0.0
                db_df["nvss_e_S1.4"]=0.0
                db_df["nvss_nvss_snr"]=0.0
                # db_df["askap_nvss_snr"]=0.0
            else:
                db_df["sumss_name"]="N/A"
                db_df["sumss_St"]=0.0
                db_df["sumss_e_St"]=0.0
                db_df["sumss_sumss_snr"]=0.0
                # db_df["askap_sumss_snr"]=0.0
        db_df=db_df[full_list].sort_values(by=["pipelinetag", "askap_cat_int_flux_ratio"], ascending=[False, False])
        db_df.columns=db_col_names
        values = {'sumss_iflux': 0, 'sumss_iflux_e': 0, 'sumss_name': "N/A", 'sumss_snr': 0, 'nvss_iflux': 0, 'nvss_iflux_e': 0, 'nvss_name': "N/A", 'nvss_snr': 0,
            "nvss_askap_snr":0, "askap_nvss_snr":0, "sumss_askap_snr":0, "askap_sumss_snr":0, "askap_non_conv_flux":0.0, "askap_non_conv_scaled_flux":0.0, "askap_non_conv_d2d":0.0}
        db_df.fillna(value=values, inplace=True)
        db_df.to_sql("images_largeratio", engine, if_exists="append", index=False)
        
        #ASKAP not seen =  This does not use a crossmatch_df so the values are different
        
        filter_list=["name",
                    "ra", "dec", 
                    "catalog_Mosaic_rms", "int_flux", "err_int_flux", 
                    "master_scaled_to_catalog", "non_conv_int_flux", "master_scaled_non_conv_askap_iflux", "non_conv_d2d",
                    "measured_askap_local_rms", "measured_preconv_askap_local_rms", "measured_catalog_peak_flux","rms", 
                    "rms_preconv", "askap_snr",
                    "cat_snr", "postage_stamp", "pipelinetag", "survey_used", "nearest_sources"]
        
        
        # if dualmode:
        #     filter_list=["name", "ra", "dec", "postage_stamp", "int_flux", "err_int_flux", "askap_scaled_to_sumss", "askap_scaled_to_nvss", "askap_snr",
        #     "cat_snr", "cat_peak_flux", "pipelinetag", "survey"]
        # else:
        #     if basecat=="sumss":
        #         filter_list=["name", "ra", "dec", "postage_stamp", "int_flux", "err_int_flux", "askap_scaled_to_sumss", "askap_snr",
        #         "cat_snr", "cat_peak_flux", "pipelinetag", "survey"]
        #     else:
        #         filter_list=["name", "ra", "dec", "postage_stamp", "int_flux", "err_int_flux", "askap_scaled_to_nvss", "askap_snr",
        #         "cat_snr", "cat_peak_flux", "pipelinetag", "survey"]
        
        db_df=self.transients_not_matched_askap_should_see_df.filter(filter_list, axis=1).sort_values(by=["pipelinetag", "cat_snr"], ascending=[False, False])
        db_df.columns=["master_name", "master_ra", "master_dec", "master_catalog_Mosaic_rms", "askap_int_flux", "askap_err_int_flux", 
        "master_scaled_askap_iflux", "askap_non_conv_int_flux", "master_scaled_non_conv_askap_iflux", "master_askap_non_conv_d2d", "measured_askap_local_rms", "measured_preconv_askap_local_rms", "measured_catalog_peak_flux", "askap_rms", "askap_rms_preconv", "askap_askap_snr", "master_askap_to_cat_scaled_snr",
        "postage_stamp", "pipelinetag", "survey_used", "nearest_sources"]
        db_df["image_id"]=image_id
        # db_df["match_id"]=[i for i in range(match_id, match_id+len(db_df.index))]
        #this is askap only so we need to add the respective rows
        db_df["sumss_name"] = "N/A"
        db_df["nvss_name"] = "N/A"
        db_df["sumss_St"] = 0.0
        db_df["sumss_e_St"] = 0.0
        db_df["nvss_S1.4"] = 0.0
        db_df["nvss_e_S1.4"] = 0.0
        db_df["master_catflux"] = 0.0
        db_df["measured_askap_peak_flux"] =0.0
        db_df["measured_preconv_askap_peak_flux"] = 0.0
        # db_df["measured_askap_local_rms"]=0.0
        # db_df["measured_preconv_askap_local_rms"]=0.0
        db_df["askap_cat_int_flux_ratio"] = 0.0
        db_df["master_scaled_catalog_iflux"] = 0.0
        db_df["sumss_sumss_snr"] = 0.0
        db_df["nvss_nvss_snr"] = 0.0
        db_df["master_cat_snr"] = 0.0
        db_df["master_cat_to_askap_scaled_snr"] = 0.0
        # if "cat_peak_flux" not in db_df.columns:
        #     db_df["cat_peak_flux"]=0.0
        # if not dualmode:
        #     if basecat=="sumss":
        #         db_df["askap_scaled_to_nvss"]=0.0
        #         db_df["nvss_snr"]=0.0
        #     else:
        #         db_df["askap_scaled_to_sumss"]=0.0
        #         db_df["sumss_snr"]=0.0
        # db_df=db_df[["image_id",]+filter_list]
        db_df["postage_stamp"]=[os.path.join("media/{}/stamps/{}".format(image_id, i)) for i in db_df["postage_stamp"].values]
        # db_df["pipelinetag"]="N/A"
        db_df["usertag"]="N/A"
        db_df["userreason"]="N/A"
        db_df["checkedby"]="N/A"
        db_df=db_df[full_list]
        db_df.columns=db_col_names
        values = {'askap_scaled_to_sumss': 0, 'askap_scaled_to_nvss': 0, "sumss_snr":0, "nvss_snr":0, "measured_catalog_flux":0.0, "askap_non_conv_flux":0.0,"askap_non_conv_scaled_flux":0.0, "askap_non_conv_d2d":0.0}
        db_df.fillna(value=values, inplace=True)
        db_df.to_sql("images_askapnotseen", engine, if_exists="append", index=False)
        
    def inject_good_db(self, image_id, sumss, nvss, db_engine="postgresql", 
            db_username="postgres", db_host="localhost", db_port="5432", db_database="postgres", max_separation=45.0, dualmode=False, basecat="sumss"):        
        engine = sqlalchemy.create_engine('{}://{}@{}:{}/{}'.format(db_engine, db_username, db_host, db_port, db_database))
        # match_id=1
        # result = engine.execute("SELECT match_id FROM images_goodmatch")
        # try:
        #     match_id+=int(result.fetchall()[-1][-1])
        # except:
        #     match_id=1
        
        #create the good images url
        # if basecat=="sumss":
        #     name_col="sumss_name"
        # else:
        #     name_col="nvss_name"
        
        full_list=["image_id", "master_name", "sumss_name", "nvss_name", 
                    "master_ra", "master_dec","sumss_St", "sumss_e_St", "nvss_S1.4", "nvss_e_S1.4", 
                    "master_catflux","master_catalog_Mosaic_rms", "askap_int_flux", "askap_err_int_flux", 
                    "master_scaled_askap_iflux", "askap_non_conv_int_flux", "master_scaled_non_conv_askap_iflux", "master_askap_non_conv_d2d", "measured_askap_peak_flux", "measured_preconv_askap_peak_flux", 
                    "measured_askap_local_rms","measured_preconv_askap_local_rms",
                    "measured_catalog_peak_flux","askap_cat_int_flux_ratio", "master_scaled_catalog_iflux","askap_rms", 
                    "askap_rms_preconv", "sumss_sumss_snr", "nvss_nvss_snr", "master_cat_snr", "master_cat_to_askap_scaled_snr", "askap_askap_snr",
                    "master_askap_to_cat_scaled_snr", "postage_stamp", "pipelinetag", "usertag", "userreason", "checkedby", "survey_used", "nearest_sources"]
            
        db_col_names = ["image_id", #need to add
                        # "match_id", #need to add
                        "master_name", #done
                        "sumss_name", #done
                        "nvss_name", #done
                        "ra", #done
                        "dec", #done
                        "sumss_iflux", #done
                        "sumss_iflux_e", #done
                        "nvss_iflux", #done
                        "nvss_iflux_e", #done
                        "catalog_iflux", #done
                        "catalog_rms", #done
                        "askap_iflux", #need to add
                        "askap_iflux_e", #need to add
                        "askap_scale_flux", #need to add
                        "askap_non_conv_flux",
                        "askap_non_conv_scaled_flux",
                        "askap_non_conv_d2d",
                        "measured_askap_flux", #done
                        "measured_askap_flux_2", #done
                        "measured_askap_local_rms",
                        "measured_askap_local_rms_2",
                        "measured_catalog_flux", #need to add
                        "askap_cat_ratio",  #need to add 
                        "catalog_scale_flux", #done
                        "askap_rms",    #done
                        "askap_rms_2",    #done
                        "sumss_snr",    #done
                        "nvss_snr",    #done
                        "cat_snr",    #done
                        "scaled_cat_snr",    #done
                        "askap_snr",    #need to add
                        "scaled_askap_snr",    #need to add
                        "ploturl",  #done
                        "pipelinetag",  #done
                        "usertag", #need to add
                        "userreason", #need to add
                        "checkedby", #need to add
                        "survey",
                        "nearest_sources"] #done
        
        postage_stamps=[]
        for i,row in self.goodmatches_df.iterrows():
            if row["d2d"] <= max_separation:
                postage_stamps.append("{0}_{1}_GOOD_sidebyside.jpg".format(basecat.upper(), row["master_name"]))
            else:
                postage_stamps.append("{0}_{1}_BAD_sidebyside.jpg".format(basecat.upper(), row["master_name"]))
        
        self.goodmatches_df["postage_stamp"]=np.array(postage_stamps)
            
        # self.goodmatches_df["postage_stamp"]=["SUMSS_{}_GOOD_sidebyside.jpg".format(i) for i in self.goodmatches_df["sumss_name"].values]
        
        if dualmode:
            filter_list=["image_id", "master_name", "sumss_name", "nvss_name", 
                    "master_ra", "master_dec","sumss_St", "sumss_e_St", "nvss_S1.4", "nvss_e_S1.4", 
                    "master_catflux","master_catalog_Mosaic_rms", "askap_int_flux", "askap_err_int_flux", 
                    "master_scaled_askap_iflux", "askap_non_conv_int_flux","master_scaled_non_conv_askap_iflux", "master_askap_non_conv_d2d",
                    "askap_cat_int_flux_ratio", "master_scaled_catalog_iflux","askap_rms", 
                    "askap_rms_preconv", "sumss_sumss_snr", "nvss_nvss_snr", "master_cat_snr", "master_cat_to_askap_scaled_snr", "askap_askap_snr",
                    "master_askap_to_cat_scaled_snr", "postage_stamp","survey_used", "measured_askap_local_rms", "measured_preconv_askap_local_rms"]
            if sumss:
                sortby="sumss_sumss_snr"
            else:
                sortby="nvss_nvss_snr"
        elif basecat == "sumss":
            filter_list=["image_id", "master_name", "sumss_name",
                    "master_ra", "master_dec","sumss_St", "sumss_e_St",
                    "master_catflux","master_catalog_Mosaic_rms", "askap_int_flux", "askap_err_int_flux", 
                    "master_scaled_askap_iflux", "askap_non_conv_int_flux","master_scaled_non_conv_askap_iflux", "master_askap_non_conv_d2d",
                    "askap_cat_int_flux_ratio", "master_scaled_catalog_iflux","askap_rms", 
                    "askap_rms_preconv", "sumss_sumss_snr", "master_cat_snr", "master_cat_to_askap_scaled_snr", "askap_askap_snr",
                    "master_askap_to_cat_scaled_snr", "postage_stamp","survey_used", "measured_askap_local_rms", "measured_preconv_askap_local_rms"]
            sortby="sumss_sumss_snr"
        else:
            filter_list=["image_id", "master_name", "nvss_name", 
                    "master_ra", "master_dec", "nvss_S1.4", "nvss_e_S1.4", 
                    "master_catflux","master_catalog_Mosaic_rms", "askap_int_flux", "askap_err_int_flux", 
                    "master_scaled_askap_iflux", "askap_non_conv_int_flux","master_scaled_non_conv_askap_iflux", "master_askap_non_conv_d2d",
                    "askap_cat_int_flux_ratio", "master_scaled_catalog_iflux","askap_rms", 
                    "askap_rms_preconv", "nvss_nvss_snr", "master_cat_snr", "master_cat_to_askap_scaled_snr", "askap_askap_snr",
                    "master_askap_to_cat_scaled_snr", "postage_stamp","survey_used", "measured_askap_local_rms", "measured_preconv_askap_local_rms"]
            sortby="nvss_nvss_snr"
        
        db_df=self.goodmatches_df.filter(filter_list, axis=1)
        db_df["image_id"]=image_id
        db_df["postage_stamp"]=[os.path.join("media/{}/stamps/{}".format(image_id, i)) for i in db_df["postage_stamp"].values]
        pipeline_tags = ["Convolved match" if "_GOOD_" in i else "Non-convolved match" for i in db_df["postage_stamp"].values]
        db_df["pipelinetag"]= pipeline_tags
        
        #Good matches so won't have any transient stuff
        db_df["measured_askap_peak_flux"]=0.0
        db_df["measured_preconv_askap_peak_flux"]=0.0
        db_df["measured_catalog_peak_flux"]=0.0
        if "measured_askap_local_rms" not in db_df.columns:
            db_df["measured_askap_local_rms"]=0.0
            db_df["measured_preconv_askap_local_rms"]=0.0
        db_df["nearest_sources"] = "None"
        
        #Add usual other ones
        db_df["usertag"]="N/A"
        db_df["userreason"]="N/A"
        db_df["checkedby"]="N/A"
        if not dualmode:
            if sumss:
                db_df["nvss_name"]="N/A"
                db_df["nvss_S1.4"]=0.0
                db_df["nvss_e_S1.4"]=0.0
                db_df["nvss_nvss_snr"]=0.0
            else:
                db_df["sumss_name"]="N/A"
                db_df["sumss_St"]=0.0
                db_df["sumss_e_St"]=0.0
                db_df["sumss_sumss_snr"]=0.0
        # if sumss:
        #     ra_dec=["sumss__RAJ2000", "sumss__DEJ2000"]
        # else:
        #     ra_dec=["nvss__RAJ2000", "nvss__DEJ2000"]
        db_df=db_df[full_list]
        db_df.columns=db_col_names
        values = {'sumss_iflux': 0, 'sumss_iflux_e': 0, 'sumss_name': "N/A", 'sumss_snr': 0, 'nvss_iflux': 0, 'nvss_iflux_e': 0, 'nvss_name': "N/A", 'nvss_snr': 0, "askap_non_conv_flux":0.0, "askap_non_conv_scaled_flux":0.0, "askap_non_conv_d2d":0.0}
        db_df.fillna(value=values, inplace=True)
        db_df.to_sql("images_goodmatch", engine, if_exists="append", index=False)
        
        
        
