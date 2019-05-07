#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import logging
import aplpy
from astropy.coordinates import SkyCoord
from astropy import units as u
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
    def __init__(self, base_catalog, comp_catalog, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        super(crossmatch, self).__init__()
        self.base_catalog=base_catalog
        self.comp_catalog=comp_catalog
        self.performed=False 
        
    def perform_crossmatch(self, maxsep=15.):
        self.idx, self.d2d, self.d3d = self.base_catalog._crossmatch_catalog.match_to_catalog_sky(self.comp_catalog._crossmatch_catalog)
        self.matches = self.comp_catalog.df.loc[self.idx].reset_index(drop=True)
        newcolnames={}
        for c in self.matches.columns:
            newcolnames[c]="askap_{}".format(c)
        self.matches=self.matches.rename(str, columns=newcolnames)
        self.crossmatch_df = self.base_catalog.df.copy(deep=True)
        newcolnames={}
        for c in self.crossmatch_df.columns:
            newcolnames[c]="sumss_{}".format(c)
        self.crossmatch_df=self.crossmatch_df.rename(str, columns=newcolnames)
        self.crossmatch_df=self.crossmatch_df.join(self.matches)
        self.crossmatch_df["d2d"]=self.d2d.arcsec
        
        self.goodmatches_df = self.crossmatch_df[self.crossmatch_df["d2d"]<=maxsep].reset_index(drop=True)
        self.badmatches_df = self.crossmatch_df[self.crossmatch_df["d2d"]>maxsep].reset_index(drop=True)
        
        self.performed=True
        self.logger.info("Crossmatch complete.")
        
    def write_crossmatch(self, outname):
        if self.performed:
            self.crossmatch_df.to_csv(outname, index=False)
            self.logger.info("Written crossmatch dataframe to {}.".format(outname))
        else:
            self.logger.error("Need to perform a cross match first!")
        
    def calculate_ratio(self, col1, col2, output_col_name, col1_scaling=1., col2_scaling=1.):
        self.crossmatch_df[output_col_name]=(self.crossmatch_df[col1]*col1_scaling)/(self.crossmatch_df[col2]*col2_scaling)
        
    def calculate_diff(self, col1, col2, output_col_name, col1_scaling=1., col2_scaling=1.):
        self.crossmatch_df[output_col_name]=(self.crossmatch_df[col1]*col1_scaling)-(self.crossmatch_df[col2]*col2_scaling)
        
    def _plotinitial(self, panels, key, figure, image):
        panels[key]=aplpy.FITSFigure(image, figure=figure, subplot=(1,2,key+1))
        panels[key].show_grayscale()
            # panels[key].show_contour(images[n-1], colors='red', levels=[3.*12e-6, 4.*12.e-6, 8.*12.e-6, 16*12.e-6])
        panels[key].set_theme('publication')
        return panels
    
    def produce_postage_stamps(self, sumss_mosaic_dir, postage_options, radius=15./60., nprocs=1, max_separation=15., convolve=False, pre_convolve_image=None, preconolve_catalog=None):
        askap_fits = self.crossmatch_df["askap_image"].iloc[0]
        postage_stamps.crossmatch_stamps(self, askap_fits, postage_options, nprocs, sumss_mosaic_dir, radius=13./60., max_separation=max_separation, convolve=convolve,
            askap_pre_convolve_image=pre_convolve_image, askap_pre_convolve_catalog=preconolve_catalog)
        
    # def produce_postage_stamps(self, sumss_mosaic_dir, radius=13./60., max_separation=None):
    #     if max_separation!=None:
    #         postage_df=self.crossmatch_df[self.crossmatch_df["d2d"]<max_separation].reset_index(drop=True)
    #     else:
    #         postage_df=self.crossmatch_df.copy(deep=True)
    #
    #     # self.logger.info("Estimated time to completion = {:.2f} hours".format(len(postage_df.index)*6./3600.))
    #     #Minimise fits opening so first get a list of all the unique SUMSS fits files to be used
    #     sumss_fits_mosaics = postage_df["sumss_Mosaic"].unique()
    #     #For now support one ASKAP image at a time so can just take this from the first entry.
    #     askap_fits = postage_df["askap_image"].iloc[0]
    #
    #     #Create dictionary of SUMSS images and the
    #
    #     #First initialise the figure and set up the askap panel with the image.
    #     panels={}
    #     fig = plt.figure(figsize=(12, 6))
    #     panels=self._plotinitial(panels, 0, fig, askap_fits)
    #
    #     #Load the ellipses onto the askap image
    #     panels[0].show_ellipses(postage_df["askap_ra"],postage_df["askap_dec"],postage_df["askap_b"]/3600.,
    #         postage_df["askap_a"]/3600., angle=postage_df["askap_pa"], layer="ASKAP Sources", color="#1f77b4")
    #     panels[0].show_ellipses(postage_df["sumss__RAJ2000"],postage_df["sumss__DEJ2000"],postage_df["sumss_MinAxis"]/3600.,
    #         postage_df["sumss_MajAxis"]/3600., angle=postage_df["sumss_PA"], layer="SUMSS Sources", color="#d62728")
    #     # panels[0].set_title("ASKAP")
    #
    #     #Now start the SUMSS loop
    #     for s_image in sumss_fits_mosaics:
    #         #Get the actual fits file
    #         self.logger.debug("SUMSS image: {}".format(s_image))
    #         s_image_path=os.path.join(sumss_mosaic_dir, s_image+".FITS")
    #         #Filter the dataframe such that only the SUMSS sources are present
    #         filtered_cross_matches=postage_df[postage_df["sumss_Mosaic"]==s_image].reset_index(drop=True)
    #         #Generate the base SUMSS panel
    #         panels=self._plotinitial(panels, 1, fig, s_image_path)
    #         # panels[1].set_title("SUMSS")
    #         #Add the sources
    #         panels[1].show_ellipses(postage_df["askap_ra"],postage_df["askap_dec"],postage_df["askap_b"]/3600.,
    #             postage_df["askap_a"]/3600., angle=postage_df["askap_pa"], layer="ASKAP Sources", color="#1f77b4", label="ASKAP Sources")
    #         panels[1].show_ellipses(postage_df["sumss__RAJ2000"],postage_df["sumss__DEJ2000"],postage_df["sumss_MinAxis"]/3600.,
    #             postage_df["sumss_MajAxis"]/3600., angle=postage_df["sumss_PA"], layer="SUMSS Sources", color="#d62728", label="SUMSS Sources")
    #         panels[1].axis_labels.hide()
    #         panels[1].tick_labels.hide()
    #         #Now begin the main loop per source
    #         for i, row in filtered_cross_matches.iterrows():
    #             panels[0].set_title("ASKAP "+row["askap_name"])
    #             panels[1].set_title("SUMSS "+row["sumss_name"])
    #             #Centre each image on the ASKAP coordinates for clarity
    #             recentre_ra=row["askap_ra"]
    #             recentre_dec=row["askap_dec"]
    #             for p in panels:
    #                 panels[p].recenter(recentre_ra, recentre_dec, radius)
    #             panels[0].show_circles([recentre_ra], [recentre_dec], 120./3600., color='C1', label="ASKAP source", layer="ASKAP Source")
    #             panels[1].show_circles([row["sumss__RAJ2000"]], [row["sumss__DEJ2000"]],120./3600., color='C9', label="SUMSS source", layer="SUMSS Source")
    #
    #             sep_text=plt.text(0.02, 0.02, "Distance Separation = {:.2f} arcsec".format(row["d2d"]), transform=plt.gcf().transFigure)
    #             ratio_text=plt.text(0.8, 0.02, "Int. Flux Ratio ASKAP/SUMSS = {:.2f}".format(row["askap_sumss_int_flux_ratio"]), transform=plt.gcf().transFigure)
    #
    #             #Figure name
    #             # plt.title(row["sumss_name"])
    #             figname = "SUMSS_{}_sidebyside.png".format(row["sumss_name"])
    #
    #             custom_lines = [Line2D([0], [0], color='#1f77b4'),
    #                             Line2D([0], [0], color='#d62728'),
    #                             Line2D([0], [0], color='C1'),
    #                             Line2D([0], [0], color='C9')]
    #             plt.gca().legend(custom_lines, ["ASKAP Sources", "SUMSSS Sources", "Matched ASKAP", "Matched SUMSS"])
    #
    #             plt.savefig(figname, bbox_inches="tight")
    #
    #             self.logger.info("Saved figure {}.".format(figname))
    #
    #             panels[0].remove_layer("ASKAP Source")
    #             panels[1].remove_layer("SUMSS Source")
    #
    #             sep_text.set_visible(False)
    #             ratio_text.set_visible(False)
    #
    #         #I think this will clear the SUMSS one, must be a more specific way.
    #         plt.gca().remove()
    #
    #     plt.close()
        
        
    def transient_search(self, max_separation=15.0, askap_sumss_snr_thresh=5.0, large_flux_thresh=3.0, pre_conv_crossmatch=None, image_beam_maj=45., image_beam_min=45., image_sumss_beam_maj=45., image_sumss_beam_min=45.):
        # Stage 1 - Define SUMSS sources with no match to ASKAP.
        # Stage 2 - Find those that have large flux ratios.
        # Stage 3 - Find ASKAP sources above the SUMSS threshold that have not been matched.
        
        self.logger.info("Starting transient analysis")
        self.logger.info("Creating matched and non-matched lists")
        # First sort out the matched and non-matched
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
            moved_names = rows["sumss_name"].tolist()
        #Check the pre-convolved catalogue (if available) to see if a match is found there
        if pre_conv_crossmatch !=None:
            self.logger.info("Checking pre-convolved catalogue to check for matches that have been convolved out")
            preconv_df=pre_conv_crossmatch.crossmatch_df
            no_matches_names = no_matches["sumss_name"].tolist()
            preconv_df_no_matches = preconv_df[preconv_df["sumss_name"].isin(no_matches_names)].reset_index(drop=True)
            preconv_df_have_match = preconv_df_no_matches[preconv_df_no_matches["d2d"]<max_separation].reset_index(drop=True)
            if len(preconv_df_have_match.index)>0:
                self.logger.info("{} no match sources have a match in the pre-convolved image.".format(len(preconv_df_have_match.index)))
                preconv_matches_names = preconv_df_have_match["sumss_name"].tolist()
                #Need to move these ones that do have a match to good matches
                rows_to_move = no_matches[no_matches["sumss_name"].isin(preconv_matches_names)]
                matches = matches.append(rows_to_move, ignore_index=True)
                no_matches.drop(rows_to_move.index, inplace=True)
                #We need to know which ones were moved to get the right image - this is a bit of a hacky fix unfortunately
                # conv_moved_names = rows["sumss_name"].tolist()
                #Update the goodmatches for injection later
                self.goodmatches_df=matches
            else:
                self.logger.info("No further matches found in the convolved image.")
                   
                
        # Stage 1
        # Consider all SUMSS sources with matches > max_sep to be 'not matched'
        self.logger.info("Finding SUMSS sources with no ASKAP matches...")
        #For moved double matched sources the tag needs to be changed from BAD -> GOOD to load the correct image.
        if len(moved_names)==0:
            no_match_postage_stamps=["SUMSS_{}_BAD_sidebyside.jpg".format(val) for i,val in no_matches["sumss_name"].iteritems()]
        else:
            no_match_postage_stamps=[]
            for i,val in no_matches["sumss_name"].iteritems():
                if val not in moved_names:
                    no_match_postage_stamps.append("SUMSS_{}_BAD_sidebyside.jpg".format(val))
                else:
                    no_match_postage_stamps.append("SUMSS_{}_GOOD_sidebyside.jpg".format(val))
    
        no_matches["postage_stamp"]=no_match_postage_stamps
        
        #Now check the sources to assign guesstimate type
        #First check bright source proximity
        try:
            self.logger.info("Loading SUMSS bright sources")
            sumss_bright_sources_df=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/sumss_bright_sources.csv"), delimiter=";")
            sumss_bright_tomatch = SkyCoord(ra=sumss_bright_sources_df["_RAJ2000"].values*u.deg, dec=sumss_bright_sources_df["_DEJ2000"].values*u.deg)
        
            pipeline_tags=[]
        
            for i, row in no_matches.iterrows():
                sumss_target=SkyCoord(ra=row["sumss__RAJ2000"]*u.deg, dec=row["sumss__DEJ2000"]*u.deg)
                seps = sumss_target.separation(sumss_bright_tomatch)
                if np.min(seps.arcminute) <= 10.:
                    pipeline_tags.append("Likely artefact (bright source)")
                    #If it passes the bright source, then check if it is extremely extended
                elif (row["sumss_MajAxis"] > (1.75 * image_sumss_beam_maj)) or (row["sumss_MinAxis"] > (1.75 * image_sumss_beam_min)):
                    pipeline_tags.append("Likely artefact (extended)")
                else:
                    pipeline_tags.append("Candidate")
                # print d2d_sumss.deg[min_index]
                # print min_index
            no_matches["pipelinetag"]=pipeline_tags
        except:
            self.logger.error("SUMSS bright source data cannot be found!")
            self.logger.info("Skipping bright source check.")
            no_matches["pipelinetag"]="N/A"
        
        self.transients_no_matches_df=no_matches
        no_matches.to_csv("transients_sumss_sources_no_askap_match.csv", index=False)
        self.logger.info("Written 'transients_sumss_sources_no_askap_match.csv'.")
        
        # Stage 2
        # Here we want actual matches but with 'large' flux ratios
        self.logger.info("Finding matches with large integrated flux ratios...")
        median_match_flux_ratio=matches["askap_sumss_int_flux_ratio"].median()
        std_match_flux_ratio=matches["askap_sumss_int_flux_ratio"].std()
        #For now define sources with 'large' ratios as being more than median +/- std.
        large_ratios=matches[(matches["askap_sumss_int_flux_ratio"]<(median_match_flux_ratio-(large_flux_thresh*std_match_flux_ratio))) | 
            (matches["askap_sumss_int_flux_ratio"]>(median_match_flux_ratio+(large_flux_thresh*std_match_flux_ratio)))].reset_index(drop=True)
        large_ratios_postage_stamps=[]
        for i,row in large_ratios.iterrows():
            if row["d2d"] <= max_separation:
                large_ratios_postage_stamps.append("SUMSS_{}_GOOD_sidebyside.jpg".format(row["sumss_name"]))
            else:
                large_ratios_postage_stamps.append("SUMSS_{}_BAD_sidebyside.jpg".format(row["sumss_name"]))
        large_ratios["postage_stamp"] = large_ratios_postage_stamps
        self.transients_large_ratios_df=large_ratios
        large_ratios.to_csv("transients_sumss_sources_matched_large_flux_ratio.csv", index=False)
        self.logger.info("Written 'transients_sumss_sources_matched_large_flux_ratio.csv'.")
        
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
        not_matched_askap_sources_should_see=not_matched_askap_sources[not_matched_askap_sources["sumss_snr"]>=askap_sumss_snr_thresh].reset_index(drop=True)
        # find_sumss_image_and_flux(not_matched_askap_sources_should_see)
        askapnotseen_postage_stamps=["transient_askapnotseen_ASKAP_{}_sidebyside.jpg".format(val) for i,val in not_matched_askap_sources_should_see["name"].iteritems()]
        not_matched_askap_sources_should_see["postage_stamp"]=askapnotseen_postage_stamps
        
        #Perform checks for these 'transients'
        #1. Doubles
        #2. Extended/Diffuse
        #3. Edge - need to work out how to do
        pipeline_tags=[]
        #Get the askap catalogue to match to
        askap_sources_tomatch = self.comp_catalog._crossmatch_catalog
        for i,row in not_matched_askap_sources_should_see.iterrows():
            askap_target=SkyCoord(ra=row["ra"]*u.deg, dec=row["dec"]*u.deg)
            seps = askap_target.separation(askap_sources_tomatch)
            if sorted(seps.arcsecond)[1] <= 3.*45.:
                pipeline_tags.append("Likely double")
            elif (row["a"] > (1.5 * image_beam_maj)) or (row["b"] > (1.5 * image_beam_min)):
                pipeline_tags.append("Likely diffuse/extended")
            else:
                pipeline_tags.append("Candidate")
            
        not_matched_askap_sources_should_see["pipelinetag"]=pipeline_tags
        
        self.transients_not_matched_askap_should_see_df=not_matched_askap_sources_should_see
        not_matched_askap_sources_should_see.to_csv("transients_askap_sources_no_match_should_be_seen.csv", index=False)
        self.logger.info("Written 'transients_askap_sources_no_match_should_be_seen.csv'.")
        
    def inject_transients_db(self, image_id, db_engine="postgresql", 
            db_username="postgres", db_host="localhost", db_port="5432", db_database="postgres"):
        engine = sqlalchemy.create_engine('{}://{}@{}:{}/{}'.format(db_engine, db_username, db_host, db_port, db_database))
        
        #First do sumss no match
        match_id=1
        result = engine.execute("SELECT match_id FROM images_sumssnomatch")
        try:
            match_id+=int(result.fetchall()[-1][-1])
        except:
            match_id=1
        #Better to control the order
        # plots_columns=["flux_ratio_image_view", "position_offset", "source_counts", "flux_ratios", "flux_ratios_from_centre"]
        # self.logger.info("Image run assigned id {}".format(image_id))
        # conn = psycopg2.connect("host=localhost dbname=RACS user=aste7152 port=5434")
        db_df=self.transients_no_matches_df.filter(["sumss_name", "sumss__RAJ2000", "sumss__DEJ2000","sumss_St", "sumss_e_St", "askap_int_flux", "askap_err_int_flux", 
            "sumss_sumss_snr", "postage_stamp", "pipelinetag"], axis=1).sort_values(by=["sumss_sumss_snr"], ascending=False)
        db_df["image_id"]=image_id
        db_df["match_id"]=[i for i in range(match_id, match_id+len(db_df.index))]
        db_df=db_df[["image_id", "match_id", "sumss_name", "sumss__RAJ2000", "sumss__DEJ2000", "sumss_St", "sumss_e_St", "askap_int_flux", "askap_err_int_flux", "sumss_sumss_snr", "postage_stamp", "pipelinetag"]]
        db_df["postage_stamp"]=[os.path.join("media/{}/stamps/{}".format(image_id, i)) for i in db_df["postage_stamp"].values]
        # db_df["pipelinetag"]="N/A"
        db_df["usertag"]="N/A"
        db_df["userreason"]="N/A"
        db_df["checkedby"]="N/A"
        db_df.columns=["image_id", "match_id", "sumss_name", "ra", "dec", "sumss_iflux", "sumss_iflux_e", "askap_iflux", "askap_iflux_e", "sumss_snr","ploturl", "pipelinetag", "usertag", "userreason", "checkedby"]
        # newpostage_sta
        
        # tempdf=pd.DataFrame([[image_id, self.imagename, "ASKAP RACS image", self.centre.ra.degree,
            # self.centre.dec.degree, datestamp, self.image]+plots_values], columns=["image_id", "name", "description", "ra", "dec", "runtime", "url"]+plots_columns)
        # cols=db_df.columns.tolist()
        db_df.to_sql("images_sumssnomatch", engine, if_exists="append", index=False)
        
        #LargeRatio
        match_id = 1
        result = engine.execute("SELECT match_id FROM images_largeratio")
        try:
            match_id+=int(result.fetchall()[-1][-1])
        except:
            match_id=1
        db_df=self.transients_large_ratios_df.filter(["sumss_name", "sumss__RAJ2000", "sumss__DEJ2000","sumss_St", "sumss_e_St", "askap_int_flux", "askap_err_int_flux", "sumss_sumss_snr",
            "askap_sumss_int_flux_ratio", "postage_stamp"], axis=1).sort_values(by=["askap_sumss_int_flux_ratio"], ascending=False)
        db_df["image_id"]=image_id
        db_df["match_id"]=[i for i in range(match_id, match_id+len(db_df.index))]
        db_df=db_df[["image_id", "match_id", "sumss_name", "sumss__RAJ2000", "sumss__DEJ2000", "sumss_St", "sumss_e_St", "askap_int_flux", "askap_err_int_flux", "sumss_sumss_snr", "askap_sumss_int_flux_ratio", "postage_stamp"]]
        db_df["postage_stamp"]=[os.path.join("media/{}/stamps/{}".format(image_id, i)) for i in db_df["postage_stamp"].values]
        pipeline_tags = ["Convolved match" if "_GOOD_" in i else "Non-convolved match" for i in db_df["postage_stamp"].values]
        db_df["pipelinetag"]=pipeline_tags
        db_df["usertag"]="N/A"
        db_df["userreason"]="N/A"
        db_df["checkedby"]="N/A"
        db_df.columns=["image_id", "match_id", "sumss_name", "ra", "dec", "sumss_iflux", "sumss_iflux_e", "askap_iflux", "askap_iflux_e", "sumss_snr","askap_sumss_ratio", "ploturl", "pipelinetag", "usertag","userreason", "checkedby"]
        db_df.to_sql("images_largeratio", engine, if_exists="append", index=False)
        
        #ASKAP not seen
        match_id=1
        result = engine.execute("SELECT match_id FROM images_askapnotseen")
        try:
            match_id+=int(result.fetchall()[-1][-1])
        except:
            match_id=1
        
        db_df=self.transients_not_matched_askap_should_see_df.filter(["name", "ra", "dec", "postage_stamp", "int_flux", "err_int_flux", 
            "sumss_snr", "sumss_peak_flux", "pipelinetag"], axis=1).sort_values(by=["sumss_snr"], ascending=False)
        db_df["image_id"]=image_id
        db_df["match_id"]=[i for i in range(match_id, match_id+len(db_df.index))]
        if "sumss_peak_flux" not in db_df.columns:
            db_df["sumss_peak_flux"]=0.0
        db_df=db_df[["image_id", "match_id", "name", "ra", "dec", "int_flux", "err_int_flux", "sumss_snr", "sumss_peak_flux", "postage_stamp", "pipelinetag"]]
        db_df["postage_stamp"]=[os.path.join("media/{}/stamps/{}".format(image_id, i)) for i in db_df["postage_stamp"].values]
        # db_df["pipelinetag"]="N/A"
        db_df["usertag"]="N/A"
        db_df["userreason"]="N/A"
        db_df["checkedby"]="N/A"
        db_df.columns=["image_id", "match_id", "askap_name", "ra", "dec", "askap_iflux", "askap_iflux_e", "sumss_snr", "sumss_flux", "ploturl", "pipelinetag", "usertag", "userreason", "checkedby"]
        db_df.to_sql("images_askapnotseen", engine, if_exists="append", index=False)
        
    def inject_good_db(self, image_id, db_engine="postgresql", 
            db_username="postgres", db_host="localhost", db_port="5432", db_database="postgres", max_separation=45.0):        
        engine = sqlalchemy.create_engine('{}://{}@{}:{}/{}'.format(db_engine, db_username, db_host, db_port, db_database))
        match_id=1
        result = engine.execute("SELECT match_id FROM images_goodmatch")
        try:
            match_id+=int(result.fetchall()[-1][-1])
        except:
            match_id=1
        
        #create the good images url
        postage_stamps=[]
        for i,row in self.goodmatches_df.iterrows():
            if row["d2d"] <= max_separation:
                postage_stamps.append("SUMSS_{}_GOOD_sidebyside.jpg".format(row["sumss_name"]))
            else:
                postage_stamps.append("SUMSS_{}_BAD_sidebyside.jpg".format(row["sumss_name"]))
        
        self.goodmatches_df["postage_stamp"]=np.array(postage_stamps)
            
        # self.goodmatches_df["postage_stamp"]=["SUMSS_{}_GOOD_sidebyside.jpg".format(i) for i in self.goodmatches_df["sumss_name"].values]
        
        db_df=self.goodmatches_df.filter(["sumss_name", "sumss__RAJ2000", "sumss__DEJ2000", "sumss_St", "sumss_e_St", "askap_int_flux", "askap_err_int_flux", "sumss_sumss_snr", "postage_stamp"], axis=1)
        db_df["image_id"]=image_id
        db_df["match_id"]=[i for i in range(match_id, match_id+len(db_df.index))]
        db_df=db_df[["image_id", "match_id", "sumss_name", "sumss__RAJ2000", "sumss__DEJ2000", "sumss_St", "sumss_e_St", "askap_int_flux", "askap_err_int_flux", "sumss_sumss_snr", "postage_stamp"]]
        db_df["postage_stamp"]=[os.path.join("media/{}/stamps/{}".format(image_id, i)) for i in db_df["postage_stamp"].values]
        pipeline_tags = ["Convolved match" if "_GOOD_" in i else "Non-convolved match" for i in db_df["postage_stamp"].values]
        db_df["pipelinetag"]= pipeline_tags
        db_df["usertag"]="N/A"
        db_df["userreason"]="N/A"
        db_df["checkedby"]="N/A"
        db_df.columns=["image_id", "match_id", "sumss_name", "ra", "dec", "sumss_iflux", "sumss_iflux_e", "askap_iflux", "askap_iflux_e", "sumss_snr","ploturl", "pipelinetag", "usertag", "userreason", "checkedby"]
        db_df.to_sql("images_goodmatch", engine, if_exists="append", index=False)
        
        
        
