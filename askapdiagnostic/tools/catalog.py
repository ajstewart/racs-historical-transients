#!/usr/bin/env python

import logging
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import pkg_resources
import pandas as pd
import os
from askapdiagnostic.tools.fitsimage import askapimage

class Catalog(object):
    """docstring for survey"""
    def __init__(self, df, survey_name, ref_name, ra_col="_RAJ2000", dec_col="_DEJ2000", flux_col="St", frequency=846e6, add_name_col=False, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self.df = df
        self.survey_name=survey_name
        self.ref_name=ref_name
        self.ra_col=ra_col
        self.dec_col=dec_col
        self.flux_col=flux_col
        self.frequency=frequency
        self._gen_catalog_for_crossmatch()
        if add_name_col:
            self._add_name_col()
        super(Catalog, self).__init__()

        
    def _gen_catalog_for_crossmatch(self):
        cat_ra=self.df[self.ra_col].values
        cat_dec=self.df[self.dec_col].values
        
        self._crossmatch_catalog = SkyCoord(ra=cat_ra*u.degree, dec=cat_dec*u.degree)
        
    def _add_name_col(self):
        self.df["name"]=self._crossmatch_catalog.to_string('hmsdms')
        self.df["name"]=self.df["name"].str.replace(" ", "_")
        
    def add_telescope_beam_columns(self, mode, manual_bmaj=45.0, manual_bmin=45.0):
        allowed_modes=["sumss", "nvss", "manual"]
        if mode not in allowed_modes:
            self.logger.error("Telescope beam mode not recongised.")
            return
        if mode=="sumss":
            self.df["telescope_bmaj"]=self._calculate_sumss_beam_bmaj()
            self.df["telescope_bmin"]=45.0
        elif mode=="nvss":
            self.df["telescope_bmaj"]=45.0
            self.df["telescope_bmin"]=45.0
        else:
            self.df["telescope_bmaj"]=manual_bmaj
            self.df["telescope_bmin"]=manual_bmin
     
    def _sumss_rms(self,dec):
        if dec>-50:
            return 0.002
        else:
            return 0.0012
    
    def _calculate_sumss_beam_bmaj(self):
        img_sumss_bmaj = 45.*(1./np.sin(np.deg2rad(np.abs(self.df[self.dec_col]))))
        return img_sumss_bmaj
    
    def add_single_val_col(self, colname, value, clobber=False):
        if colname in self.df.columns:
            if not clobber:
                self.logger.error("Column {} already present in Catalog and clobber = False, will not overwrite.".format(colname))
                return
            else:
                self.logger.warning("Overwriting column {} with new value {}.".format(colname, value))
        self.df[colname] = value
        
    def get_col(self, columname):
        return self.df[columname].values
        
    def remove_extended(self, threshold=1.2, ellipse_a="MajAxis", ellipse_b="MinAxis", beam_a=45., beam_b=45., ellipse_unit="arcsec", sumss_psf=False, nvss_psf=False):
        if sumss_psf:
            sumss_bmaj = 45.*(1./np.sin(np.deg2rad(np.abs(self.df[self.dec_col]))))
            sumss_bmin = 45.
            self.sumss_no_ext_cat = self.df[(self.df[ellipse_a] <= threshold*sumss_bmaj) & 
                (self.df[ellipse_b] <= threshold*sumss_bmin)].reset_index(drop=True)
            self.sumss_ext_cat = self.df[(self.df[ellipse_a] > threshold*sumss_bmaj) | 
                (self.df[ellipse_b] > threshold*sumss_bmin)].reset_index(drop=True)
            return self.sumss_no_ext_cat
        elif nvss_psf:
            nvss_bmaj = 45.
            nvss_bmin = 45.
            self.nvss_no_ext_cat = self.df[(self.df[ellipse_a] <= threshold*nvss_bmaj) & 
                (self.df[ellipse_b] <= threshold*nvss_bmin)].reset_index(drop=True)
            self.nvss_ext_cat = self.df[(self.df[ellipse_a] > threshold*nvss_bmaj) | 
                (self.df[ellipse_b] > threshold*nvss_bmin)].reset_index(drop=True)
            return self.nvss_no_ext_cat
        else:
            # askap_beam = 45.*45.*(1./np.sin(np.deg2rad(np.abs(self.df[self.dec_col]))))
            self.askap_no_ext_cat=self.df[(self.df[ellipse_a] <= threshold*beam_a) & (self.df[ellipse_b] <= threshold*beam_b)].reset_index(drop=True)
            self.askap_ext_cat=self.df[(self.df[ellipse_a] > threshold*beam_a) | (self.df[ellipse_b] > threshold*beam_b)].reset_index(drop=True)
            return self.askap_no_ext_cat
            
    def askap_remove_out_of_sumss_bounds(self, max_dec):
        #add two SUMSS beam width as a buffer
        max_dec+=2.*(45./3600.)
        self.logger.info("Removing ASKAP sources beyond the SUMSS border.")
        before = len(self.df.index)
        self.df = self.df[self.df[self.dec_col] < max_dec].reset_index=True
        after = len(self.df.index)
        self.logger("{} sources removed leaving {}.".format(before-after, after))
        
    def add_distance_from_pos(self, position, label="centre"):
        self.df["distance_from_{}".format(label)]=position.separation(self._crossmatch_catalog).deg
     
    def write_ann(self, name="", color="GREEN", ellipse_a="MajAxis", ellipse_b="MinAxis", ellipse_pa="PA", ellipse_unit="arcsec"):
        conversions={"arcsec":3600., "arcmin":60., "deg":1.}
        if name=="":
            name=self.survey_name+".ann"
        catalog = SkyCoord(ra=self.df[self.ra_col], dec=self.df[self.dec_col], unit=(u.deg, u.deg))
        
        with open(name, 'w') as f:
            f.write("COORD W\n")
            f.write("PA STANDARD\n")
            f.write("COLOR {}\n".format(color))
            f.write("FONT hershey14\n")
            for i,val in enumerate(catalog):
                f.write("ELLIPSE {} {} {} {} {}\n".format(val.ra.deg, val.dec.deg, 
                self.df[ellipse_a].iloc[i]/conversions[ellipse_unit], self.df[ellipse_b].iloc[i]/conversions[ellipse_unit], self.df[ellipse_pa].iloc[i]))
                
        self.logger.info("Wrote annotation file {}.".format(name))
        
    def add_sumss_sn(self, flux_col="int_flux", dec_col="dec", flux_scaling=1., use_image_rms=False):
        if use_image_rms:
            if self.survey_name.lower()=="sumss":
                self.df["sumss_snr"]=(flux_scaling*self.df[flux_col])/self.df["Mosaic_rms"]
            else:
                #Dealing with an ASKAP catalog, some may have been matched to NVSS images
                snrs=[]
                for i,row in self.df.iterrows():
                    if row["catalog_Mosaic"].startswith("J"):
                        snrs.append((flux_scaling*row[flux_col])/row["catalog_Mosaic_rms"])
                    else:
                        self.logger.warning("{} was matched to a non-SUMSS image, falling back to Dec based RMS.".format(row["name"]))
                        snrs.append((flux_scaling*row[flux_col])/self._sumss_rms(row[dec_col]))
                self.df["sumss_snr"]=snrs
        else:
            self.df["sumss_snr"]=(flux_scaling*self.df[flux_col])/self.df[dec_col].apply(self._sumss_rms)
        
    def add_nvss_sn(self, flux_col="S1.4", dec_col="_DEJ2000", flux_scaling=1., use_image_rms=False):
        if use_image_rms:
            if self.survey_name.lower()=="nvss":
                self.df["nvss_snr"]=(flux_scaling*self.df[flux_col])/self.df["Mosaic_rms"]
            else:
                #Dealing with an ASKAP catalog, some may have been matched to SUMSS image
                snrs=[]
                for i,row in self.df.iterrows():
                    if row["catalog_Mosaic"].startswith("J"):
                        self.logger.warning("{} was matched to a non-NVSS image, falling back to standard RMS.".format(row["name"]))
                        snrs.append((flux_scaling*row[flux_col])/0.0005)
                    else:
                        snrs.append((flux_scaling*row[flux_col])/row["catalog_Mosaic_rms"])
                self.df["nvss_snr"]=snrs
        else:
            self.df["nvss_snr"]=(flux_scaling*self.df[flux_col])/0.0005
        
    def _add_askap_sn(self):
        try:
            self.df["askap_snr"]=self.df["int_flux"]/self.df["local_rms"]
        except:
            self.logger.error("Adding ASKAP SNR not supported for this dataframe.")
            
    def add_manual_sn(self, rms, sn_col_name, flux_col="St", flux_scaling=1.):
        self.df[sn_col_name]=(flux_scaling*self.df[flux_col])/rms
    
    def _scale_flux(self, freq1, freq2, flux_col, si=-0.8):
        return ((freq1/freq2)**(si))*self.df[flux_col]
        
    def calculate_scaled_flux(self, label, frequency=-1, to_catalog="none", flux_col="int_flux", si=-0.8):
        to_catalog=to_catalog.lower()
        frequencies={
            "sumss":843.e6,
            'nvss':1400.e6
        }
        
        if to_catalog!="none":
            if to_catalog not in frequencies:
                self.logger.error("Catalog '{}' not recongised for flux scaling".format(to_catalog))
                return
            else:
                this_freq=frequencies[to_catalog]
                # label=to_catalog
        else:
            this_freq=frequency
            label=label
        
            if this_freq==-1:
                self.logger.error("No catalog or frequency defined for flux scaling.")
                return
                
        
        self.df["{}_scaled_to_{}".format(self.survey_name, label)]=self._scale_flux(this_freq, self.frequency, flux_col, si=si)
        self.logger.info("Fluxes scaled to {} catalogue frequency ({} MHz)".format(label, this_freq/1.e6))
        
    def _add_sumss_mosaic_info(self, sumss_mosaic_dir):
        if self.survey_name.lower()!="sumss":
            self.logger.error("Catalog is not defined as a SUMSS catalog. Cannot add SUMSS mosaic directory information.")
            return
        try:
            self.logger.info("Loading SUMSS image data.")
            sumss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/sumss_images_info.csv"))
            sumss_centres = SkyCoord(ra=sumss_mosaic_data["center-ra"].values*u.deg, dec=sumss_mosaic_data["center-dec"].values*u.deg)
        except:
            self.logger.error("SUMSS mosaic data cannot be found!")
        full_paths = [os.path.join(sumss_mosaic_dir, i+".FITS") for i in self.df["Mosaic"].tolist()]
        rms_values = [sumss_mosaic_data[sumss_mosaic_data["image"]==i+".FITS"].iloc[0]["rms"] for i in self.df["Mosaic"].tolist()]
        self.logger.debug("RMS values {}".format(rms_values))
        self.logger.info("Adding SUMSS mosaic full path information.")
        self.df["Mosaic_path"]= full_paths
        self.logger.info("Adding SUMSS mosaic rms information.")
        self.df["Mosaic_rms"]= rms_values
        
    def _find_nvss_mosaic_info(self, nvss_mosaic_dir):
        if self.survey_name.lower()!="nvss":
            self.logger.error("Catalog is not defined as a NVSS catalog. Cannot add NVSS mosaic directory information.")
            return
        try:
            self.logger.info("Loading NVSS image data.")
            nvss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/nvss_images_info.csv"))
            nvss_centres = SkyCoord(ra=nvss_mosaic_data["center-ra"].values*u.deg, dec=nvss_mosaic_data["center-dec"].values*u.deg)
        except:
            self.logger.error("NVSS mosaic data cannot be found!")
        images = []
        images_full_path = []
        images_rms = []
        for i,row in self.df.iterrows():
            image,rms=self._find_matching_image(row[self.ra_col], row[self.dec_col], nvss_centres, nvss_mosaic_data, rms=True)
            images.append(image)
            images_full_path.append(os.path.join(nvss_mosaic_dir, image))
            images_rms.append(rms)
        self.logger.info("Adding NVSS mosaic information.")
        self.df["Mosaic"] = images
        self.df["Mosaic_path"] = images_full_path
        self.df["Mosaic_rms"] = images_rms
            
    
    def _image_check(self, ra, dec, image, sumss_mosaic_dir="", nvss_mosaic_dir=""):
        if image.startswith("J"):
            full_image = os.path.join(sumss_mosaic_dir, image)
            sumss=True
        elif image.startswith("C"):
            full_image = os.path.join(nvss_mosaic_dir, image)
            sumss=False
        else:
            self.logger.error("Cannot perform image check as image type unrecongised.")
            return True
            
        fits_image=askapimage(full_image)
        fits_image.load_wcs()
        fits_image.load_fits_data()
        naned=self._check_for_nan_image(ra, dec, fits_image.wcs, fits_image.data, sumss=sumss)
        
        #Switch this around to make it more logical - True return = check passed.
        if naned:
            return False
        else:
            return True
   
    def _find_matching_image(self,ra, dec, centres, images_data, rms=False, check=False, attempts=3, sumss_mosaic_dir="", nvss_mosaic_dir=""):
        target=SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        seps = target.separation(centres)
        min_index = np.argmin(seps.deg)
        # print d2d_sumss.deg[min_index]
        # print min_index
        # find the next index if needed
        sorted_seps = sorted(seps.deg)
        backup_indexes = [list(seps.deg).index(i) for i in sorted_seps[1:attempts]]
        image = images_data.iloc[min_index]["image"]
        if check:
            check_result = self._image_check(ra, dec, image, sumss_mosaic_dir=sumss_mosaic_dir, nvss_mosaic_dir=nvss_mosaic_dir)
                # print image
                # print askap_target.to_string('hmsdms')    
            
            if not check_result:
                self.logger.warning("Coordinates out of range for closest image (Coord: {} Image: {}).".format(target.to_string("hmsdms"), image))
                if attempts>1:
                    self.logger.warning("Will try for {} more attempts to find an image.".format(attempts-1))
                    found=False
                    for i in backup_indexes:
                        this_image = images_data.iloc[i]["image"]
                        check_result = self._image_check(ra, dec, this_image, sumss_mosaic_dir=sumss_mosaic_dir, nvss_mosaic_dir=nvss_mosaic_dir)
                        if check_result:
                            self.logger.info("Match found!")
                            image = this_image
                            self.logger.info("New image: {}".format(image))
                            min_index = i
                            found=True
                            break
                        else:
                            continue
                    if not found:
                        self.logger.error("No other suitable image was found.")
                        self.logger.error("No more attempts to be made. Image for {} will be out of range".format(target.to_string('hmsdms')))
                else:
                    self.error("No more attempts to be made. Image for {} will be out of range".format(target.to_string('hmsdms')))
                    
        rms = images_data.iloc[min_index]["rms"]
        # print image
        # print askap_target.to_string('hmsdms')
        if rms:
            return image, rms
        else:
            return image
        
    
    def _check_for_nan_image(self, ra, dec, img_wcs, img_data, num_pixels=10, sumss=False):
        #Works with ASKAP Image for now
        pixels=img_wcs.wcs_world2pix(ra, dec, 1)
        y,x = pixels
        y = int(y)
        x = int(x)
        if y < 0 or x < 0:
            self.logger.error("Pixels out of range - returning nan.")
            return True
        self.logger.debug("RA: {}, Dec:{}".format(ra,dec))
        self.logger.debug("Pixels: {}, {}".format(y,x))
        #Search 10 pixels around, see if in nan is in there
        row_idx = np.array([range(x-num_pixels, x+num_pixels+1)])
        col_idx = np.array([range(y-num_pixels, y+num_pixels+1)])
        try:
            if sumss:
                data_selection=img_data[row_idx[:, None], col_idx]
            else:
                data_selection=img_data[0,0,row_idx[:, None], col_idx]
        except:
            self.logger.error("Measuring failed - returning nan.")
            data_selection=[np.nan]
        self.logger.debug(data_selection)
        return np.isnan(data_selection).any()
    
    def find_matching_catalog_image(self, sumss_mosaic_dir="", nvss_mosaic_dir="", sumss=True, nvss=False, dualmode=False):
        if sumss:
            try:
                self.logger.info("Loading SUMSS image data.")
                sumss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/sumss_images_info.csv"))
                sumss_centres = SkyCoord(ra=sumss_mosaic_data["center-ra"].values*u.deg, dec=sumss_mosaic_data["center-dec"].values*u.deg)
                if not dualmode:
                    centers_to_use = sumss_centres
                    mosaic_data_to_use = sumss_mosaic_data
                    path_dir = sumss_mosaic_dir
            except:
                self.logger.error("SUMSS mosaic data cannot be found!")
                return
    
        if nvss:
                try:
                    self.logger.info("Loading NVSS image data.")
                    nvss_mosaic_data=pd.read_csv(pkg_resources.resource_filename(__name__, "../data/nvss_images_info.csv"))
                    nvss_centres = SkyCoord(ra=nvss_mosaic_data["center-ra"].values*u.deg, dec=nvss_mosaic_data["center-dec"].values*u.deg)
                    if not dualmode:
                        centers_to_use = nvss_centres
                        mosaic_data_to_use = nvss_mosaic_data
                        path_dir = nvss_mosaic_dir
                except:
                    sel.flogger.error("NVSS mosaic data cannot be found!")
                    
        images = []
        images_full_path = []
        images_rms = []
        survey = []
        
        for i,row in self.df.iterrows():
            if dualmode:
                if row[self.dec_col]<=-30:
                    centers_to_use = sumss_centres
                    mosaic_data_to_use = sumss_mosaic_data
                    path_dir = sumss_mosaic_dir
                    survey.append("sumss")
                else:
                    centers_to_use = nvss_centres
                    mosaic_data_to_use = nvss_mosaic_data
                    path_dir = nvss_mosaic_dir
                    survey.append("nvss")
            else:
                if nvss:
                    centers_to_use = nvss_centres
                    mosaic_data_to_use = nvss_mosaic_data
                    path_dir = nvss_mosaic_dir
                    survey.append("nvss")
                elif sumss:
                    centers_to_use = sumss_centres
                    mosaic_data_to_use = sumss_mosaic_data
                    path_dir = sumss_mosaic_dir
                    survey.append("sumss")
                    
                    
            image,rms = self._find_matching_image(row[self.ra_col], row[self.dec_col], centers_to_use, mosaic_data_to_use, rms=True, check=True, sumss_mosaic_dir=sumss_mosaic_dir,
                nvss_mosaic_dir=nvss_mosaic_dir)
            self.logger.debug("{} matched to catalog mosaic {}.".format(row["name"], image))
            images.append(image)
            images_full_path.append(os.path.join(path_dir, image))
            images_rms.append(rms)
           
        self.df["catalog_Mosaic"] = images
        self.df["catalog_Mosaic_path"] = images_full_path
        self.df["catalog_Mosaic_rms"] = images_rms
        self.df["survey_used"] = survey
        

    def _merge_askap_non_convolved_catalogue(self, non_conv_cat, sumss=False, nvss=False):
        """
        This is only to be used with an ASKAP catalogue and is a quick fix to get a match between each convolved ASKAP source, and it's equivalent
        in the non-convolved catalogue.
        Essentially a crossmatch 'lite'.
        """
        idx, d2d, d3d = self._crossmatch_catalog.match_to_catalog_sky(non_conv_cat._crossmatch_catalog)
        matches = non_conv_cat.df.loc[idx].reset_index(drop=True)
        #going to copy over RA, Dec, iflux, iflux_err, a, b and pa
        cols_to_copy = ["name", non_conv_cat.ra_col, non_conv_cat.dec_col, "int_flux", "err_int_flux", "a", "b", "pa"]
        if sumss:
            cols_to_copy+=["askap_scaled_to_sumss"]
        if nvss:
            cols_to_copy+=["askap_scaled_to_nvss"]
        matches = matches.filter(cols_to_copy, axis=1)
        newnames={}
        for i in cols_to_copy:
            newnames[i] = "non_conv_{}".format(i)
        matches=matches.rename(str, columns=newnames)
        # matches["non_conv_name"] = matches['non_conv_name'].astype(basestring)
        matches.index=range(len(matches.index))
        self.df=self.df.merge(matches, left_index=True, right_index=True, how="left")
        self.df["non_conv_d2d"]=d2d.arcsec
        self.logger.info("Non-convolved cross-matches added to ASKAP catalog.")
    


