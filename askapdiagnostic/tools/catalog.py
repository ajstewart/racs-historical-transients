#!/usr/bin/env python

import logging
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np

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
        
    def add_telescope_beam_columns(self, mode, manual=False, manual_bmaj=45.0, manual_bmin=45.0):
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
        
    def add_sumss_sn(self, flux_col="int_flux", dec_col="dec", flux_scaling=1.):
        self.df["sumss_snr"]=(flux_scaling*self.df[flux_col])/self.df[dec_col].apply(self._sumss_rms)
        
    def add_nvss_sn(self, flux_col="S1.4", dec_col="_DEJ2000", flux_scaling=1.):
        self.df["nvss_snr"]=(flux_scaling*self.df[flux_col])/0.5
        
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
                label=to_catalog
        else:
            this_freq=frequency
            label=label
        
            if this_freq==-1:
                self.logger.error("No catalog or frequency defined for flux scaling.")
                return
                
        
        self.df["{}_scaled_to_{}".format(self.survey_name, label)]=self._scale_flux(this_freq, self.frequency, flux_col, si=si)
        self.logger.info("Fluxes scaled to {} catalogue frequency ({} MHz)".format(label, this_freq/1.e6))
    


