#!/usr/bin/env python

import logging
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np

class Catalog(object):
    """docstring for survey"""
    def __init__(self, df, surveyref, ra_col="_RAJ2000", dec_col="_DEJ2000", add_name_col=False, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self.df = df
        self.surveyref=surveyref
        self.ra_col=ra_col
        self.dec_col=dec_col
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
        
    def _sumss_rms(self,dec):
        if dec>-50:
            return 0.002
        else:
            return 0.0012
    
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
        
    def remove_extended(self, threshold=1.2, ellipse_a="MajAxis", ellipse_b="MinAxis", beam_a=45., beam_b=45., ellipse_unit="arcsec", sumss_psf=False):
        if sumss_psf:
            sumss_bmaj = 45.*(1./np.sin(np.deg2rad(np.abs(self.df[self.dec_col]))))
            sumss_bmin = 45.
            self.sumss_no_ext_cat = self.df[(self.df[ellipse_a] <= threshold*sumss_bmaj) & 
                (self.df[ellipse_b] <= threshold*sumss_bmin)].reset_index(drop=True)
            self.sumss_ext_cat = self.df[(self.df[ellipse_a] > threshold*sumss_bmaj) | 
                (self.df[ellipse_b] > threshold*sumss_bmin)].reset_index(drop=True)
            return self.sumss_no_ext_cat
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
            name=self.surveyref+".ann"
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
        
    def _add_askap_sn(self):
        try:
            self.df["snr"]=self.df["int_flux"]/self.df["local_rms"]
        except:
            self.logger.error("Adding ASKAP SNR not supported for this dataframe.")
        
        
    



        
        
    
        
        
