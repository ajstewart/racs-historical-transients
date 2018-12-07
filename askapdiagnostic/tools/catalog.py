#!/usr/bin/env python

import logging
import astropy.units as u
from astropy.coordinates import SkyCoord

class catalog(object):
    """docstring for survey"""
    def __init__(self, df, ra_col="_RAJ2000", dec_col="_DEJ2000", logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self.df = df
        self.ra_col=ra_col
        self.dec_col=dec_col
        self._gen_catalog_for_crossmatch()
        super(catalog, self).__init__()

        
    def _gen_catalog_for_crossmatch(self):
        cat_ra=self.df[self.ra_col].values
        cat_dec=self.df[self.dec_col].values
        
        self.crossmatch_catalog = SkyCoord(ra=cat_ra*u.degree, dec=cat_dec*u.degree)
        
        
    def get_col(self, columname):
        return self.df[columname].values
        
        
    



        
        
    
        
        
