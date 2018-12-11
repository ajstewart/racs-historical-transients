#!/usr/bin/env python

import logging


class crossmatch(object):
    """docstring for crossmatch"""
    def __init__(self, base_catalog, comp_catalog, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        super(crossmatch, self).__init__()
        self.base_catalog=base_catalog
        self.comp_catalog=comp_catalog
        self.performed=False
        
        
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
            newcolnames[c]="sumss_{}".format(c)
        self.crossmatch_df=self.crossmatch_df.rename(str, columns=newcolnames)
        self.crossmatch_df=self.crossmatch_df.join(self.matches)
        self.crossmatch_df["d2d"]=self.d2d.arcsec
        self.performed=True
        self.logger.info("Crossmatch complete.")
        
    def write_crossmatch(self, outname):
        if self.performed:
            self.crossmatch_df.to_csv(outname, index=False)
            self.logger.info("Written crossmatch dataframe to {}.".format(outname))
        else:
            self.logger.error("Need to perform a cross match first!")
        
        
        
        
    
        
        
