#!/usr/bin/env python

import logging



class crossmatch(object):
    """docstring for crossmatch"""
    def __init__(self, base_catalog, comp_catalog, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        super(crossmatch, self).__init__()
        self.base_catalog=base_catalog
        self.comp_catalog=comp_catalog
        
        
    def perform_crossmatch(self):
        self.idx, self.d2d, self.d3d = self.base_catalog.match_to_catalog_sky(self.comp_catalog)
        self.matches = self.comp_catalog[self.idx]
        self.logger.info("Crossmatch complete.")
        
        
        
    
        
        
