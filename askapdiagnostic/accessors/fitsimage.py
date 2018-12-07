#!/usr/bin/env python

from astropy import wcs
from astropy.io import fits
import os
import subprocess
import logging
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import pandas as pd

class askapimage(object):
    """docstring for fitsimage"""
    def __new__(cls,arg,readinfo=False,logger=None):
        if not os.path.isfile(arg):
            raise ValueError("Dataset does not seem to exist, check path!")
            return None
        else:
            return super(askapimage, cls).__new__(cls)
            
    def __init__(self, arg, logger=None, readinfo=False):
        super(askapimage, self).__init__()
        self.logger = logger or logging.getLogger(__name__)
        self.image = arg
        self.imagename=self.image.split("/")[-1]
        if readinfo:
            self.load_all()
                
    def load_wcs(self):
        # Load the FITS hdulist using astropy.io.fits
        with fits.open(self.image) as hdulist:

        # Parse the WCS keywords in the primary HDU
            self.wcs = wcs.WCS(hdulist[0].header, naxis=2)

        # Print out the "name" of the WCS, as defined in the FITS header
        # print(w.wcs.name)

    def load_position_dimensions(self):
        if not self.wcs:
            self.load_wcs()
        if not self.header:
            self.load_header()
        self.size_x=self.header["NAXIS1"]
        self.size_y=self.header["NAXIS2"]
        central_coords=[[self.size_x/2.,self.size_y/2.]]
        center = self.wcs.wcs_pix2world(np.array(central_coords, np.float_), 1)
        self.center=SkyCoord(ra=center[0][0]*u.degree, dec=center[0][1]*u.degree)
        if self.size_x >=self.size_y:
            radius_pixels=[[self.size_x,self.size_y/2.]]
        else:
            radius_pixels=[[self.size_x/2.,self.size_y]]
        outskirt = self.wcs.wcs_pix2world(np.array(radius_pixels, np.float_), 1)
        outskirt_coords=SkyCoord(ra=outskirt[0][0]*u.degree, dec=outskirt[0][1]*u.degree)
        self.radius=self.center.separation(outskirt_coords)
        self.logger.info("Image Center: {}".format(self.center.to_string('hmsdms')))
        self.logger.info("Radius: {} (deg)".format(self.radius.degree))
    
    def load_fits_data(self):
        with fits.open(self.image) as hdul:
            self.data = hdul[0].data
        
    def load_fits_header(self):
        with fits.open(self.image) as hdul:
            self.header = hdul[0].header
            
    def load_all(self):
        with fits.open(self.image) as hdul:
            self.header = hdul[0].header
            self.data = hdul[0].data
            self.wcs = wcs.WCS(self.header, naxis=2)
            self.load_position_dimensions()
            
    def find_sources(self, sf="aegean", options={}):
        if sf=="aegean":
            command="aegean "
            for opt in sorted(options):
                if opt=="table":
                    continue
                command+="--{} {} ".format(opt, options[opt])
            command+="--table {} ".format(self.imagename.replace(".fits", ".csv"))
            command+=self.image
            subprocess.call(command, shell=True)
            
            
    def get_sumss_catalogue(self, boundary_value="nan", write_ann=False):
        cat_ids={"SUMSS":"VIII/81B/sumss212"}
        if not self.wcs:
            self.load_wcs()
        if not self.header:
            self.load_header()
        if not self.center:
            self.load_position_dimensions()
        if np.abs(self.size_x-self.size_y) > max([self.size_x,self.size_y])*0.05:
            self.logger.warning("Non square image detected! Will pad the Vizier search area by 1.2.")
            pad=1.2
        else:
            pad=1.0
        # if not self.data:
        #     self.load_data()
        v = Vizier(columns=["_r", "_RAJ2000","_DEJ2000", "**"])
        v.ROW_LIMIT=-1
        self.logger.info("Querying SUMSS using Vizier...")
        sumss_result=v.query_region(self.center, radius=self.radius*pad, catalog="SUMSS")
        sumss_result=sumss_result[cat_ids["SUMSS"]].to_pandas()
        self.raw_sumss_sources=sumss_result
        self.logger.info("SUMSS sources obtained.")
        self.logger.info("Filtering SUMSS sources to only those within the image area...")
        self.sumss_sources=self._filter_sumss_catalogue(boundary_value=boundary_value)
        if write_ann:
            self._write_ann(self.raw_sumss_sources, self.imagename.replace(".fits", "_raw_sumss_sources.ann" ), "_RAJ2000", "_DEJ2000", "RED")
            self._write_ann(self.sumss_sources, self.imagename.replace(".fits", "_sumss_sources.ann" ), "_RAJ2000", "_DEJ2000", "GREEN")
        return self.sumss_sources
        
        
    def _filter_sumss_catalogue(self, boundary_value="nan"):
        self.logger.info("{} sources before filtering.".format(len(self.raw_sumss_sources.index)))
        #First gather all the coordinates of the SUMSS sources in to an array
        world_coords=[]
        for i, row in self.raw_sumss_sources.iterrows():
            world_coords.append([row["_RAJ2000"],row["_DEJ2000"]])
        
        #Convert these to pixel coordinates on the ASKAP image.
        pixel_coords = self.wcs.wcs_world2pix(np.array(world_coords, np.float_), 1)
        
        #Prepare a record of which ones are equal to the boundary value (default nan)
        nan_mask=[]
        
        #Loop over pixel values checking the value from the image data array.
        for i in pixel_coords.astype(np.int):
            ra=i[0]
            dec=i[1]
            #Get rid of any negative ones, these can wrap around if not handled.
            if (ra < 0) or (dec < 0):
                isnan=False
            elif (ra > self.size_x) or (dec > self.size_y):
                isnan=False
            elif boundary_value=="nan":
                try:
                    if len(self.data.shape)>2:
                        isnan = np.isnan(self.data[0,0,i[1], i[0]])
                    else:
                        isnan = np.isnan(self.data[i[1], i[0]])
                #Some might be outside the image boundary
                except:
                    isnan = False
            else:
                try:
                    if len(self.data.shape)>2:
                        num_nonzero=np.count_nonzero(self.data[0,0,i[1], i[0]])
                    else:
                        num_nonzero=np.count_nonzero(self.data[i[1], i[0]])
                    if num_nonzero==0:
                        isnan=False
                    else:
                        isnan=True
                except:
                    isnan=False
            
            nan_mask.append(isnan)

        #Create new DF
        # mask_sumss=pd.DataFrame().reindex_like(self.raw_sumss_sources)

        # for i, row in self.raw_sumss_sources.iterrows():
            # if nan_mask[i]==False:
                # mask_sumss.iloc[i]=row
                    
        if boundary_value=="nan":
            nan_mask=~np.array(nan_mask)
        mask_sumss=self.raw_sumss_sources[nan_mask]

        mask_sumss=mask_sumss.dropna(how="all")
        
        self.logger.info("{} sources remain after filtering.".format(len(mask_sumss.index)))
        
        return mask_sumss
        
    def _write_ann(self, df,name, ra_col, dec_col, color):
        ra=df[ra_col].values
        dec=df[dec_col].values
        
        catalog = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
        
        with open(name, 'w') as f:
            f.write("COORD W\n")
            f.write("PA STANDARD\n")
            f.write("COLOR {}\n".format(color))
            f.write("FONT hershey14\n")
            for i in catalog:
                f.write("CIRCLE {0} {1} 0.0111111\n".format(i.ra.deg, i.dec.deg))
                
        self.logger.info("Wrote annotation file {}.".format(name))
                
    def write_sumss_sources(self):
        name=self.imagename.replace(".fits", "_sumss_comp.csv")
        self.sumss_sources.to_csv(name, sep=",", index=False)
        self.logger.info("Wrote SUMSS sources to {}.".format(name))
            
                 
            
        
