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
import datetime
import sqlalchemy
import psycopg2
from askapdiagnostic.plotting import plots
from askapdiagnostic.tools import utils
import bdsf

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
        self.original_name = None
        self.non_convolved = None
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
        centre = self.wcs.wcs_pix2world(np.array(central_coords, np.float_), 1)
        self.centre=SkyCoord(ra=centre[0][0]*u.degree, dec=centre[0][1]*u.degree)
        if self.size_x >=self.size_y:
            radius_pixels=[[self.size_x,self.size_y/2.]]
        else:
            radius_pixels=[[self.size_x/2.,self.size_y]]
        outskirt = self.wcs.wcs_pix2world(np.array(radius_pixels, np.float_), 1)
        outskirt_coords=SkyCoord(ra=outskirt[0][0]*u.degree, dec=outskirt[0][1]*u.degree)
        self.radius=self.centre.separation(outskirt_coords)
        self.logger.info("Image Centre: {}".format(self.centre.to_string('hmsdms')))
        self.logger.info("Radius: {} (deg)".format(self.radius.degree))
    
    def load_fits_data(self):
        with fits.open(self.image) as hdul:
            self.data = hdul[0].data
        
    def load_fits_header(self):
        with fits.open(self.image) as hdul:
            self.header = hdul[0].header
        try:
            self.freq = self.header["RESTFRQ"]
            self.logger.info("Frequency = {} MHz".format(self.freq/1.e6))
        except:
            self.logger.warning("Frequency of image couldn't be determined.")
        try:
            self.bmaj = float(self.header["BMAJ"])
            self.bmin = float(self.header["BMIN"])
            self.bpa = float(self.header["BPA"])
            self.logger.info("Beam = {:.2f}\" x {:.2f}\" ({:.2f} deg)".format(self.bmaj*3600., self.bmin*3600., self.bpa))
        except:
            self.logger.warning("Beam information could not be determined.")
         
    def load_all(self):
        self.load_fits_header()
        with fits.open(self.image) as hdul:
            self.data = hdul[0].data
            self.wcs = wcs.WCS(self.header, naxis=2)
            self.load_position_dimensions()
            
    def calculate_sumss_beam(self):
        self.img_sumss_bmaj = 45.*(1./np.sin(np.deg2rad(np.abs(self.centre.dec.degree))))
        self.img_sumss_bmin = 45.
    
    def _get_default_selavy_options(self):
        mapnames=self.imagename.reaplce(".fits", "")
        selavy_options={"Selavy.imagetype":"fits",
        "Selavy.spectralTermsFromTaylor":"true",
        "Selavy.findSpectralTerms":"[false, false]",
        "Selavy.nsubx":1,
        "Selavy.nsuby":1,
        "Selavy.overlapx":0,
        "Selavy.overlapy":0,
        # Detection threshold
        "Selavy.snrCut":5,
        "Selavy.flagGrowth":"true",
        "Selavy.growthThreshold":3,
        #
        "Selavy.VariableThreshold":"true",
        "Selavy.VariableThreshold.boxSize":50,
        "Selavy.VariableThreshold.ThresholdImageName":"detThresh.{}".format(mapnames),
        "Selavy.VariableThreshold.NoiseImageName":"noiseMap.{}".format(mapnames),
        "Selavy.VariableThreshold.AverageImageName":"meanMap.{}".format(mapnames),
        "Selavy.VariableThreshold.SNRimageName":"snrMap.{0}".format(mapnames),
        "Selavy.VariableThreshold.imagetype":"fits",
        # Selavy.Weights.weightsImage                     = field_RACS_test2_1.05_2030-68A.altpb.weight.fits
        # Selavy.Weights.weightsCutoff                    = 0.04
        #
        "Selavy.Fitter.doFit":"true",
        "Selavy.Fitter.fitTypes":"[full]",
        "Selavy.Fitter.numGaussFromGuess":"true",
        "Selavy.Fitter.maxReducedChisq":10.0,
        "Selavy.Fitter.imagetype":"fits",
        #
        "Selavy.threshSpatial":5,
        "Selavy.flagAdjacent":"true",
        #
        "Selavy.minPix":3,
        "Selavy.minVoxels":3,
        "Selavy.minChannels":1,
        "Selavy.sortingParam":"-pflux",
        #
        # Not performing RM Synthesis for this case
        "Selavy.RMSynthesis":"false"}
        return selavy_options
        
    def _write_selavy_parset(self, parset_out, settings_to_write):
        with open(parset_out, 'w') as f:
            f.write("# selavy parset for {}\n\n".format(self.imagename))
            for i in sorted(settings_to_write):
                f.write("{0:<40} = {1}\n".format(i, settings_to_write[i]))
    
    def find_sources(self, sf="aegean", options={}):
        outfile=self.imagename.replace(".fits", "_comp.csv")
        if sf=="aegean":
            command="aegean "
            for opt in sorted(options):
                if opt=="table":
                    continue
                command+="--{} {} ".format(opt, options[opt])
            command+="--table {} ".format(self.imagename.replace(".fits", ".csv"))
            command+=self.image+" "
            command+="> "+self.imagename.replace(".fits", "_aegean.log")
            self.logger.debug("Aegean command: {}".format(command))
            subprocess.call(command, shell=True)

        elif sf=="pybdsf":
            self.logger.info("Running source finding with PyBDSF...")
            kwargs = options
            # procimg_kwargs.update(options)
            pybdsf_img = bdsf.process_image(self.image, **kwargs)
            pybdsf_img.write_catalog(outfile=outfile+".original", catalog_type="gaul", format="csv")
            utils.pybdsf2aegean(outfile+".original", outfile)
            self.logger.info("PyBDSF output converted to aegean format and placed in {}.".format(outfile))
            

        elif sf=="selavy":
            self.logger.info("Running source finding with Selavy...")
            selavy_parset_name = self.imagename.replace(".fits", ".selavy.in")
            #Add the image to the parset
            selavy_outputname=self.imagename.replace(".fits", ".txt")
            options["Selavy.image"]=self.image
            options["Selavy.resultsFile"]=selavy_outputname
            self._write_selavy_parset(self, selavy_parset_name, options)
            self.logger.info("Written Selavy parset file {}.".format(selavy_parset_name))
            if "Selavy.nsubx" in options:
                nsubx=int(options["Selavy.nsubx"])
            if "Selavy.nsuby" in options:
                nsuby=int(options["Selavy.nsuby"])
                
            num_procs = nsubx * nsuby +1
            
            command = "mpirun -n {} selavy -c {} > ".format(num_procs, selavy_parset_name)+self.imagename.replace(".fits", "_selavy.log")
            subprocess.call(command, shell=True)
            
            utils.selavy2aegean(selavy_outputname.replace(".txt", ".components.txt"), outfile)
            self.logger.info("Selavy output converted to aegean format and placed in {}.".format(outfile))
            
            
            
    def get_sumss_catalogue(self, boundary_value="nan"):
        cat_ids={"SUMSS":"VIII/81B/sumss212"}
        if not self.wcs:
            self.load_wcs()
        if not self.header:
            self.load_header()
        if not self.centre:
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
        sumss_result=v.query_region(self.centre, radius=self.radius*pad, catalog="SUMSS")
        sumss_result=sumss_result[cat_ids["SUMSS"]].to_pandas()
        self.raw_sumss_sources=sumss_result
        self.logger.info("SUMSS sources obtained.")
        self.logger.info("Filtering SUMSS sources to only those within the image area...")
        self.sumss_sources=self._filter_sumss_catalogue(self.raw_sumss_sources, boundary_value=boundary_value)
        return self.sumss_sources
        
    def get_nvss_catalogue(self, boundary_value="nan"):
        cat_ids={"NVSS":"VIII/65/nvss"}
        if not self.wcs:
            self.load_wcs()
        if not self.header:
            self.load_header()
        if not self.centre:
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
        self.logger.info("Querying NVSS using Vizier...")
        nvss_result=v.query_region(self.centre, radius=self.radius*pad, catalog="NVSS")
        nvss_result=nvss_result[cat_ids["NVSS"]].to_pandas()
        self.raw_nvss_sources=nvss_result
        self.logger.info("NVSS sources obtained.")
        self.logger.info("Filtering NVSS sources to only those within the image area...")
        self.nvss_sources=self._filter_sumss_catalogue(self.raw_nvss_sources, boundary_value=boundary_value)
        return self.sumss_sources
        
        
    def _filter_sumss_catalogue(self, tofilter, boundary_value="nan"):
        self.logger.info("{} sources before filtering.".format(len(tofilter.index)))
        #First gather all the coordinates of the SUMSS sources in to an array
        world_coords=[]
        for i, row in tofilter.iterrows():
            world_coords.append([row["_RAJ2000"],row["_DEJ2000"]])
        
        #Convert these to pixel coordinates on the ASKAP image.
        pixel_coords = self.wcs.wcs_world2pix(np.array(world_coords, np.float_), 1)
        
        #Prepare a record of which ones are equal to the boundary value (default nan)
        nan_mask=[]
        self.logger.debug("Size X: {}, Size Y: {}".format(self.size_x, self.size_y))
        #Loop over pixel values checking the value from the image data array.
        for i in pixel_coords.astype(np.int):
            self.logger.debug(i)
            ra=i[0]
            dec=i[1]
            #Get rid of any negative ones, these can wrap around if not handled.
            if (ra < 0) or (dec < 0):
                isnan=True
            elif (ra > self.size_x) or (dec > self.size_y):
                isnan=True
            elif boundary_value=="nan":
                try:
                    if len(self.data.shape)>2:
                        isnan = np.isnan(self.data[0,0,i[1], i[0]])
                    else:
                        isnan = np.isnan(self.data[i[1], i[0]])
                #Some might be outside the image boundary
                except:
                    isnan = True
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
            self.logger.debug("is nan: {}".format(isnan))
            nan_mask.append(isnan)

        #Create new DF
        # mask_sumss=pd.DataFrame().reindex_like(self.raw_sumss_sources)

        # for i, row in self.raw_sumss_sources.iterrows():
            # if nan_mask[i]==False:
                # mask_sumss.iloc[i]=row
                    
        if boundary_value=="nan":
            nan_mask=~np.array(nan_mask)
        mask_filter=tofilter[nan_mask]

        mask_catalog=mask_filter.dropna(how="all").reset_index(drop=True)
        
        self.logger.info("{} sources remain after filtering.".format(len(mask_catalog.index)))
        
        return mask_catalog
        
                
    def write_sumss_sources(self):
        name=self.imagename.replace(".fits", "_sumss_comp.csv")
        self.sumss_sources.to_csv(name, sep=",", index=False)
        self.logger.info("Wrote SUMSS sources to {}.".format(name))
        
    def write_nvss_sources(self):
        name=self.imagename.replace(".fits", "_nvss_comp.csv")
        self.nvss_sources.to_csv(name, sep=",", index=False)
        self.logger.info("Wrote NVSS sources to {}.".format(name))
        
    def inject_db(self, datestamp=datetime.datetime.utcnow(), user="unknown", description="", db_engine="postgresql", 
            db_username="postgres", db_host="localhost", db_port="5432", db_database="postgres"):
        image_id=1
        engine = sqlalchemy.create_engine('{}://{}@{}:{}/{}'.format(db_engine, db_username, db_host, db_port, db_database))
        result = engine.execute("SELECT id FROM images_image")
        try:
            image_id+=int(result.fetchall()[-1][0])
        except:
            image_id=1
        #Better to control the order
        plots_columns=["flux_ratio_image_view", "position_offset", "source_counts", "flux_ratios", "flux_ratios_from_centre", "askap_overlay", "sumss_overlay"]
        plots_values=["media/{}/{}".format(image_id, self.plots[i]) for i in plots_columns]
        self.logger.info("Image run assigned id {}".format(image_id))
        # conn = psycopg2.connect("host=localhost dbname=RACS user=aste7152 port=5434")
        if self.original_name is not None:
            thisimagename=self.original_name
        else:
            thisimagename=self.imagename
        tempdf=pd.DataFrame([[image_id, thisimagename, description, self.centre.ra.degree, 
            self.centre.dec.degree, datestamp, self.image]+plots_values+[user,self.total_askap_sources, self.total_sumss_sources, self.rms]], columns=["image_id", "name", "description", 
                "ra", "dec", "runtime", "url"]+plots_columns+["runby", "number_askap_sources", "number_sumss_sources", "rms"])
        tempdf.to_sql("images_image", engine, if_exists="append", index=False)
        return image_id
        
    def inject_processing_db(self, image_id, output, askap_cat_file, sumss_source_cat, askap_ext_thresh, sumss_ext_thresh, max_separation, aegean_sigmas, db_engine="postgresql", 
            db_username="postgres", db_host="localhost", db_port="5432", db_database="postgres"):
        engine = sqlalchemy.create_engine('{}://{}@{}:{}/{}'.format(db_engine, db_username, db_host, db_port, db_database))
        settings_columns=["image_id", "output_dir", "askap_csv", "sumss_csv", "askap_ext_thresh", "sumss_ext_thresh", "max_separation", "aegean_det_sigma", "aegean_fill_sigma"]
        settings_data=[image_id, output, askap_cat_file, sumss_source_cat, askap_ext_thresh, sumss_ext_thresh, max_separation]+aegean_sigmas
        tempdf=pd.DataFrame([settings_data], columns=settings_columns)
        tempdf.to_sql("images_processingsettings", engine, if_exists="append", index=False)
        
    def create_overlay_plot(self, overlay_cat, overlay_cat_label="sources", overlay_cat_2=None, overlay_cat_label_2=None, sumss=False):
        thename=plots.image_sources_overlay(self.image, self.imagename, overlay_cat, overlay_cat_label=overlay_cat_label, overlay_cat_2=overlay_cat_2, overlay_cat_label_2=overlay_cat_label_2, sumss=sumss)
        return thename
        
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
            med = np.median(arr, axis=axis)
            std = np.std(arr, axis=axis)
            if np.any(ncount-ncount_new > xtol*ncount):
                self.logger.warning("xtol reached {}; breaking at iteration {}".format(1-1.*ncount_new/ncount, niter+1))
                break
            if np.any(ncount_old-ncount_new < ftol*ncount_old):
                self.logger.warning("ftol reached {}; breaking at iteration {}".format(1-1.*ncount_new/ncount_old, niter+1) )
                break
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

    def get_rms_clipping(self, max_iter=10, sigma=4):
        fln=fits.open(self.image)
        rawdata=np.nan_to_num(fln[0].data)
        angle=fln[0].header['crval1']
        bscale=fln[0].header['bscale']
        rawdata=rawdata.squeeze()
        rawdata=rawdata*bscale
        while len(rawdata) < 20:
            rawdata = rawdata[0]
        X,Y = np.shape(rawdata)
        # rawdata = rawdata[Y/6:5*Y/6,X/6:5*X/6]
        rawdata = rawdata[2*Y/6:4*Y/6,2*X/6:4*X/6]
        orig_raw = rawdata
        med, std, mask = self._Median_clip(rawdata, full_output=True, ftol=0.0, max_iter=max_iter, sigma=sigma)
        rawdata[mask==False] = med
        self.logger.info("{0} estimate rms: {1} Jy".format(self.imagename,std))
        fln.close()
        self.rms = std
        
    def weight_crop(self, weight_image, weight_value):
        with fits.open(weight_image) as weight:
            weight_max=np.max(weight[0].data)
       
            mask=np.where(weight[0].data > weight_max * weight_value, False, True)

        with fits.open(self.image) as fitsimage:
            fitsdata=fitsimage[0].data
        
            fitsdata[mask]=np.nan
            fitsimage[0].data=fitsdata
        
            weight_cropped_output=self.imagename.replace(".fits", ".weightcrop_{}.fits".format(weight_value))
        
            fitsimage.writeto(weight_cropped_output, overwrite=True)
        
        self.weight_cropped_image = weight_cropped_output
        
        return os.path.abspath(weight_cropped_output)
        
    def convolve2sumss(self):
        
        convolve_outname = self.imagename.replace(".fits", ".convolve2sumss_casa.fits")
        script="""#!/usr/bin/env python

import sys
import numpy as np
import glob

def calc_sumss_beam(dec):
    sumss_beam_maj=45.*(1./np.sin(np.deg2rad(np.abs(dec))))
    return sumss_beam_maj

def convolve_img(fitsfile):
    dec = {0}
    sumss_beam=calc_sumss_beam(float(dec))
    outname=fitsfile.replace(".fits", ".convolve2sumss_casa.image")
    imsmooth(imagename=fitsfile, major="{{}}arcsec".format(sumss_beam), minor="45arcsec", pa="0deg", outfile=outname, targetres=True)
    exportfits(imagename=outname, fitsimage=outname.replace(".image", ".fits"))
    

convolve_img("{1}")""".format(self.centre.dec.deg, self.image)

        with open("convolve2casa.py", "w") as f:
            f.write(script)
            
        cmd = "casa --nologger --nogui -c convolve2casa.py"
        try:
            self.logger.info("Running CASA to convolve the image")
            p = subprocess.call(cmd, shell=True)
            success=True
        except subprocess.CalledProcessError as e:
            self.logger.error(e)
            self.logger.error("There was a problem with {}!".format(thisloggermsg))
            success=False
        except KeyboardInterrupt:
            self.logger.error("Keyboard Interrupt caught, will terminate")
            p.terminate()
            success=False
        
        if success:
            with fits.open(convolve_outname, mode="update") as convolved_img:
                convolved_img[0].header["RESTFRQ"]=self.freq
                
            return convolve_outname
            
        else:
            return "Failed"
        
        
        
        
            
                 
            
        
