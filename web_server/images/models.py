# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models

# Create your models here.

class Image(models.Model):
    image_id = models.IntegerField()
    name = models.CharField(max_length=100, unique=False, default='name')
    description = models.CharField(max_length=100, default="ASKAP image")
    ra = models.DecimalField(max_digits=10, decimal_places=7)
    dec = models.DecimalField(max_digits=10, decimal_places=7)
    runtime = models.DateTimeField(auto_now=False, auto_now_add=True)
    url = models.URLField(max_length=300)
    matched_to = models.CharField(max_length=40, default="SUMSS")
    flux_ratio_image_view = models.CharField(max_length=200, unique=False, default="plot")
    position_offset = models.CharField(max_length=200, unique=False, default="plot")
    source_counts = models.CharField(max_length=200, unique=False, default="plot")
    flux_ratios = models.CharField(max_length=200, unique=False, default="plot")
    flux_ratios_from_centre = models.CharField(max_length=200, unique=False, default="plot")
    askap_overlay = models.CharField(max_length=200, unique=False, default="plot")
    sumss_overlay = models.CharField(max_length=200, unique=False, default="plot")
    nvss_overlay = models.CharField(max_length=200, unique=False, default="plot")
    runby = models.CharField(max_length=20, unique=False, default="unknown")
    number_askap_sources=models.IntegerField(default=0)
    number_sumss_sources=models.IntegerField(default=0)
    number_nvss_sources=models.IntegerField(default=0)
    rms=models.DecimalField(max_digits=20, decimal_places=17, default=0.0)
    transients_noaskapmatchtocatalog_total=models.IntegerField(default=0)
    transients_noaskapmatchtocatalog_candidates=models.IntegerField(default=0)
    transients_nocatalogmatchtoaskap_total=models.IntegerField(default=0)
    transients_nocatalogmatchtoaskap_candidates=models.IntegerField(default=0)
    transients_largeratio_total=models.IntegerField(default=0)
    transients_goodmatches_total=models.IntegerField(default=0)
    
    def __str__(self):
        return self.name
    
class Processingsettings(models.Model):
    image_id = models.IntegerField()
    output_dir = models.URLField(max_length=300)
    askap_csv = models.URLField(max_length=300)
    sumss_csv = models.URLField(max_length=300)
    nvss_csv = models.URLField(max_length=300, default="")
    askap_ext_thresh = models.DecimalField(max_digits=4, decimal_places=2, default=0.0)
    sumss_ext_thresh = models.DecimalField(max_digits=4, decimal_places=2, default=0.0)
    nvss_ext_thresh = models.DecimalField(max_digits=4, decimal_places=2, default=0.0)
    max_separation = models.DecimalField(max_digits=4, decimal_places=2)
    aegean_det_sigma = models.DecimalField(max_digits=4, decimal_places=2)
    aegean_fill_sigma = models.DecimalField(max_digits=4, decimal_places=2)

    
class Askapsrc(models.Model):
    askap_source_id = models.IntegerField()
    name = models.CharField(max_length=50, unique=False)
    ra = models.DecimalField(max_digits=10, decimal_places=7)
    dec = models.DecimalField(max_digits=10, decimal_places=7)
    int_flux = models.DecimalField(max_digits=10, decimal_places=7)
    
    def __str__(self):
        return self.name
    
    
class Sumsssrc(models.Model):
    sumss_source_id = models.IntegerField()
    name = models.CharField(max_length=50, unique=False)
    ra = models.DecimalField(max_digits=10, decimal_places=7)
    dec = models.DecimalField(max_digits=10, decimal_places=7)
    int_flux = models.DecimalField(max_digits=10, decimal_places=7)
    
    def __str__(self):
        return self.name
    
    
class Sumssnomatch(models.Model):
    image_id = models.IntegerField()
    match_id = models.IntegerField()
    master_name = models.CharField(max_length=50, unique=False, default="source")
    sumss_name = models.CharField(max_length=50, unique=False, default="sumss source")
    nvss_name = models.CharField(max_length=50, unique=False, default="nvss source")
    ra = models.DecimalField(max_digits=10, decimal_places=7)
    dec = models.DecimalField(max_digits=10, decimal_places=7)
    sumss_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    sumss_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    nvss_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    nvss_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_scale_flux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    sumss_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    nvss_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    cat_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    scaled_cat_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    askap_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    ploturl = models.CharField(max_length=200, unique=False, default="plot")
    pipelinetag = models.CharField(max_length=50, unique=False, default="N/A")
    usertag = models.CharField(max_length=30, unique=False, default="N/A")
    userreason = models.CharField(max_length=30, unique=False, default="N/A")
    checkedby = models.CharField(max_length=20, unique=False, default="N/A")
    survey = models.CharField(max_length=10, unique=False, default="N/A")
    
    def __str__(self):
        return self.sumss_name
        
        
class Largeratio(models.Model):
    image_id = models.IntegerField()
    match_id = models.IntegerField()
    master_name = models.CharField(max_length=50, unique=False, default="source")
    sumss_name = models.CharField(max_length=50, unique=False, default="sumss source")
    nvss_name = models.CharField(max_length=50, unique=False, default="nvss source")
    ra = models.DecimalField(max_digits=10, decimal_places=7)
    dec = models.DecimalField(max_digits=10, decimal_places=7)
    sumss_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    sumss_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    nvss_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    nvss_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_scale_flux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_cat_ratio = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    sumss_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    nvss_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    cat_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    scaled_cat_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    askap_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    ploturl = models.CharField(max_length=200, unique=False, default="plot")
    pipelinetag = models.CharField(max_length=50, unique=False, default="N/A")
    usertag = models.CharField(max_length=30, unique=False, default="N/A")
    userreason = models.CharField(max_length=30, unique=False, default="N/A")
    checkedby = models.CharField(max_length=20, unique=False, default="N/A")
    survey = models.CharField(max_length=10, unique=False, default="N/A")
    
    def __str__(self):
        return self.sumss_name
        

class Askapnotseen(models.Model):
    image_id = models.IntegerField()
    match_id = models.IntegerField()
    askap_name = models.CharField(max_length=50, unique=False, default="askap source")
    ra = models.DecimalField(max_digits=10, decimal_places=7)
    dec = models.DecimalField(max_digits=10, decimal_places=7)
    askap_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_scale_flux_sumss = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_scale_flux_nvss = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    scaled_cat_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    cat_flux = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    # nvss_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    # nvss_flux = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    ploturl = models.CharField(max_length=200, unique=False, default="plot")
    pipelinetag = models.CharField(max_length=50, unique=False, default="N/A")
    usertag = models.CharField(max_length=30, unique=False, default="N/A")
    userreason = models.CharField(max_length=30, unique=False, default="N/A")
    checkedby = models.CharField(max_length=20, unique=False, default="N/A")
    survey = models.CharField(max_length=10, unique=False, default="N/A")
    
    def __str__(self):
        return self.sumss_name
  
class Goodmatch(models.Model):
    image_id = models.IntegerField()
    match_id = models.IntegerField()
    master_name = models.CharField(max_length=50, unique=False, default="source")
    sumss_name = models.CharField(max_length=50, unique=False, default="sumss source")
    nvss_name = models.CharField(max_length=50, unique=False, default="nvss source")
    ra = models.DecimalField(max_digits=10, decimal_places=7)
    dec = models.DecimalField(max_digits=10, decimal_places=7)
    sumss_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    sumss_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    nvss_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    nvss_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    askap_scale_flux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    sumss_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    nvss_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    ploturl = models.CharField(max_length=200, unique=False, default="plot")
    pipelinetag = models.CharField(max_length=50, unique=False, default="N/A")
    usertag = models.CharField(max_length=30, unique=False, default="N/A")
    userreason = models.CharField(max_length=30, unique=False, default="N/A")
    checkedby = models.CharField(max_length=20, unique=False, default="N/A")
    survey = models.CharField(max_length=10, unique=False, default="N/A")
    
    def __str__(self):
        return self.sumss_name  
        
class Query(models.Model):
    transient_type = models.CharField(max_length=50)
    user_tag = models.CharField(max_length=50)
    user = models.CharField(max_length=50)

    # def get_absolute_url(self):
    #     return reverse("queries:detail", kwargs={"id": self.id})
        #return "/queries/%s/" %(self.id)
