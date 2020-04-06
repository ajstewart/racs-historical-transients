# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models
from django.contrib.auth.models import User
from django.db.models.signals import post_save
from django.dispatch import receiver

# Create your models here.

class Image(models.Model):
    # image_id = models.IntegerField()
    unique_tag = models.CharField(max_length=100, unique=False, default='name')
    name = models.CharField(max_length=100, unique=False, default='name')
    description = models.CharField(max_length=100, default="ASKAP image")
    ra = models.DecimalField(max_digits=10, decimal_places=7)
    dec = models.DecimalField(max_digits=10, decimal_places=7)
    runtime = models.DateTimeField(auto_now=False, auto_now_add=True)
    url = models.URLField(max_length=300)
    url_2 = models.URLField(max_length=300, default="N/A")
    matched_to = models.CharField(max_length=40, default="SUMSS")
    runby = models.CharField(max_length=20, unique=False, default="unknown")
    number_askap_sources=models.IntegerField(default=0)
    number_sumss_sources=models.IntegerField(default=0)
    number_nvss_sources=models.IntegerField(default=0)
    rms=models.DecimalField(max_digits=20, decimal_places=17, default=0.0)
    transients_master_total=models.IntegerField(default=0)
    claimed_by=models.CharField(max_length=20, unique=False, default="Unclaimed")
    number_candidates_checked = models.IntegerField(default=0)

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
    flux_ratio_image_view = models.CharField(max_length=200, unique=False, default="plot")
    position_offset = models.CharField(max_length=200, unique=False, default="plot")
    source_counts = models.CharField(max_length=200, unique=False, default="plot")
    flux_ratios = models.CharField(max_length=200, unique=False, default="plot")
    flux_ratios_from_centre = models.CharField(max_length=200, unique=False, default="plot")
    askap_overlay = models.CharField(max_length=200, unique=False, default="plot")
    sumss_overlay = models.CharField(max_length=200, unique=False, default="plot")
    nvss_overlay = models.CharField(max_length=200, unique=False, default="plot")


class Crossmatches(models.Model):
    image_id = models.IntegerField()
    # match_id = models.IntegerField()
    master_name = models.CharField(max_length=50, unique=False, default="askap source")
    askap_name = models.CharField(max_length=50, unique=False, default="catalog source")
    catalog_name = models.CharField(max_length=50, unique=False, default="catalog source")
    ra = models.DecimalField(max_digits=10, decimal_places=7)
    dec = models.DecimalField(max_digits=10, decimal_places=7)
    gal_l = models.DecimalField(max_digits=10, decimal_places=7)
    gal_b = models.DecimalField(max_digits=10, decimal_places=7)
    catalog_iflux = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    catalog_iflux_e = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    catalog_rms = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    # catalog_local_rms = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    catalog_mosaic = models.CharField(max_length=50, unique=False, default="catalog mosaic")
    askap_image = models.CharField(max_length=250, unique=False, default="askap_image")
    askap_image_2 = models.CharField(max_length=250, unique=False, default="askap_image_2")
    askap_iflux = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    askap_iflux_e = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    askap_scale_flux = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    askap_scale_flux_e = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    askap_non_conv_flux = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    askap_non_conv_flux_e = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    askap_non_conv_scaled_flux = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    askap_non_conv_scaled_flux_e = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    askap_non_conv_d2d = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    askap_rms = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    askap_rms_2 = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    d2d_askap_centre = models.DecimalField("Distance to ASKAP Image Centre", max_digits=20, decimal_places=12, default=0)
    ratio = models.DecimalField("Flux Ratio", max_digits=20, decimal_places=12, default=0)
    ratio_e = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    ratio_catalog_flux = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    ratio_catalog_flux_err = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    ratio_askap_flux = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    ratio_askap_flux_err = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    ploturl = models.CharField(max_length=200, unique=False, default="plot")
    pipelinetag = models.CharField("Pipeline Tag", max_length=50, unique=False, default="N/A")
    usertag = models.CharField("User Tag", max_length=30, unique=False, default="N/A")
    userreason = models.CharField("User Reason", max_length=30, unique=False, default="N/A")
    checkedby = models.CharField("Checked By", max_length=20, unique=False, default="N/A")
    survey = models.CharField(max_length=10, unique=False, default="N/A")
    nearest_sources = models.CharField(max_length=250, unique=False, default="")
    transient_type = models.CharField(max_length=50, unique=False, default="")
    aegean_rms_used = models.CharField(max_length=6, unique=False, default="False")
    inflated_convolved_flux = models.CharField("Flux Convolved Error", max_length=6, unique=False, default="False")
    vs = models.DecimalField("Vs", max_digits=20, decimal_places=12, default=0)
    m = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    m_abs = models.DecimalField(max_digits=20, decimal_places=12, default=0)
    using_preconv = models.CharField(max_length=6, unique=False, default="False")
    d2d_nn_askap_cat = models.DecimalField("Distance to nearest ASKAP neighbour", max_digits=20, decimal_places=12, default=0)
    flagged = models.BooleanField(default=False)

    def __str__(self):
        return self.catalog_name

class Query(models.Model):
    transient_type = models.CharField(max_length=50)
    user_tag = models.CharField(max_length=50)
    user = models.CharField(max_length=50)

class Profile(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE)
    last_query = models.TextField(max_length=200, blank=True)

@receiver(post_save, sender=User)
def create_user_profile(sender, instance, created, **kwargs):
    if created:
        Profile.objects.create(user=instance)

@receiver(post_save, sender=User)
def save_user_profile(sender, instance, **kwargs):
    instance.profile.save()
