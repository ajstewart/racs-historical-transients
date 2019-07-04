#!/usr/bin/env python

import django_tables2 as tables
from .models import Image, Sumssnomatch, Largeratio, Goodmatch, Askapnotseen, Transients
from django_tables2.utils import A
from templatetags import units

class RAColumn(tables.Column):
    def render(self, value):
        return units.deg_to_hms(value)
    def value(self, value):
        return value

class DecColumn(tables.Column):
    def render(self, value):
        return units.deg_to_dms(value)
    def value(self, value):
        return value
        
class RMSColumn(tables.Column):
    def render(self, value):
        return "{:.2f}".format(units.jy_to_mjy(value))
        
class FloatColumn(tables.Column):
    def render(self, value):
        return "{:.2f}".format(value)
        
class CapitalColumn(tables.Column):
    def render(self, value):
        return value.upper()
    

class ImageTable(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    name = tables.LinkColumn('image_detail', args=[A('pk')], orderable=True,)
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    # rms = RMSColumn(verbose_name="RMS (mJy)")
    matched_to = tables.Column(verbose_name= 'Matched To')
    number_askap_sources = tables.Column(verbose_name= '# ASKAP Sources')
    number_sumss_sources = tables.Column(verbose_name= '# SUMSS Sources')
    number_nvss_sources = tables.Column(verbose_name= '# NVSS Sources')
    transients_noaskapmatchtocatalog_total = tables.Column(verbose_name= '# No ASKAP match to Catalog')
    transients_nocatalogmatchtoaskap_total = tables.Column(verbose_name= '# No Catalog match to ASKAP')
    transients_largeratio_total = tables.Column(verbose_name= '# Large Ratio')
    transients_goodmatches_total = tables.Column(verbose_name= '# Good Matches')
    runby = tables.Column(verbose_name= 'Run By')
    class Meta:
        model = Image
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id", "name", "description", "ra", "dec", "matched_to", "number_askap_sources", "number_sumss_sources", "number_nvss_sources", "transients_noaskapmatchtocatalog_total",
            "transients_nocatalogmatchtoaskap_total", "transients_largeratio_total", "transients_goodmatches_total", "runtime", "runby" )
        attrs = {"th":{"bgcolor":"#EBEDEF"},}
        
        
class SumssNoMatchListTable(tables.Table):
    export_formats = ['csv',]
    id = tables.Column(verbose_name= 'ID')
    image_id = tables.LinkColumn('image_detail', args=[A('image_id'),], orderable=True, verbose_name= 'Img. ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "noaskapmatchtocatalog", A('id')], orderable=True, verbose_name= 'Source Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    # ra_decimal = table.Column(verbose_name= 'RA dec', accessor=A('ra'))
    # askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    catalog_iflux = FloatColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    cat_snr = tables.Column(verbose_name= 'Catalog SNR')
    scaled_cat_snr = tables.Column(verbose_name= 'Scaled ASKAP SNR')
    survey = CapitalColumn(verbose_name= 'Ref Survey')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    class Meta:
        model = Sumssnomatch
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id", "image_id", "master_name", "ra", "dec", "catalog_iflux", "cat_snr", "scaled_cat_snr", "survey", "pipelinetag", "usertag", "userreason", "checkedby")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}
    
class LargeRatioListTable(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    image_id = tables.LinkColumn('image_detail', args=[A('image_id'),], orderable=True, verbose_name= 'Img. ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "largeratio", A('id')], orderable=True, verbose_name= 'Source Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    catalog_iflux = FloatColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    cat_snr = tables.Column(verbose_name= 'Catalog SNR')
    askap_cat_ratio = FloatColumn(verbose_name= 'ASKAP / Catalog')
    survey = CapitalColumn(verbose_name= 'Ref Survey')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    class Meta:
        model = Largeratio
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","image_id", "master_name", "ra", "dec", "askap_iflux", "catalog_iflux", "cat_snr", "askap_cat_ratio", "survey", "pipelinetag", "usertag", "userreason", "checkedby")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}   
        
class GoodMatchListTable(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    image_id = tables.LinkColumn('image_detail', args=[A('image_id'),], orderable=True, verbose_name= 'Img. ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "goodmatch", A('id')], orderable=True, verbose_name= 'Source Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    survey = CapitalColumn(verbose_name= 'Survey')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    catalog_iflux = FloatColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    sumss_snr = tables.Column(verbose_name= 'SUMSS SNR')
    nvss_snr = tables.Column(verbose_name= 'NVSS SNR')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    class Meta:
        model = Goodmatch
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id", "image_id", "master_name", "ra", "dec", "askap_iflux", "catalog_iflux", "survey", "sumss_snr", "nvss_snr", "pipelinetag", "usertag", "userreason", "checkedby")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}   
        
class AskapNotSeenListTable(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    image_id = tables.LinkColumn('image_detail', args=[A('image_id'),], orderable=True, verbose_name= 'Img. ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "nocatalogmatchtoaskap", A('id')], orderable=True, verbose_name= 'ASKAP Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    askap_snr = tables.Column(verbose_name= 'Local ASKAP SNR')
    # catalog_iflux = tables.Column(verbose_name= 'Cat. Int. Flux (mJy)')
    scaled_askap_snr = tables.Column(verbose_name= 'Scaled Catalog SNR')
    survey = CapitalColumn(verbose_name= 'Ref Survey')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    class Meta:
        model = Askapnotseen
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","image_id","master_name", "ra", "dec", "askap_iflux", "askap_snr", "scaled_askap_snr", "survey", "pipelinetag", "usertag", "userreason", "checkedby")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}   
    
    
    
class CrossmatchDetailFluxTable(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    master_name = tables.Column(verbose_name = 'Name')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    askap_scale_flux = RMSColumn(verbose_name= 'Scaled ASKAP Int. Flux (mJy)')
    # measured_askap_local_rms = RMSColumn(verbose_name= 'Local RMS (mJy)')
    # measured_askap_local_rms_2 = RMSColumn(verbose_name= 'Non-convolved local RMS (mJy)')
    catalog_iflux = RMSColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    # catalog_scale_flux = FloatColumn(verbose_name= 'Scaled Cat. Int. Flux (mJy)')
    ratio = FloatColumn(verbose_name= 'Int. Flux Ratio')
    askap_non_conv_flux = RMSColumn(verbose_name= 'Non-convolved Int. Flux (mJy)')
    d2d_askap_centre = FloatColumn(verbose_name = "Distance from ASKAP Centre (deg)")
    # askap_non_conv_scaled_flux = RMSColumn(verbose_name= 'Scaled Non-convolved Int. Flux (mJy)')
    # askap_non_conv_d2d = FloatColumn(verbose_name= 'Distance to ASKAP Convolved Source (arcsec)')
    survey = CapitalColumn(verbose_name= 'Survey Used')

    class Meta:
        model = Sumssnomatch
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","master_name", "askap_iflux", "askap_scale_flux", "catalog_iflux", 
        "ratio", "askap_non_conv_flux", "d2d_askap_centre", "survey")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}   
        
        
class NearestSourceDetailFluxTable(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "transients", A('id')], orderable=True, verbose_name= 'Name')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    askap_scale_flux = RMSColumn(verbose_name= 'Scaled ASKAP Int. Flux (mJy)')
    catalog_iflux = RMSColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    ratio = FloatColumn(verbose_name= 'ASKAP / Cat Int. Flux Ratio')
    survey = CapitalColumn(verbose_name= 'Survey Used')

    class Meta:
        model = Sumssnomatch
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","master_name", "askap_iflux", "askap_scale_flux", "catalog_iflux", "ratio", "survey")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}  
        
        
class TransientTable(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    image_id = tables.LinkColumn('image_detail', args=[A('image_id'),], orderable=True, verbose_name= 'Img. ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "transients", A('id')], orderable=True, verbose_name= 'Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    ratio = FloatColumn(verbose_name='Freq. Scaled Ratio')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    # askap_snr = tables.Column(verbose_name= 'Local ASKAP SNR')
    catalog_iflux = RMSColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    # scaled_askap_snr = tables.Column(verbose_name= 'Scaled Catalog SNR')
    survey = CapitalColumn(verbose_name= 'Ref Survey')
    transient_type = tables.Column(verbose_name= 'Type')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    class Meta:
        model = Transients
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","image_id","master_name", "ra", "dec", "askap_iflux","catalog_iflux", "ratio", "survey", "transient_type", "pipelinetag", "usertag", "userreason", "checkedby")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}  
  
    
    # image_id = models.IntegerField()
    # match_id = models.IntegerField()
    # sumss_name = models.CharField(max_length=50, unique=False, default="sumss source")
    # ra = models.DecimalField(max_digits=10, decimal_places=7)
    # dec = models.DecimalField(max_digits=10, decimal_places=7)
    # sumss_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    # sumss_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    # askap_iflux = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    # askap_iflux_e = models.DecimalField(max_digits=20, decimal_places=3, default=0)
    # sumss_snr = models.DecimalField(max_digits=20, decimal_places=2, default=0)
    # ploturl = models.CharField(max_length=200, unique=False, default="plot")
    # pipelinetag = models.CharField(max_length=30, unique=False, default="N/A")
    # usertag = models.CharField(max_length=30, unique=False, default="N/A")
    # checkedby = models.CharField(max_length=20, unique=False, default="N/A")


