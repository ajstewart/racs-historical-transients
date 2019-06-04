#!/usr/bin/env python

import django_tables2 as tables
from .models import Image, Sumssnomatch, Largeratio, Goodmatch, Askapnotseen
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
        return "{:.3f}".format(units.jy_to_mjy(value))
    

class ImageTable(tables.Table):
    image_id = tables.Column(verbose_name= 'ID')
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
        fields = ("image_id", "name", "description", "ra", "dec", "matched_to", "number_askap_sources", "number_sumss_sources", "number_nvss_sources", "transients_noaskapmatchtocatalog_total",
            "transients_nocatalogmatchtoaskap_total", "transients_largeratio_total", "transients_goodmatches_total", "runtime", "runby" )
        attrs = {"th":{"bgcolor":"#EBEDEF"},}
        
        
class SumssNoMatchListTable(tables.Table):
    export_formats = ['csv',]
    match_id = tables.Column(verbose_name= 'ID')
    image_id = tables.LinkColumn('image_detail', args=[A('image_id'),], orderable=True, verbose_name= 'Img. ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "noaskapmatchtocatalog", A('match_id')], orderable=True, verbose_name= 'Source Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    cat_snr = tables.Column(verbose_name= 'Catalog SNR')
    scaled_cat_snr = tables.Column(verbose_name= 'Scaled ASKAP SNR')
    survey = tables.Column(verbose_name= 'Ref Survey')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    class Meta:
        model = Sumssnomatch
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("match_id", "image_id", "master_name", "ra", "dec", "cat_snr", "scaled_cat_snr", "survey", "pipelinetag", "usertag", "userreason", "checkedby")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}
    
class LargeRatioListTable(tables.Table):
    match_id = tables.Column(verbose_name= 'ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "largeratio", A('match_id')], orderable=True, verbose_name= 'Source Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    cat_snr = tables.Column(verbose_name= 'Catalog SNR')
    askap_cat_ratio = tables.Column(verbose_name= 'ASKAP / Catalog')
    survey = tables.Column(verbose_name= 'Ref Survey')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    class Meta:
        model = Largeratio
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("match_id", "master_name", "ra", "dec", "cat_snr", "askap_cat_ratio", "survey", "pipelinetag", "usertag", "userreason", "checkedby")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}   
        
class GoodMatchListTable(tables.Table):
    match_id = tables.Column(verbose_name= 'ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "goodmatch", A('match_id')], orderable=True, verbose_name= 'Source Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    survey = tables.Column(verbose_name= 'Survey')
    sumss_snr = tables.Column(verbose_name= 'SUMSS SNR')
    nvss_snr = tables.Column(verbose_name= 'NVSS SNR')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    class Meta:
        model = Goodmatch
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("match_id", "master_name", "ra", "dec", "survey", "sumss_snr", "nvss_snr", "pipelinetag", "usertag", "userreason", "checkedby")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}   
        
class AskapNotSeenListTable(tables.Table):
    match_id = tables.Column(verbose_name= 'ID')
    askap_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "nocatalogmatchtoaskap", A('match_id')], orderable=True, verbose_name= 'ASKAP Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    askap_snr = tables.Column(verbose_name= 'Local ASKAP SNR')
    scaled_cat_snr = tables.Column(verbose_name= 'Scaled Catalog SNR')
    survey = tables.Column(verbose_name= 'Ref Survey')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    class Meta:
        model = Askapnotseen
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("match_id", "askap_name", "ra", "dec", "askap_snr", "scaled_cat_snr", "survey", "pipelinetag", "usertag", "userreason", "checkedby")
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


