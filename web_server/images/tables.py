#!/usr/bin/env python

import django_tables2 as tables
from .models import Image, Crossmatches
from django_tables2.utils import A
from images.templatetags import units

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
    def value(self,value):
        return value

class FloatCoordColumn(tables.Column):
    def render(self, value):
        return "{:.3f}".format(value)
    def value(self,value):
        return value

class CapitalColumn(tables.Column):
    def render(self, value):
        return value.upper()


class ImageTable(tables.Table):
    export_formats = ['csv',]
    id = tables.Column(verbose_name= 'ID')
    name = tables.LinkColumn('image_detail', args=[A('pk')], orderable=True,)
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    # rms = RMSColumn(verbose_name="RMS (mJy)")
    matched_to = tables.Column(verbose_name= 'Matched To')
    number_askap_sources = tables.Column(verbose_name= '# ASKAP Sources')
    number_sumss_sources = tables.Column(verbose_name= '# SUMSS Sources')
    number_nvss_sources = tables.Column(verbose_name= '# NVSS Sources')
    transients_master_candidates_total = tables.Column(verbose_name= '# Candidate Sources')
    # transients_noaskapmatchtocatalog_total = tables.Column(verbose_name= '# No ASKAP match to Catalog')
    # transients_nocatalogmatchtoaskap_total = tables.Column(verbose_name= '# No Catalog match to ASKAP')
    # transients_largeratio_total = tables.Column(verbose_name= '# Large Ratio')
    # transients_goodmatches_total = tables.Column(verbose_name= '# Good Matches')
    runby = tables.Column(verbose_name= 'Run By')
    # claimed_by = tables.Column(verbose_name= 'Claimed By')
    class Meta:
        model = Image
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id", "name", "description", "ra", "dec", "matched_to", "number_askap_sources", "number_sumss_sources",
         "number_nvss_sources", "transients_master_candidates_total", "runtime", "runby")
        attrs = {"th":{"bgcolor":"#EBEDEF"},}
        row_attrs = {
                'complete': lambda record: 'true' if (record.number_candidates_checked >= record.transients_master_candidates_total) else 'false'
                }


class CrossmatchDetailTable(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    image_id = tables.LinkColumn('image_detail', args=[A('image_id'),], orderable=True, verbose_name= 'Img. ID')
    # master_name = tables.Column(verbose_name = 'Name')
    ra = RAColumn(verbose_name="RA")
    dec = DecColumn(verbose_name="Dec")
    gal_l = FloatColumn(verbose_name="Gal. l (deg)")
    gal_b = FloatColumn(verbose_name="Gal. b (deg)")
    d2d_askap_centre = FloatColumn(verbose_name = "Distance from ASKAP Centre (deg)")
    d2d_nn_askap_cat = FloatColumn(verbose_name = "Distance from ASKAP neighbour (arcsec)")
    survey = CapitalColumn(verbose_name= 'Survey Used')
    catalog_mosaic = tables.Column(verbose_name = 'Mosaic')
    transient_type = tables.Column(verbose_name = 'Candidate Type')

    class Meta:
        model = Crossmatches
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","image_id", "ra", "dec", "gal_l", "gal_b", "d2d_askap_centre","d2d_nn_askap_cat", "survey", "catalog_mosaic", "transient_type")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}

class CrossmatchDetailFluxTable(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    askap_scale_flux = RMSColumn(verbose_name= 'Scaled ASKAP Int. Flux (mJy)')
    askap_non_conv_flux = RMSColumn(verbose_name= 'Non-convolved Int. Flux (mJy)')
    # measured_askap_local_rms = RMSColumn(verbose_name= 'Local RMS (mJy)')
    # measured_askap_local_rms_2 = RMSColumn(verbose_name= 'Non-convolved local RMS (mJy)')
    # catalog_iflux = RMSColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    ratio_askap_flux = RMSColumn(verbose_name= 'Ratio ASKAP Int. Flux (mJy)')
    ratio_askap_flux_err = RMSColumn(verbose_name= 'Ratio ASKAP Int. Flux Error (mJy)')
    ratio_catalog_flux = RMSColumn(verbose_name= 'Ratio Cat. Int. Flux (mJy)')
    ratio_catalog_flux_err = RMSColumn(verbose_name= 'Ratio Cat. Int. Flux Error (mJy)')
    # catalog_scale_flux = FloatColumn(verbose_name= 'Scaled Cat. Int. Flux (mJy)')
    ratio = FloatColumn(verbose_name= 'Int. Flux Ratio')
    ratio_e = FloatColumn(verbose_name= 'Int. Flux Ratio Error')
    # askap_non_conv_scaled_flux = RMSColumn(verbose_name= 'Scaled Non-convolved Int. Flux (mJy)')
    # askap_non_conv_d2d = FloatColumn(verbose_name= 'Distance to ASKAP Convolved Source (arcsec)')
    # survey = CapitalColumn(verbose_name= 'Survey Used')
    aegean_rms_used = tables.Column(verbose_name = "3x RMS Used")
    inflated_convolved_flux = tables.Column(verbose_name = "Conv. Ratio Error")
    vs = FloatColumn(verbose_name= 'Vs')
    m = FloatColumn(verbose_name= '|m|')
    using_preconv = tables.Column(verbose_name = "Using Non-convolved info")

    class Meta:
        model = Crossmatches
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","aegean_rms_used", "inflated_convolved_flux", "using_preconv", "askap_iflux", "askap_scale_flux","askap_non_conv_flux",
        "ratio_askap_flux", "ratio_askap_flux_err", "ratio_catalog_flux", "ratio_catalog_flux_err", "ratio", "ratio_e",
        "vs", "m")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}


class NearestSourceDetailFluxTable(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "crossmatches", A('id')], orderable=True, verbose_name= 'Name')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    askap_scale_flux = RMSColumn(verbose_name= 'Scaled ASKAP Int. Flux (mJy)')
    catalog_iflux = RMSColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    ratio = FloatColumn(verbose_name= 'ASKAP / Cat Int. Flux Ratio')
    survey = CapitalColumn(verbose_name= 'Survey Used')

    class Meta:
        model = Crossmatches
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","master_name", "askap_iflux", "askap_scale_flux", "catalog_iflux", "ratio", "survey")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}


class AssocFluxTable(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    image_id = tables.LinkColumn('image_detail', args=[A('image_id'),], orderable=True, verbose_name= 'Img. ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "crossmatches", A('id')], orderable=True, verbose_name= 'Name')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    askap_scale_flux = RMSColumn(verbose_name= 'Scaled ASKAP Int. Flux (mJy)')
    catalog_iflux = RMSColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    ratio = FloatColumn(verbose_name= 'ASKAP / Cat Int. Flux Ratio')
    d2d_askap_centre = FloatColumn(verbose_name = "Distance from ASKAP Centre (deg)")
    survey = CapitalColumn(verbose_name= 'Survey Used')

    class Meta:
        model = Crossmatches
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","image_id", "master_name", "askap_iflux", "askap_scale_flux", "catalog_iflux", "ratio","d2d_askap_centre", "survey")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}


class CrossmatchTable(tables.Table):
    export_formats = ['csv',]
    id = tables.Column(verbose_name= 'ID')
    image_id = tables.LinkColumn('image_detail', args=[A('image_id'),], orderable=True, verbose_name= 'Img. ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "crossmatches", A('id')], orderable=True, verbose_name= 'Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    ratio = FloatColumn(verbose_name='Freq. Scaled Ratio')
    vs = FloatColumn(verbose_name='Vs')
    m = FloatColumn(verbose_name='m')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    # askap_snr = tables.Column(verbose_name= 'Local ASKAP SNR')
    catalog_iflux = RMSColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    # scaled_askap_snr = tables.Column(verbose_name= 'Scaled Catalog SNR')
    d2d_askap_centre = FloatColumn(verbose_name="Dist. from ASKAP Centre (deg)")
    survey = CapitalColumn(verbose_name= 'Ref Survey')
    transient_type = tables.Column(verbose_name= 'Type')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    askap_image = tables.Column()
    catalog_mosaic = tables.Column()

    def before_render(self, request):
        self.columns.hide("d2d_askap_centre")
        self.columns.hide("askap_image")
        self.columns.hide("catalog_mosaic")

    class Meta:
        model = Crossmatches
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","image_id","master_name", "ra", "dec", "askap_iflux","catalog_iflux", "ratio", "vs", "m",
        "d2d_askap_centre", "survey", "transient_type", "pipelinetag", "usertag", "userreason", "checkedby")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}
        row_attrs = {
                'highlight': lambda record: 'true' if (record.ratio >= 2.0 and record.pipelinetag == "Candidate" and record.checkedby != "N/A") else 'false'
                }

class CrossmatchTableAll(tables.Table):
    id = tables.Column(verbose_name= 'ID')
    image_id = tables.LinkColumn('image_detail', args=[A('image_id'),], orderable=True, verbose_name= 'Img. ID')
    master_name = tables.LinkColumn('crossmatch_detail', args=[A('image_id'), "crossmatches", A('id')], orderable=True, verbose_name= 'Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    ratio = FloatColumn(verbose_name='Freq. Scaled Ratio')
    vs = FloatColumn(verbose_name='Vs')
    m = FloatColumn(verbose_name='m')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    # askap_snr = tables.Column(verbose_name= 'Local ASKAP SNR')
    catalog_iflux = RMSColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    # scaled_askap_snr = tables.Column(verbose_name= 'Scaled Catalog SNR')
    d2d_askap_centre = FloatColumn(verbose_name="Dist. from ASKAP Centre (deg)")
    survey = CapitalColumn(verbose_name= 'Ref Survey')
    transient_type = tables.Column(verbose_name= 'Type')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    askap_image = tables.Column()
    catalog_mosaic = tables.Column()

    def before_render(self, request):
        self.columns.hide("d2d_askap_centre")
        self.columns.hide("askap_image")
        self.columns.hide("catalog_mosaic")


    class Meta:
        model = Crossmatches
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","image_id","master_name", "ra", "dec", "askap_iflux","catalog_iflux", "ratio", "vs", "m", "d2d_askap_centre", "survey", "transient_type", "pipelinetag", "usertag", "userreason", "checkedby")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}
        row_attrs = {
                'highlight': lambda record: 'true' if (record.ratio >= 2.0 and record.pipelinetag == "Candidate") else 'false'
                }

class CrossmatchQueryTable(tables.Table):
    export_formats = ['csv',]
    id = tables.Column(verbose_name= 'ID')
    image_id = tables.LinkColumn('image_detail', args=[A('image_id'),], orderable=True, verbose_name= 'Img. ID')
    master_name = tables.LinkColumn('crossmatch_detail_query', args=[A('id')], orderable=True, verbose_name= 'Name')
    ra = RAColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'RA' )
    dec = DecColumn(attrs={"td":{"style":"white-space:nowrap;"}}, verbose_name= 'Dec')
    gal_l = FloatColumn(verbose_name="Gal. l (deg)")
    gal_b = FloatColumn(verbose_name="Gal. b (deg)")
    ratio = FloatColumn(verbose_name='Freq. Scaled Ratio')
    vs = FloatColumn(verbose_name='Vs')
    m = FloatColumn(verbose_name='m')
    askap_iflux = RMSColumn(verbose_name= 'ASKAP Int. Flux (mJy)')
    askap_non_conv_flux = RMSColumn(verbose_name= 'ASKAP Raw Int. Flux (mJy)')
    # askap_snr = tables.Column(verbose_name= 'Local ASKAP SNR')
    catalog_iflux = RMSColumn(verbose_name= 'Cat. Int. Flux (mJy)')
    # scaled_askap_snr = tables.Column(verbose_name= 'Scaled Catalog SNR')
    d2d_askap_centre = FloatColumn(verbose_name="Dist. from ASKAP Centre (deg)")
    survey = CapitalColumn(verbose_name= 'Ref Survey')
    transient_type = tables.Column(verbose_name= 'Type')
    pipelinetag = tables.Column(verbose_name= 'Pipeline Tag')
    usertag = tables.Column(verbose_name= 'User Tag')
    userreason = tables.Column(verbose_name= 'User Reason')
    checkedby = tables.Column(verbose_name= 'Checked By')
    askap_image = tables.Column()
    catalog_mosaic = tables.Column()

    def before_render(self, request):
        self.columns.hide("d2d_askap_centre")
        self.columns.hide("askap_image")
        self.columns.hide("catalog_mosaic")
        self.columns.hide("askap_non_conv_flux")

    class Meta:
        model = Crossmatches
        template_name = 'django_tables2/bootstrap4.html'
        fields = ("id","image_id","master_name", "ra", "dec", "gal_l", "gal_b",
        "askap_iflux","catalog_iflux", "ratio", "vs", "m", "d2d_askap_centre",
        "survey", "transient_type", "pipelinetag", "usertag", "userreason", "checkedby")
        attrs = {"th":{"bgcolor":"#EBEDEF"}}
        row_attrs = {
                'highlight': lambda record: 'true' if (record.checkedby != "N/A") else 'false'
                }