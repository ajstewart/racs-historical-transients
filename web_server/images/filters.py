from .models import Crossmatches
from .tables import CrossmatchTable
import django_filters
from django_filters.widgets import RangeWidget
from crispy_forms.bootstrap import InlineField
from django import forms
from django_filters.filters import Filter
from decimal import Decimal
from django.contrib.auth.models import User

class MyRangeWidget(RangeWidget):
    def __init__(self, from_attrs=None, to_attrs=None, attrs=None):
        super(MyRangeWidget, self).__init__(attrs)

        if from_attrs:
            self.widgets[0].attrs.update(from_attrs)
        if to_attrs:
            self.widgets[1].attrs.update(to_attrs)

class ConeSearchWidget(forms.widgets.MultiWidget):
    def __init__(self, attrs=None, dt=None, mode=0):

        _CHOICES = (
        ('arcsec', 'arcsec'),
        ('arcmin', 'arcmin'),
        ('deg', 'deg'),
        )


        _widgets = (
            forms.widgets.NumberInput(attrs={"placeholder":"RA (deg)", "step":0.0000000001, "style":"margin:5px", 'style':'max-width: 8em'}),
            forms.widgets.NumberInput(attrs={"placeholder":"Dec (deg)", "step":0.0000000001, "style":"margin:5px", 'style':'max-width: 8em'}),
            forms.widgets.NumberInput(attrs={"placeholder":"Radius", "step":0.0000000001, "style":"margin:5px", 'style':'max-width: 8em'}),
            forms.widgets.Select(attrs={"placeholder":"Radius", "style":"margin:5px"}, choices = _CHOICES)
            )
        super(ConeSearchWidget, self).__init__(_widgets, attrs)


class ConeSearchFilter(Filter):
    # field_class = forms.DecimalField

    def filter(self, qs, value):

        return super()

class MyRangeFilter(Filter):
    field_class = django_filters.fields.RangeField

    def filter(self, qs, value):
        if value:
            if value.start is not None and value.stop is not None:
                value_start = value.start/Decimal(1000.)
                value_stop = value.stop/Decimal(1000.)
                lookup = '%s__range' % self.field_name
                return self.get_method(qs)(**{lookup: (value_start, value_stop)})
            else:
                if value.start is not None:
                    value_start = value.start/Decimal(1000.)
                    qs = self.get_method(qs)(**{'%s__gte' % self.field_name: value_start})
                if value.stop is not None:
                    value_stop = value.stop/Decimal(1000.)
                    qs = self.get_method(qs)(**{'%s__lte' % self.field_name: value_stop})
            if self.distinct:
                qs = qs.distinct()
        return qs


class TransientFilter(django_filters.FilterSet):

    def __init__(self, *args, **kwargs):
        super(TransientFilter, self).__init__(*args, **kwargs)
        if self.data == {}:
            if not self.queryset.image_crossmatch:
                self.queryset = self.queryset.none()

    def cone_search_filter(self, queryset, name, value):
        #Don't run if not all value are entered.
        if '' in value or None in value:
            return queryset
        else:
            ra=float(value[0])
            dec=float(value[1])
            radius=float(value[2])
            unit = value[3]
            scale_factors = {"arcsec":3600., "arcmin":60.}
            if unit in scale_factors:
                radius/=scale_factors[unit]

            return queryset.extra(where=["q3c_radial_query(ra, dec, {:.6f}, {:.6f}, {:.6f})".format(ra,dec,radius)])

    def sort_by_filter(self, queryset, name, value):
        print(value)
        if value=='' or value==None:
            return queryset
        else:
            return queryset.order_by(value)

    source_id = django_filters.NumberFilter(field_name ="id", label="Source ID")
    image_id = django_filters.NumberFilter(field_name ="image_id", label="Image ID")

    d2d__gt = django_filters.RangeFilter(field_name = 'd2d_askap_centre', widget=MyRangeWidget(from_attrs={'placeholder': 'Min', 'style':'max-width: 8em'}, to_attrs={'placeholder':'Max', 'style':'max-width: 8em'}))
    ratio__gt = django_filters.RangeFilter(field_name = 'ratio', widget=MyRangeWidget(from_attrs={'placeholder': 'Min', 'style':'max-width: 8em'}, to_attrs={'placeholder':'Max', 'style':'max-width: 8em'}))

    vs__gt = django_filters.RangeFilter(field_name = 'vs', widget=MyRangeWidget(from_attrs={'placeholder': 'Min', 'style':'max-width: 8em'}, to_attrs={'placeholder':'Max', 'style':'max-width: 8em'}))
    m_abs__gt = django_filters.RangeFilter(field_name = 'm_abs', widget=MyRangeWidget(from_attrs={'placeholder': 'Min', 'style':'max-width: 8em'}, to_attrs={'placeholder':'Max', 'style':'max-width: 8em'}), label="m")

    askap_iflux__gt = MyRangeFilter(field_name = 'askap_iflux', widget=MyRangeWidget(from_attrs={'placeholder': 'Min (mJy)', 'style':'max-width: 8em'}, to_attrs={'placeholder':'Max (mJy)', 'style':'max-width: 8em'}), label="ASKAP Integrated Flux")
    catalog_iflux__gt = MyRangeFilter(field_name = 'catalog_iflux', widget=MyRangeWidget(from_attrs={'placeholder': 'Min (mJy)', 'style':'max-width: 8em'}, to_attrs={'placeholder':'Max (mJy)', 'style':'max-width: 8em'}), label="Catalogue Integrated Flux")

    pipelinetag_choices = (
        ('Good', 'Good'),
        ('Created convolved source', 'Created convolved source'),
        ('Edge of ASKAP image', 'Edge of ASKAP image'),
        ('Extended source', 'Extended source'),
        ('Failed', 'Failed'),
        ('Inflated convolved flux error', 'Inflated convolved flux error'),
        ('Large island source', 'Large island source'),
        ('Likely artefact (bright source)', 'Likely artefact (bright source)'),
        ('Likely bright extended structure', 'Likely bright extended structure'),
        ('Likely double', 'Likely double/multi source'),
        )

    pipelinetag = django_filters.ChoiceFilter(field_name = 'pipelinetag', choices=pipelinetag_choices)

    flagged = django_filters.BooleanFilter(field_name='flagged')

    convolved_error_choices = (
        ('True', 'True'),
        ('False', 'False'),
        )

    inflated_convolved_flux = django_filters.ChoiceFilter(field_name = 'inflated_convolved_flux', choices=convolved_error_choices)

    try:
        users = list(User.objects.values())
    except:
        users = []

    users = [(i["username"],i["username"]) for i in sorted(users, key=lambda users: users['username'])]+[("N/A", "Unchecked")]

    checkedby = django_filters.ChoiceFilter(field_name = 'checkedby', choices=users)

    transient_type_choices = (
        ('Match', 'Match'),
        ('No askap match', 'No ASKAP match'),
        ('No catalog match', 'No catalogue match'),
        )

    transient_type = django_filters.MultipleChoiceFilter(choices=transient_type_choices, widget=forms.CheckboxSelectMultiple)

    survey_choices = (
        ('sumss', 'SUMSS'),
        ('nvss', 'NVSS'),
        )

    survey = django_filters.MultipleChoiceFilter(choices=survey_choices, widget=forms.CheckboxSelectMultiple)

    usertag_choices = (
        ('investigate', 'Investigate'),
        ('transient', 'Transient'),
        ('help', 'Help'),
        ('problem', 'Problem'),
        ('ignore', 'Ignore'),
        )

    usertag = django_filters.ChoiceFilter(field_name = 'usertag', choices=usertag_choices)

    userreason_choices = (
        ('askap artefact', 'ASKAP artefact'),
        ('bright source', 'Bright source'),
        ('catalogue artefact', 'Catalogue artefact'),
        ('convolved flux error', 'Convolved flux error'),
        ('edge of field', 'Edge of field'),
        ('extended', 'Extended'),
        ('is match', 'Is match'),
        ('investigate', 'Investigate'),
        ('multiple', 'Multiple'),
        ('non-detection', 'Non-detection'),
        ('not sure', 'Not sure'),
        ('src. ext. error', 'Source extraction error'),
        ('too noisy', 'Too noisy'),
        )

    userreason = django_filters.ChoiceFilter(field_name = 'userreason', choices=userreason_choices)

    ra__gt = django_filters.RangeFilter(field_name = 'ra', widget=MyRangeWidget(from_attrs={'placeholder': 'Min (deg)', 'style':'max-width: 8em'}, to_attrs={'placeholder':'Max (deg)', 'style':'max-width: 8em'}), label="Right Ascension")
    dec__gt = django_filters.RangeFilter(field_name = 'dec', widget=MyRangeWidget(from_attrs={'placeholder': 'Min (deg)', 'style':'max-width: 8em'}, to_attrs={'placeholder':'Max (deg)', 'style':'max-width: 8em'}), label="Declination")

    gal_l__gt = django_filters.RangeFilter(field_name = 'gal_l', widget=MyRangeWidget(from_attrs={'placeholder': 'Min (deg)', 'style':'max-width: 8em'}, to_attrs={'placeholder':'Max (deg)', 'style':'max-width: 8em'}), label="Galactic l")
    gal_b__gt = django_filters.RangeFilter(field_name = 'gal_b', widget=MyRangeWidget(from_attrs={'placeholder': 'Min (deg)', 'style':'max-width: 8em'}, to_attrs={'placeholder':'Max (deg)', 'style':'max-width: 8em'}), label="Galactic b")

    d2d_nn_askap_cat__gt = django_filters.RangeFilter(field_name = 'd2d_nn_askap_cat',
        widget=MyRangeWidget(from_attrs={'placeholder': 'Min (arcsec)', 'style':'max-width: 8em'}, to_attrs={'placeholder':'Max (arcsec)', 'style':'max-width: 8em'}), label="Distance to nearest ASKAP neighbour.")

    cone_search = ConeSearchFilter(field_name = 'cone_search', widget=ConeSearchWidget(), method="cone_search_filter", label="Cone Search")

    sortby_choices = (
    ('', '---------'),
    ('ratio', 'Ratio (asc)'),
    ('-ratio', 'Ratio (desc)'),
    ('askap_iflux', 'ASKAP Int. Flux (asc)'),
    ('-askap_iflux', 'ASKAP Int. Flux (desc)'),
    ('catalog_iflux', 'Catalogue Int. Flux (asc)'),
    ('-catalog_iflux', 'Catalogue Int. Flux (desc)'),
    ('id', "ID")
    )

    sort_by = ConeSearchFilter(field_name = 'sort_by', widget=forms.widgets.Select(choices=sortby_choices), method="sort_by_filter", label="Sort By")

    class Meta:
        model = Crossmatches
        fields = []