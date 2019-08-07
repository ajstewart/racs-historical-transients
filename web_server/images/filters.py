from .models import Transients
from .tables import TransientTable
import django_filters
from django_filters.widgets import RangeWidget
from crispy_forms.bootstrap import InlineField
from django import forms

class MyRangeWidget(RangeWidget):
    def __init__(self, from_attrs=None, to_attrs=None, attrs=None):
        super(MyRangeWidget, self).__init__(attrs)

        if from_attrs:
            self.widgets[0].attrs.update(from_attrs)
        if to_attrs:
            self.widgets[1].attrs.update(to_attrs)


class TransientFilter(django_filters.FilterSet):
    
    def __init__(self, *args, **kwargs):
        super(TransientFilter, self).__init__(*args, **kwargs)
        if self.data == {}:
            self.queryset = self.queryset.none()
        # self.form.fields["transient_type"].widget.attrs['class'] == 
    
    d2d__gt = django_filters.RangeFilter(name = 'd2d_askap_centre', widget=MyRangeWidget(from_attrs={'placeholder': 'Min'}, to_attrs={'placeholder':'Max'}))
    ratio__gt = django_filters.RangeFilter(name = 'ratio', widget=MyRangeWidget(from_attrs={'placeholder': 'Min'}, to_attrs={'placeholder':'Max'}))
    
    pipelinetag_choices = (
        ('Candidate', 'Candidate'),
        )
    
    pipelinetag = django_filters.ChoiceFilter(name = 'pipelinetag', choices=pipelinetag_choices)
    
    convolved_error_choices = (
        ('True', 'True'),
        ('False', 'False'),
        )
    
    inflated_convolved_flux = django_filters.ChoiceFilter(name = 'inflated_convolved_flux', choices=convolved_error_choices)
    
    
    transient_type_choices = (
        ('Good match', 'Good match'),
        ('No askap match', 'No ASKAP match'),
        ('No catalog match', 'No catalogue match'),
        )
    
    transient_type = django_filters.MultipleChoiceFilter(choices=transient_type_choices, widget=forms.CheckboxSelectMultiple)
    
    survey_choices = (
        ('SUMSS', 'SUMSS'),
        ('NVSS', 'NVSS'),
        )
    
    survey = django_filters.MultipleChoiceFilter(choices=survey_choices, widget=forms.CheckboxSelectMultiple)
    
    usertag_choices = (
        ('investigate', 'Investigate'),
        )
    
    usertag = django_filters.ChoiceFilter(name = 'usertag', choices=usertag_choices)
    
    class Meta:
        model = Transients
        fields = []
        