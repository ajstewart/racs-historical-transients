#!/usr/bin/env python

from django import forms
from .models import Query

from django.contrib.auth.models import User

users = list(User.objects.values())

users = [("all", "all")]+[(i["username"],i["username"]) for i in sorted(users)]

tag_options=[("help", "Help"),
("ignore", "Ignore"),
("investigate", "Investigate"),
("problem", "Problem")]

transient_options=[("sumssnomatch", "No ASKAP match to Catalog"),
("askapnotseen", "No Catalog match to ASKAP"),
("largeratio", "Large flux ratio"),
{"transients", "Transients"}]

class TagForm(forms.ModelForm):
    class Meta:
            model = Query
            fields = ["transient_type", "user_tag", "user"]
            
    def __init__(self, *args, **kwargs):
            super(TagForm, self).__init__(*args, **kwargs)
            self.fields['transient_type'] = forms.ChoiceField(choices=transient_options)

            self.fields['user_tag'] = forms.ChoiceField(choices=tag_options)
            
            self.fields['user'] = forms.ChoiceField(choices=users)
    
