#!/usr/bin/env python

from django import forms
from .models import Query

from django.contrib.auth.models import User

try:
    users = list(User.objects.values())
except:
    users = []

users = [("all", "all")]+[(i["username"],i["username"]) for i in sorted(users, key=lambda users: users['username'])]

tag_options=[("investigate", "Investigate"),
("help", "Help"),
("ignore", "Ignore"),
("problem", "Problem")]

transient_options=[
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

