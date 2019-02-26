# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render

# Create your views here.
from django.http import HttpResponse

from .models import Image, Sumssnomatch, Largeratio, Askapnotseen, Goodmatch, Processingsettings

def home(request):
    images = Image.objects.all()
    return render(request, 'home.html', {'images': images})
    
def image_detail(request, pk):
    image = Image.objects.get(pk=pk)
    processing_options = Processingsettings.objects.get(image_id=pk)
    return render(request, 'image_detail.html', {'image':image, "processing_options":processing_options})
    
def sumssnomatch(request,pk):
    sumssnomatch_sources = Sumssnomatch.objects.all().filter(image_id=pk).order_by("match_id")
    image = Image.objects.get(pk=pk)
    return render(request, 'sumss_nomatch.html', {'sumssnomatch_sources':sumssnomatch_sources, 'image':image})
    
def largeratio(request,pk):
    largeratio_sources = Largeratio.objects.all().filter(image_id=pk).order_by("match_id")
    image = Image.objects.get(pk=pk)
    return render(request, 'large_ratio.html', {'largeratio_sources':largeratio_sources, 'image':image})
    
def askapnotseen(request,pk):
    askapnotseen_sources = Askapnotseen.objects.all().filter(image_id=pk).order_by("match_id")
    image = Image.objects.get(pk=pk)
    return render(request, 'askap_notseen.html', {'askapnotseen_sources':askapnotseen_sources, 'image':image})
    
def goodmatch(request,pk):
    # goodmatch_sources = Goodmatch.objects.all()
    goodmatch_sources = Goodmatch.objects.all().filter(image_id=pk).order_by("match_id")
    image = Image.objects.get(pk=pk)
    return render(request, 'good_match.html', {'goodmatch_sources':goodmatch_sources, 'image':image})
    
def crossmatch_detail(request,pk,querytype,cross_id):
    object_from_query={"sumssnomatch":Sumssnomatch,
                        "largeratio":Largeratio,
                        "askapnotseen":Askapnotseen,
                        "goodmatch":Goodmatch}
    title ={"sumssnomatch":"No ASKAP Match to SUMSS",
                        "largeratio":"Large Ratio",
                        "askapnotseen":"No SUMSS Match to ASKAP",
                        "goodmatch":"Good Matches"}
    html ={"sumssnomatch":"sumssnomatch",
                        "largeratio":"largeratio",
                        "askapnotseen":"askapnotseen",
                        "goodmatch":"goodmatch"}
    allsources = object_from_query[querytype].objects.all().filter(image_id=pk)
    min_id = min(list(allsources.values_list('match_id', flat=True)))
    max_id = max(list(allsources.values_list('match_id', flat=True)))
    crossmatch_source = object_from_query[querytype].objects.get(image_id=pk, match_id=cross_id)
    image = Image.objects.get(pk=pk)
    title_to_use = title[querytype]
    url_to_use = html[querytype]
    return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':querytype, 
                                                        'title':title_to_use, 'type_url':url_to_use, 'max_id':max_id, 'min_id':min_id},)
                                                        
def crossmatch_commit(request,pk,querytype,cross_id):
    username = request.GET['username']
    usertag = request.GET['usertag']
    object_from_query={"sumssnomatch":Sumssnomatch,
                        "largeratio":Largeratio,
                        "askapnotseen":Askapnotseen,
                        "goodmatch":Goodmatch}
    title ={"sumssnomatch":"No ASKAP Match to SUMSS",
                        "largeratio":"Large Ratio",
                        "askapnotseen":"No SUMSS Match to ASKAP",
                        "goodmatch":"Good Matches"}
    html ={"sumssnomatch":"sumssnomatch",
                        "largeratio":"largeratio",
                        "askapnotseen":"askapnotseen",
                        "goodmatch":"goodmatch"}
    allsources = object_from_query[querytype].objects.all().filter(image_id=pk)
    min_id = min(list(allsources.values_list('match_id', flat=True)))
    max_id = max(list(allsources.values_list('match_id', flat=True)))
    crossmatch_source = object_from_query[querytype].objects.get(image_id=pk, match_id=cross_id)
    if crossmatch_source.checkedby.lower()=="n/a":
        crossmatch_source.checkedby=username
        crossmatch_source.usertag=usertag
        crossmatch_source.save()
        saved=True
    else:
        saved=False
    image = Image.objects.get(pk=pk)
    title_to_use = title[querytype]
    url_to_use = html[querytype]
    return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':querytype, 
                                                        'title':title_to_use, 'type_url':url_to_use, 'max_id':max_id, 'min_id':min_id, "saved":saved},)
    