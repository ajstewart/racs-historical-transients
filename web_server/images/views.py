# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render
from django_tables2 import RequestConfig

# Create your views here.
from django.http import HttpResponse

from .models import Image, Sumssnomatch, Largeratio, Askapnotseen, Goodmatch, Processingsettings
from .tables import ImageTable, SumssNoMatchListTable, LargeRatioListTable, GoodMatchListTable, AskapNotSeenListTable

def home(request):
    images = Image.objects.all()
    table = ImageTable(images)
    RequestConfig(request).configure(table)
    return render(request, 'home.html', {'images': images, 'table':table})
    
def image_detail(request, pk):
    image = Image.objects.get(pk=pk)
    processing_options = Processingsettings.objects.get(image_id=pk)
    return render(request, 'image_detail.html', {'image':image, "processing_options":processing_options})
    
def sumssnomatch(request,pk):
    sumssnomatch_sources = Sumssnomatch.objects.all().filter(image_id=pk).order_by("match_id")
    image = Image.objects.get(pk=pk)
    table = SumssNoMatchListTable(sumssnomatch_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    return render(request, 'sumss_nomatch.html', {'sumssnomatch_sources':sumssnomatch_sources, 'image':image, "table":table})
    
def largeratio(request,pk):
    largeratio_sources = Largeratio.objects.all().filter(image_id=pk).order_by("match_id")
    image = Image.objects.get(pk=pk)
    table = LargeRatioListTable(largeratio_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    return render(request, 'large_ratio.html', {'largeratio_sources':largeratio_sources, 'image':image, "table":table})
    
def askapnotseen(request,pk):
    askapnotseen_sources = Askapnotseen.objects.all().filter(image_id=pk).order_by("match_id")
    image = Image.objects.get(pk=pk)
    table = AskapNotSeenListTable(askapnotseen_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    return render(request, 'askap_notseen.html', {'askapnotseen_sources':askapnotseen_sources, 'image':image, "table":table})
    
def goodmatch(request,pk):
    # goodmatch_sources = Goodmatch.objects.all()
    goodmatch_sources = Goodmatch.objects.all().filter(image_id=pk).order_by("match_id")
    image = Image.objects.get(pk=pk)
    table = GoodMatchListTable(goodmatch_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    return render(request, 'good_match.html', {'goodmatch_sources':goodmatch_sources, 'image':image, "table":table})
    
def crossmatch_detail(request,pk,querytype,cross_id):
    object_from_query={"noaskapmatchtosumss":Sumssnomatch,
                        "largeratio":Largeratio,
                        "nosumssmatchtoaskap":Askapnotseen,
                        "goodmatch":Goodmatch}
    title ={"noaskapmatchtosumss":"No ASKAP Match to SUMSS",
                        "largeratio":"Large Ratio",
                        "nosumssmatchtoaskap":"No SUMSS Match to ASKAP",
                        "goodmatch":"Good Matches"}
    html ={"noaskapmatchtosumss":"noaskapmatchtosumss",
                        "largeratio":"largeratio",
                        "nosumssmatchtoaskap":"nosumssmatchtoaskap",
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
    userreason = request.GET['userreason']
    object_from_query={"noaskapmatchtosumss":Sumssnomatch,
                        "largeratio":Largeratio,
                        "nosumssmatchtoaskap":Askapnotseen,
                        "goodmatch":Goodmatch}
    title ={"noaskapmatchtosumss":"No ASKAP Match to SUMSS",
                        "largeratio":"Large Ratio",
                        "nosumssmatchtoaskap":"No SUMSS Match to ASKAP",
                        "goodmatch":"Good Matches"}
    html ={"noaskapmatchtosumss":"noaskapmatchtosumss",
                        "largeratio":"largeratio",
                        "nosumssmatchtoaskap":"nosumssmatchtoaskap",
                        "goodmatch":"goodmatch"}
    allsources = object_from_query[querytype].objects.all().filter(image_id=pk)
    min_id = min(list(allsources.values_list('match_id', flat=True)))
    max_id = max(list(allsources.values_list('match_id', flat=True)))
    crossmatch_source = object_from_query[querytype].objects.get(image_id=pk, match_id=cross_id)
    if crossmatch_source.checkedby.lower()=="n/a":
        crossmatch_source.checkedby=username
        crossmatch_source.usertag=usertag
        crossmatch_source.userreason=userreason
        crossmatch_source.save()
        saved=True
        updated=False
    elif crossmatch_source.checkedby == username:
        crossmatch_source.checkedby=username
        crossmatch_source.usertag=usertag
        crossmatch_source.userreason=userreason
        crossmatch_source.save()
        updated=True
        saved=False
    else:
        saved=False
        updated=False
    image = Image.objects.get(pk=pk)
    title_to_use = title[querytype]
    url_to_use = html[querytype]
    return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':querytype, 
                                                        'title':title_to_use, 'type_url':url_to_use, 'max_id':max_id, 'min_id':min_id, "saved":saved, "updated":updated},)
    