# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render
from django.shortcuts import redirect
from django_tables2 import RequestConfig
# from django_tables2.config import RequestConfig
from django_tables2.export.export import TableExport

# Create your views here.
from django.http import HttpResponse

from .models import Image, Sumssnomatch, Largeratio, Askapnotseen, Goodmatch, Processingsettings, Query
from .tables import ImageTable, SumssNoMatchListTable, LargeRatioListTable, GoodMatchListTable, AskapNotSeenListTable
from .forms import TagForm

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
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('noaskapmatchtosumss_image{}.{}'.format(pk, export_format))
    return render(request, 'sumss_nomatch.html', {'sumssnomatch_sources':sumssnomatch_sources, 'image':image, "table":table})
    
def largeratio(request,pk):
    largeratio_sources = Largeratio.objects.all().filter(image_id=pk).order_by("match_id")
    image = Image.objects.get(pk=pk)
    table = LargeRatioListTable(largeratio_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('largeratio_image{}.{}'.format(pk, export_format))
    return render(request, 'large_ratio.html', {'largeratio_sources':largeratio_sources, 'image':image, "table":table})
    
def askapnotseen(request,pk):
    askapnotseen_sources = Askapnotseen.objects.all().filter(image_id=pk).order_by("match_id")
    image = Image.objects.get(pk=pk)
    table = AskapNotSeenListTable(askapnotseen_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('nosumssmatchtoaskap_image{}.{}'.format(pk, export_format))
    return render(request, 'askap_notseen.html', {'askapnotseen_sources':askapnotseen_sources, 'image':image, "table":table})
    
def goodmatch(request,pk):
    # goodmatch_sources = Goodmatch.objects.all()
    goodmatch_sources = Goodmatch.objects.all().filter(image_id=pk).order_by("match_id")
    image = Image.objects.get(pk=pk)
    table = GoodMatchListTable(goodmatch_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('goodmatch_image{}.{}'.format(pk, export_format))
    return render(request, 'good_match.html', {'goodmatch_sources':goodmatch_sources, 'image':image, "table":table})
    
def query_queries(request):
    if request.method == "POST":
        form = TagForm(request.POST or None)
        if form.is_valid():
            data={
            "transient_type":form.cleaned_data["transient_type"],
            "user_tag":form.cleaned_data["user_tag"],
            "user":form.cleaned_data["user"]
            }
            # return render(request, 'search_results', context)
            return redirect("search_results", transient_type=data["transient_type"], user_tag=data["user_tag"], user=data["user"])
    else:
        form = TagForm(initial={'user': 'all'})
        context={"form":form}

    return render(request, 'search.html', context)
    
def search_results(request, transient_type, user_tag, user):
    if transient_type == "sumssnomatch":
        sumssnomatch_sources = Sumssnomatch.objects.all().filter(usertag=user_tag).order_by("match_id")
        if user != "all":
            sumssnomatch_sources = sumssnomatch_sources.filter(checkedby=user)
        table = SumssNoMatchListTable(sumssnomatch_sources)
        RequestConfig(request, paginate={'per_page': 100}).configure(table)
        friendly_type = "No ASKAP Match to SUMSS"
    elif transient_type == "askapnotseen":
        askapnotseen_sources = Askapnotseen.objects.all().filter(usertag=user_tag).order_by("match_id")
        if user != "all":
            askapnotseen_sources = askapnotseen_sources.filter(checkedby=user)
        table = AskapNotSeenListTable(askapnotseen_sources)
        RequestConfig(request, paginate={'per_page': 100}).configure(table)
        friendly_type = "No SUMSS Match to ASKAP"
    else:
        largeratio_sources = Largeratio.objects.all().filter(usertag=user_tag).order_by("match_id")
        if user != "all":
            largeratio_sources = largeratio_sources.filter(checkedby=user)
        table = LargeRatioListTable(largeratio_sources)
        RequestConfig(request, paginate={'per_page': 100}).configure(table)
        friendly_type = "Large Flux Ratio"
        
        
    # context = {
    #     'queries': queries
    # }
    # a=queries
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('search_results_{}_{}_{}.{}'.format(friendly_type.lower().replace(" ", "-"), user_tag, user, export_format))
        
    return render(request, "search_results.html", {"transient_type":friendly_type, "user_tag":user_tag, "table":table, "user":user})

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
    total=max_id-min_id+1
    crossmatch_source = object_from_query[querytype].objects.get(image_id=pk, match_id=cross_id)
    image = Image.objects.get(pk=pk)
    title_to_use = title[querytype]
    url_to_use = html[querytype]
    simbad_query="http://simbad.u-strasbg.fr/simbad/sim-coo?CooEpoch=2000&Coord={}d{}d&Radius.unit=arcmin&CooEqui=2000&CooFrame=FK5&Radius=10".format(crossmatch_source.ra, crossmatch_source.dec)
    return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':querytype, 
                                                        'title':title_to_use, 'type_url':url_to_use, 'max_id':max_id, 'min_id':min_id, 'total':total, "saved":False, "updated":False, "conflict":False,
                                                    'simbad_query':simbad_query},)
                                                        
def crossmatch_commit(request,pk,querytype,cross_id):
    user = request.user
    username = user.get_username()
    # username = request.GET['username']
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
    total=max_id-min_id+1
    crossmatch_source = object_from_query[querytype].objects.get(image_id=pk, match_id=cross_id)
    if crossmatch_source.checkedby.lower()=="n/a":
        crossmatch_source.checkedby=username
        crossmatch_source.usertag=usertag
        crossmatch_source.userreason=userreason
        crossmatch_source.save()
        saved=True
        updated=False
        conflict=False
    elif crossmatch_source.checkedby == username:
        crossmatch_source.checkedby=username
        crossmatch_source.usertag=usertag
        crossmatch_source.userreason=userreason
        crossmatch_source.save()
        updated=True
        saved=False
        conflict = False
    else:
        if user.is_staff:
            crossmatch_source.checkedby=username
            crossmatch_source.usertag=usertag
            crossmatch_source.userreason=userreason
            crossmatch_source.save()
            updated=True
            saved=False
            conflict = False
        else:
            saved=False
            updated=False
            conflict=True
    image = Image.objects.get(pk=pk)
    title_to_use = title[querytype]
    url_to_use = html[querytype]
    return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':querytype, 
                                                        'title':title_to_use, 'type_url':url_to_use, 'max_id':max_id, 'min_id':min_id, 'total':total, "saved":saved, "updated":updated, "conflict":conflict},)

