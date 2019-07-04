# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render
from django.shortcuts import redirect
from django_tables2 import RequestConfig
# from django_tables2.config import RequestConfig
from django_tables2.export.export import TableExport

# Create your views here.
from django.http import HttpResponse

from .models import Image, Sumssnomatch, Largeratio, Askapnotseen, Goodmatch, Processingsettings, Query, Transients
from .tables import ImageTable, SumssNoMatchListTable, LargeRatioListTable, GoodMatchListTable, AskapNotSeenListTable, CrossmatchDetailFluxTable, NearestSourceDetailFluxTable, TransientTable
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
    sumssnomatch_sources = Sumssnomatch.objects.all().filter(image_id=pk).filter(pipelinetag="Candidate").order_by("id")
    image = Image.objects.get(pk=pk)
    table = SumssNoMatchListTable(sumssnomatch_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('noaskapmatchtocatalog_image{}.{}'.format(pk, export_format))
    return render(request, 'sumss_nomatch.html', {'sumssnomatch_sources':sumssnomatch_sources, 'image':image, "table":table, "querytype":"noaskapmatchtocatalog"})
    
def largeratio(request,pk):
    largeratio_sources = Largeratio.objects.all().filter(image_id=pk).filter(pipelinetag__contains="Match ").order_by("id")
    image = Image.objects.get(pk=pk)
    table = LargeRatioListTable(largeratio_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('largeratio_image{}.{}'.format(pk, export_format))
    return render(request, 'large_ratio.html', {'largeratio_sources':largeratio_sources, 'image':image, "table":table, "querytype":"largeratio"})
    
def askapnotseen(request,pk):
    askapnotseen_sources = Askapnotseen.objects.all().filter(image_id=pk).filter(pipelinetag="Candidate").order_by("id")
    image = Image.objects.get(pk=pk)
    table = AskapNotSeenListTable(askapnotseen_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('nocatalogmatchtoaskap_image{}.{}'.format(pk, export_format))
    return render(request, 'askap_notseen.html', {'askapnotseen_sources':askapnotseen_sources, 'image':image, "table":table, "querytype":"nocatalogmatchtoaskap"})
    
def goodmatch(request,pk):
    # goodmatch_sources = Goodmatch.objects.all()
    goodmatch_sources = Goodmatch.objects.all().filter(image_id=pk).order_by("id")
    image = Image.objects.get(pk=pk)
    table = GoodMatchListTable(goodmatch_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('goodmatch_image{}.{}'.format(pk, export_format))
    return render(request, 'good_match.html', {'goodmatch_sources':goodmatch_sources, 'image':image, "table":table})

def transients(request,pk):
    # goodmatch_sources = Goodmatch.objects.all()
    transient_sources = Transients.objects.all().filter(image_id=pk).order_by("id")
    image = Image.objects.get(pk=pk)
    table = TransientTable(transient_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('transients_image{}.{}'.format(pk, export_format))
    return render(request, 'transients.html', {'transient_sources':transient_sources, 'image':image, "table":table, "querytype":"transients"})
    
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
        sumssnomatch_sources = Sumssnomatch.objects.all().filter(usertag=user_tag).order_by("id")
        if user != "all":
            sumssnomatch_sources = sumssnomatch_sources.filter(checkedby=user)
        table = SumssNoMatchListTable(sumssnomatch_sources)
        RequestConfig(request, paginate={'per_page': 100}).configure(table)
        friendly_type = "No ASKAP Match to Catalog"
    elif transient_type == "askapnotseen":
        askapnotseen_sources = Askapnotseen.objects.all().filter(usertag=user_tag).order_by("id")
        if user != "all":
            askapnotseen_sources = askapnotseen_sources.filter(checkedby=user)
        table = AskapNotSeenListTable(askapnotseen_sources)
        RequestConfig(request, paginate={'per_page': 100}).configure(table)
        friendly_type = "No Catalog Match to ASKAP"
    elif transient_type == "largeratio":
        largeratio_sources = Largeratio.objects.all().filter(usertag=user_tag).order_by("id")
        if user != "all":
            largeratio_sources = largeratio_sources.filter(checkedby=user)
        table = LargeRatioListTable(largeratio_sources)
        RequestConfig(request, paginate={'per_page': 100}).configure(table)
        friendly_type = "Large Flux Ratio"
    else:
        transient_sources = Transients.objects.all().filter(usertag=user_tag).order_by("id")
        if user != "all":
            transient_sources = transient_sources.filter(checkedby=user)
        table = TransientTable(transient_sources)
        RequestConfig(request, paginate={'per_page': 100}).configure(table)
        friendly_type = "Transients"
        
        
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
    object_from_query={"noaskapmatchtocatalog":Sumssnomatch,
                        "largeratio":Largeratio,
                        "nocatalogmatchtoaskap":Askapnotseen,
                        "goodmatch":Goodmatch,
                        "transients":Transients}
    title ={"noaskapmatchtocatalog":"No ASKAP Match to Catalog",
                        "largeratio":"Large Ratio",
                        "nocatalogmatchtoaskap":"No Catalog Match to ASKAP",
                        "goodmatch":"Good Matches",
                        "transients":"Transients"}
    html ={"noaskapmatchtocatalog":"noaskapmatchtocatalog",
                        "largeratio":"largeratio",
                        "nocatalogmatchtoaskap":"nocatalogmatchtoaskap",
                        "goodmatch":"goodmatch",
                        "transients":"transients"}
    allsources = object_from_query[querytype].objects.all().filter(image_id=pk)
    if querytype == "nocatalogmatchtoaskap" or querytype == "noaskapmatchtocatalog":
        thesources = allsources.filter(pipelinetag="Candidate")
    elif querytype == "largeratio":
        thesources = allsources.filter(pipelinetag__contains="Match ")
    min_id = min(list(allsources.values_list('id', flat=True)))
    max_id = max(list(allsources.values_list('id', flat=True)))
    total=max_id-min_id+1
    crossmatch_source = object_from_query[querytype].objects.get(image_id=pk, id=cross_id)
    table = CrossmatchDetailFluxTable(allsources.filter(id=cross_id))
    image = Image.objects.get(pk=pk)
    if querytype != "goodmatch":
        #Also as large ratio we can fetch the nearest sources table
        nearest_sources = crossmatch_source.nearest_sources
        nearest_sources = nearest_sources.split(",")
        print nearest_sources
        nearest_sources = Transients.objects.all().filter(image_id=pk, master_name__in=nearest_sources)
        nearest_sources_table = NearestSourceDetailFluxTable(nearest_sources)
    else:
        nearest_sources_table = ""
    title_to_use = title[querytype]
    url_to_use = html[querytype]
    simbad_query="http://simbad.u-strasbg.fr/simbad/sim-coo?CooEpoch=2000&Coord={}d{}d&Radius.unit=arcmin&CooEqui=2000&CooFrame=FK5&Radius=10".format(crossmatch_source.ra, crossmatch_source.dec)
    ned_query="https://ned.ipac.caltech.edu/conesearch?search_type=Near%20Position%20Search&coordinates={}d%20{}d&radius=2.0&in_csys=Equatorial&in_equinox=J2000.0&out_csys=Equatorial&out_equinox=J2000.0&hconst=67.8&omegam=0.308&omegav=0.692&wmap=4&corr_z=1".format(crossmatch_source.ra, crossmatch_source.dec)
    return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':querytype, 
                                                        'title':title_to_use, 'type_url':url_to_use, 'max_id':max_id, 'min_id':min_id, 'total':total, "saved":False, "updated":False, "conflict":False,
                                                    'simbad_query':simbad_query, 'ned_query':ned_query, 'table':table, 'nearest_sources_table':nearest_sources_table},)
                                                        
def crossmatch_commit(request,pk,querytype,cross_id):
    user = request.user
    username = user.get_username()
    # username = request.GET['username']
    usertag = request.GET['usertag']
    userreason = request.GET['userreason']
    object_from_query={"noaskapmatchtocatalog":Sumssnomatch,
                        "largeratio":Largeratio,
                        "nocatalogmatchtoaskap":Askapnotseen,
                        "goodmatch":Goodmatch,
                        "transients":Transients}
    title ={"noaskapmatchtocatalog":"No ASKAP Match to Catalog",
                        "largeratio":"Large Ratio",
                        "nocatalogmatchtoaskap":"No Catalog Match to ASKAP",
                        "goodmatch":"Good Matches",
                        "transients":"Transients"}
    html ={"noaskapmatchtocatalog":"noaskapmatchtocatalog",
                        "largeratio":"largeratio",
                        "nocatalogmatchtoaskap":"nocatalogmatchtoaskap",
                        "goodmatch":"goodmatch",
                        "transients":"transients"}
    allsources = object_from_query[querytype].objects.all().filter(image_id=pk)
    if querytype == "nocatalogmatchtoaskap" or querytype == "noaskapmatchtocatalog":
        thesources = allsources.filter(pipelinetag="Candidate")
    elif querytype == "largeratio":
        thesources = allsources.filter(pipelinetag__contains="Match ")
    min_id = min(list(allsources.values_list('id', flat=True)))
    max_id = max(list(allsources.values_list('id', flat=True)))
    total=max_id-min_id+1
    crossmatch_source = object_from_query[querytype].objects.get(image_id=pk, id=cross_id)
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
    table = CrossmatchDetailFluxTable(allsources.filter(id=cross_id))
    if querytype != "goodmatch":
        #Also as large ratio we can fetch the nearest sources table
        nearest_sources = crossmatch_source.nearest_sources
        nearest_sources = nearest_sources.split(",")
        print nearest_sources
        nearest_sources = Goodmatch.objects.all().filter(image_id=pk, master_name__in=nearest_sources)
        nearest_sources_table = NearestSourceDetailFluxTable(nearest_sources)
    else:
        nearest_sources_table = ""
    simbad_query="http://simbad.u-strasbg.fr/simbad/sim-coo?CooEpoch=2000&Coord={}d{}d&Radius.unit=arcmin&CooEqui=2000&CooFrame=FK5&Radius=10".format(crossmatch_source.ra, crossmatch_source.dec)
    ned_query="https://ned.ipac.caltech.edu/conesearch?search_type=Near%20Position%20Search&coordinates={}d%20%2B{}d&radius=2.0&in_csys=Equatorial&in_equinox=J2000.0&out_csys=Equatorial&out_equinox=J2000.0&hconst=67.8&omegam=0.308&omegav=0.692&wmap=4&corr_z=1".format(crossmatch_source.ra, crossmatch_source.dec)
    return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':querytype, 
                                                        'title':title_to_use, 'type_url':url_to_use, 'max_id':max_id, 'min_id':min_id, 
                                                        'total':total, "saved":saved, "updated":updated, "conflict":conflict, 'simbad_query':simbad_query, 'ned_query':ned_query, 'table':table, 'nearest_sources_table':nearest_sources_table},)


def crossmatch_quickview(request,pk,querytype):
    object_from_query={"noaskapmatchtocatalog":Sumssnomatch,
                        "largeratio":Largeratio,
                        "nocatalogmatchtoaskap":Askapnotseen,
                        "goodmatch":Goodmatch,
                        "transients":Transients}
    title ={"noaskapmatchtocatalog":"No ASKAP Match to Catalog",
                        "largeratio":"Large Ratio",
                        "nocatalogmatchtoaskap":"No Catalog Match to ASKAP",
                        "goodmatch":"Good Matches",
                        "transients":"Transients"}
    html ={"noaskapmatchtocatalog":"noaskapmatchtocatalog",
                        "largeratio":"largeratio",
                        "nocatalogmatchtoaskap":"nocatalogmatchtoaskap",
                        "goodmatch":"goodmatch",
                        "transients":"transients"}
    if querytype=="transients":
        allsources = object_from_query[querytype].objects.all().filter(image_id=pk).filter(ratio__gte=2.0).order_by("id")
    else:
        allsources = object_from_query[querytype].objects.all().filter(image_id=pk).exclude(pipelinetag="Candidate")
    image = Image.objects.get(pk=pk)
    title_to_use = title[querytype]
    total = len(list(allsources.values_list('id', flat=True)))
    return render(request, 'crossmatch_quickview.html', {"crossmatch_sources":allsources,"image":image, "querytype":querytype, "title":title_to_use, "total":total, "html":html[querytype]})
    
