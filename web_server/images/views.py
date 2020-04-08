# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib.auth.models import User
from django.shortcuts import render
from django.shortcuts import redirect, reverse
from django_tables2 import RequestConfig
# from django_tables2.config import RequestConfig
from django_tables2.export.export import TableExport

# Create your views here.
from django.http import HttpResponse

from .models import Image, Processingsettings, Query, Crossmatches
from .tables import ImageTable, CrossmatchDetailFluxTable, NearestSourceDetailFluxTable, CrossmatchTable
from .tables import ImageCrossmatchTable, CrossmatchTableAll, CrossmatchDetailTable, CrossmatchQueryTable, AssocFluxTable
from .forms import TagForm
from .filters import TransientFilter

import os
from django.conf import settings
import slack


def home(request):
    images = Image.objects.all().order_by("id")
    table = ImageTable(images)
    RequestConfig(request).configure(table)
    export_format = request.GET.get('_export', None)
    vs_thresh = settings.CANDIDATE_VS_THRESHOLD
    m_thresh = settings.CANDIDATE_M_THRESHOLD

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('images.{}'.format(export_format))

    try:
        hotkeys_off = request.GET["hotkeys_off"]
        request.session["hotkeys"]="false"
        next = request.GET["next"]
        return redirect(next)
    except:
        pass

    return render(request, 'home.html', {'images': images, 'table':table, 'vs_thresh':vs_thresh, 'm_thresh':m_thresh})

def image_detail(request, pk):
    image = Image.objects.get(pk=pk)
    processing_options = Processingsettings.objects.get(image_id=pk)
    vs_thresh = settings.CANDIDATE_VS_THRESHOLD
    m_thresh = settings.CANDIDATE_M_THRESHOLD
    num = len(Crossmatches.objects.all().filter(
        image_id=pk
    ).filter(pipelinetag="Good").filter(
        vs__gte=vs_thresh
    ).filter(m_abs__gte=m_thresh))
    flagged_num = len(Crossmatches.objects.all().filter(
        image_id=pk
    ).exclude(pipelinetag="Good").filter(
        vs__gte=vs_thresh
    ).filter(m_abs__gte=m_thresh))
    return render(request, 'image_detail.html', {'image':image, "processing_options":processing_options, "claimed_success":False, "claim_attempt":False,
        'vs_thresh':vs_thresh, 'm_thresh':m_thresh, 'num_candidates':num, 'flagged_num_candidates':flagged_num})

def image_detail_claim(request, pk):
    user = request.user
    username = user.get_username()
    image = Image.objects.get(pk=pk)
    processing_options = Processingsettings.objects.get(image_id=pk)
    vs_thresh = settings.CANDIDATE_VS_THRESHOLD
    m_thresh = settings.CANDIDATE_M_THRESHOLD
    num = len(Crossmatches.objects.all().filter(
        image_id=pk
    ).filter(pipelinetag="Good").filter(
        vs__gte=vs_thresh
    ).filter(m_abs__gte=m_thresh))
    flagged_num = len(Crossmatches.objects.all().filter(
        image_id=pk
    ).exclude(pipelinetag="Good").filter(
        vs__gte=vs_thresh
    ).filter(m_abs__gte=m_thresh))
    # username = request.GET['username']
    if image.claimed_by == "Unclaimed":
        image.claimed_by = username
        image.save()
        claimed_success = True
    else:
        claimed_success = False

    return render(request, 'image_detail.html', {'image':image, "processing_options":processing_options, "claimed_success":claimed_success, "claim_attempt":True,
        'vs_thresh':vs_thresh, 'm_thresh':m_thresh, 'num_candidates':num, 'flagged_num_candidates':flagged_num})

def image_detail_reset(request, pk):
    user = request.user
    username = user.get_username()
    image = Image.objects.get(pk=pk)
    processing_options = Processingsettings.objects.get(image_id=pk)
    vs_thresh = settings.CANDIDATE_VS_THRESHOLD
    m_thresh = settings.CANDIDATE_M_THRESHOLD
    num = len(Crossmatches.objects.all().filter(
        image_id=pk
    ).filter(pipelinetag="Good").filter(
        vs__gte=vs_thresh
    ).filter(m_abs__gte=m_thresh))
    flagged_num = len(Crossmatches.objects.all().filter(
        image_id=pk
    ).exclude(pipelinetag="Good").filter(
        vs__gte=vs_thresh
    ).filter(m_abs__gte=m_thresh))
    # username = request.GET['username']
    if user.is_staff or username == image.claimed_by:
        image.claimed_by = "Unclaimed"
        image.save()

    return render(request, 'image_detail.html', {'image':image, "processing_options":processing_options, "claimed_success":False, "claim_attempt":False,
        'vs_thresh':vs_thresh, 'm_thresh':m_thresh, 'num_candidates':num, 'flagged_num_candidates':flagged_num})

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
    transient_sources = Crossmatches.objects.all().filter(usertag=user_tag).order_by("id")
    if user != "all":
        transient_sources = transient_sources.filter(checkedby=user)
    table = CrossmatchTable(transient_sources)
    RequestConfig(request, paginate={'per_page': 100}).configure(table)
    friendly_type = "Crossmatches"

    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('search_results_{}_{}_{}.{}'.format(friendly_type.lower().replace(" ", "-"), user_tag, user, export_format))

    return render(request, "search_results.html", {"transient_type":friendly_type, "user_tag":user_tag, "table":table, "user":user})

def crossmatch_quickview(request,pk):
    image = Image.objects.get(pk=pk)
    try:
        crossmatch_sources = [Crossmatches.objects.get(id=id) for id in request.session['query_ids']]
    except:
        crossmatch_sources = []

    total = len(crossmatch_sources)
    return render(request, 'crossmatch_quickview.html', {"crossmatch_sources":crossmatch_sources, "total":total, "query":False, "image":image})

def crossmatch_quickview_query(request):
    try:
        crossmatch_sources = [Crossmatches.objects.get(id=id) for id in request.session['query_ids']]
    except:
        crossmatch_sources = []

    total = len(crossmatch_sources)
    return render(request, 'crossmatch_quickview.html', {"crossmatch_sources":crossmatch_sources, "total":total, "query":True})

def turn_on_hotkeys(request):
    request.session['hotkeys']='true'

def turn_off_hotkeys(request):
    request.session['hotkeys']='false'

def transient_query(request):
    transients = Crossmatches.objects.all()
    transients.image_crossmatch = False
    transients_filter = TransientFilter(request.GET, queryset=transients)
    table = CrossmatchQueryTable(transients_filter.qs)
    vs_thresh = settings.CANDIDATE_VS_THRESHOLD
    m_thresh = settings.CANDIDATE_M_THRESHOLD

    #save the queryset IDs to the session
    query_ids = [t.id for t in transients_filter.qs]
    request.session['query_ids'] = query_ids

    user = request.user
    user_id = user.id

    #Have to be logged in to submit a category and to own the image
    if user.is_authenticated:
        the_url = request.build_absolute_uri()
        if "?" in the_url:
            query_string = the_url.split("?")
            user_obj = User.objects.get(pk=user_id)
            user_obj.profile.last_query = query_string[-1]
            user_obj.save()

    RequestConfig(request).configure(table)
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('query_result.{}'.format(export_format))

    return render(request, 'transients_list.html', {'filter': transients_filter, 'table':table, 'vs_thresh':vs_thresh, 'm_thresh':m_thresh})

def image_crossmatch_query(request, pk):
    vs_thresh = settings.CANDIDATE_VS_THRESHOLD
    m_thresh = settings.CANDIDATE_M_THRESHOLD
    image = Image.objects.get(pk=pk)
    transients = Crossmatches.objects.all().filter(image_id=pk)
    transients.image_crossmatch = True
    transients_filter = TransientFilter(request.GET, queryset=transients)
    table = ImageCrossmatchTable(transients_filter.qs)

    #save the queryset IDs to the session
    query_ids = [t.id for t in transients_filter.qs]
    request.session['query_ids'] = query_ids

    user = request.user
    user_id = user.id

    #Have to be logged in to submit a category and to own the image
    # if user.is_authenticated:
    #     the_url = request.build_absolute_uri()
    #     if "?" in the_url:
    #         query_string = the_url.split("?")
    #         user_obj = User.objects.get(pk=user_id)
    #         user_obj.profile.last_query = query_string[-1]
    #         user_obj.save()

    RequestConfig(request).configure(table)
    export_format = request.GET.get('_export', None)

    if TableExport.is_valid_format(export_format):
        exporter = TableExport(export_format, table)
        return exporter.response('query_result.{}'.format(export_format))

    return render(request, 'image_crossmatches.html', {'filter': transients_filter, 'table':table, 'image':image, 'vs_thresh':vs_thresh, 'm_thresh':m_thresh})


def crossmatch_detail_query(request,cross_id):
    object_from_query={"crossmatches":Crossmatches,}
    title ={"crossmatches":"Crossmatches",}
    html ={"crossmatches":"crossmatches",}

    #Get the query id list to sort out next and previous
    try:
        query_ids = request.session['query_ids']
        total = len(query_ids)
        this_index = query_ids.index(int(cross_id))
        if this_index == 0:
            prev_id = -1
            next_id = query_ids[1]
        elif this_index == total-1:
            prev_id = query_ids[this_index-1]
            next_id = -1
        else:
            prev_id = query_ids[this_index-1]
            next_id = query_ids[this_index+1]
    except:
        #likely a user has clicked a link pasted somewhere
        next_id=-1
        prev_id=-1
        total = 0
        this_index = -1

    allsources = Crossmatches.objects.all()

    #See if we are turning on hotkeys
    try:
        hotkeys_on = request.GET["hotkeys_on"]
        request.session["hotkeys"]="true"
        # return render(request, 'crossmatch_detail_query.html', {"crossmatch_sources":allsources,"image":image, "querytype":querytype, "title":title_to_use, "total":total, "html":html[querytype],
        #     "subquery_type":subquery_type, "transient_filter":transient_filter})
    except:
        pass

    try:
        hotkeys_off = request.GET["hotkeys_off"]
        request.session["hotkeys"]="false"
        # return render(request, 'crossmatch_quickview.html', {"crossmatch_sources":allsources,"image":image, "querytype":querytype, "title":title_to_use, "total":total, "html":html[querytype],
        #     "subquery_type":subquery_type, "transient_filter":transient_filter})
    except:
        pass

    #Check if a navigation hotkey has been used
    try:
        go_to = request.GET['go_to']
        if go_to == "next":
            go_to_id = next_id
        elif go_to == "previous":
            go_to_id = prev_id
        return redirect('crossmatch_detail_query', cross_id=go_to_id)
    except:
        go_to = False

    crossmatch_source = Crossmatches.objects.get(id=cross_id)
    detail_table = CrossmatchDetailTable(allsources.filter(id=cross_id))
    ratio_table = CrossmatchDetailFluxTable(allsources.filter(id=cross_id))
    image = Image.objects.get(pk=crossmatch_source.image_id)
    #Also as large ratio we can fetch the nearest sources table
    nearest_sources = crossmatch_source.nearest_sources
    nearest_sources = nearest_sources.split(",")
    nearest_sources = Crossmatches.objects.all().filter(image_id=image.id, master_name__in=nearest_sources)
    nearest_sources_table = NearestSourceDetailFluxTable(nearest_sources)
    radius = 5./3600.
    ra=float(crossmatch_source.ra)
    dec=float(crossmatch_source.dec)
    possible_assoc = allsources.extra(where=["q3c_radial_query(ra, dec, {:.6f}, {:.6f}, {:.6f})".format(ra,dec,radius)]).exclude(id=cross_id)
    possible_assoc_table = AssocFluxTable(possible_assoc)
    title_to_use = "Crossmatches"
    url_to_use = "crossmatches"
    simbad_query="http://simbad.u-strasbg.fr/simbad/sim-coo?CooEpoch=2000&Coord={}d{}d&Radius.unit=arcmin&CooEqui=2000&CooFrame=FK5&Radius=10".format(crossmatch_source.ra, crossmatch_source.dec)
    ned_query="https://ned.ipac.caltech.edu/conesearch?search_type=Near%20Position%20Search&coordinates={}d%20{}d&radius=2.0&in_csys=Equatorial&in_equinox=J2000.0&out_csys=Equatorial&out_equinox=J2000.0&hconst=67.8&omegam=0.308&omegav=0.692&wmap=4&corr_z=1".format(crossmatch_source.ra, crossmatch_source.dec)
    follow_up_page = "http://ada.physics.usyd.edu.au:8015/docs/ratio_query/{}_{}_followup.html".format(crossmatch_source.image_id, crossmatch_source.id)
    admin_only = ["transient", "maybe"]
    #Check if an assign button has been used.
    try:
        assign = request.GET['assign']
        usertag = request.GET['usertag']
        userreason = request.GET['userreason']
        hotkey = request.GET['hotkey']

        user = request.user
        username = user.get_username()

        #Have to be logged in to submit a category and to own the image
        if not user.is_authenticated:
            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":False, "updated":False, "conflict":False,
                                                            'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                        'query':True, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)

        elif username!=image.claimed_by and not user.is_staff:
            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":False, "updated":False, "conflict":False,
                                                                'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                                'query':True, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)


        elif not user.is_staff and usertag in admin_only:
            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":False, "updated":False, "conflict":False,
                                                                'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                                'query':True, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)

        else:

            if crossmatch_source.checkedby.lower()=="n/a":
                crossmatch_source.checkedby=username
                crossmatch_source.usertag=usertag
                crossmatch_source.userreason=userreason
                crossmatch_source.save()
                saved=True
                updated=False
                conflict=False

                #update image checked:
                transient_sources = Crossmatches.objects.all().filter(image_id=image.id).filter(pipelinetag="Good").filter(ratio__gte=2.0).exclude(checkedby="N/A")
                total_checked = len(list(transient_sources.values_list('id', flat=True)))
                image.number_candidates_checked = total_checked
                image.save()
                if hotkey=="true":
                    if next_id != -1:
                        return redirect(reverse('crossmatch_detail_query', kwargs={"cross_id":next_id}) + '?hotkeyassign=true&update=false')
                    else:
                        return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                            'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                                        'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                                        'query':True, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)
            elif crossmatch_source.checkedby == username:
                crossmatch_source.checkedby=username
                crossmatch_source.usertag=usertag
                crossmatch_source.userreason=userreason
                crossmatch_source.save()
                updated=True
                saved=False
                conflict = False
                #update image checked:
                transient_sources = Crossmatches.objects.all().filter(image_id=image.id).filter(pipelinetag="Good").filter(ratio__gte=2.0).exclude(checkedby="N/A")
                total_checked = len(list(transient_sources.values_list('id', flat=True)))
                image.number_candidates_checked = total_checked
                image.save()
                if hotkey=="true":
                    if next_id != -1:
                        return redirect(reverse('crossmatch_detail_query', kwargs={"cross_id":next_id}) + '?hotkeyassign=true&update=true')
                    else:
                        return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                            'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                                        'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                                        'query':True, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)
            else:
                if user.is_staff:
                    if usertag not in admin_only:
                        crossmatch_source.checkedby=username
                    crossmatch_source.usertag=usertag
                    crossmatch_source.userreason=userreason
                    crossmatch_source.save()
                    # image = Image.objects.get(pk=pk)
                    #update image checked:
                    transient_sources = Crossmatches.objects.all().filter(image_id=image.id).filter(pipelinetag="Good").filter(ratio__gte=2.0).exclude(checkedby="N/A")
                    total_checked = len(list(transient_sources.values_list('id', flat=True)))
                    image.number_candidates_checked = total_checked
                    image.save()
                    updated=True
                    saved=False
                    conflict = False
                    if hotkey=="true":
                        if next_id!=-1:
                            return redirect(reverse('crossmatch_detail_query', kwargs={"cross_id":next_id}) + '?hotkeyassign=true&update=true')
                        else:
                            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                                'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                                            'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                                            'query':True, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)
                else:
                    saved=False
                    updated=False
                    conflict=True

    except:
        assign = False
        saved=False
        updated=False
        conflict=False
        hotkey=False

    #See if we are coming from a hotkey assign
    try:
        hotkeyassign = request.GET['hotkeyassign']
        update = request.GET["update"]
        if hotkeyassign == "true":
            if update == "true":
                hotkey_updated=True
            else:
                hotkey_updated=False
            prev_crossmatch_source = Crossmatches.objects.get(id=prev_id)
            hotkey_usertag = prev_crossmatch_source.usertag
            hotkey_userreason = prev_crossmatch_source.userreason
            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                            'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                            "hotkey_assign":True, "hotkey_prev_id":prev_id, "hotkey_usertag":hotkey_usertag, "hotkey_userreason":hotkey_userreason, "hotkey_updated":hotkey_updated,
                                                            'query':True, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)
    except:
        pass

    try:
        send_slack = request.GET["slack"]
        ploturl = request.GET["ploturl"]
        ploturl = "/static/"+ploturl
        settings_dir = os.path.dirname(__file__)
        PROJECT_ROOT = os.path.abspath(os.path.dirname(settings_dir))
        # print(PROJECT_ROOT)
        ploturl = PROJECT_ROOT + ploturl
        if send_slack == "true":
            username = request.user.get_username()
            client = slack.WebClient(token=settings.SLACK_TOKEN)
            url = request.build_absolute_uri().split("?")[0]
            urltogo = url + "?slack_sent=true"
            # direct users to source via image url, not query
            if "query" in url:
                url = url.replace("/query/", "/image/{}/crossmatches/".format(crossmatch_source.image_id))
            response = client.files_upload(
                channels=settings.SLACK_CHANNEL_ID,
                initial_comment="User *{}* would like to share the following source:\n\nLink: {}\n\nSee the thread for more information!".format(username, url),
                file=ploturl,
                filename=ploturl.split("/")[-1],
                title=crossmatch_source.master_name
                )
            ts = response['file']['shares']['private'][settings.SLACK_CHANNEL_ID][0]['ts']

            blocks = [
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": "Details of {}:".format(crossmatch_source.master_name),
                    }
                },
                {
                    "type": "section",
                    "fields": [
                    {
                        "type": "mrkdwn",
                        "text": "*RA, Dec:*\n{:.03f}, {:.03f}".format(crossmatch_source.ra, crossmatch_source.dec)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*Gal l, b:*\n{:.03f}, {:.03f}".format(crossmatch_source.gal_l, crossmatch_source.gal_b)
                    },
                    ]
                },
                {
                    "type": "divider"
                },
                {
                    "type": "section",
                    "fields": [
                    {
                        "type": "mrkdwn",
                        "text": "*Type:*\n{} ({})".format(crossmatch_source.transient_type, crossmatch_source.survey.upper())
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*Ratio:*\n{:.02f}".format(crossmatch_source.ratio)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*ASKAP Flux (scaled):*\n{:.02f} mJy".format(float(crossmatch_source.ratio_askap_flux)*1.e3)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*{} Flux:*\n{:.02f} mJy".format(crossmatch_source.survey.upper(), float(crossmatch_source.ratio_catalog_flux)*1.e3)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*Vs:*\n{:.02f}".format(crossmatch_source.vs)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*m:*\n{:.02f}".format(crossmatch_source.m)
                    },
                    ]
                },
                {
                    "type": "divider"
                },
                {
                    "type": "section",
                    "fields": [
                    {
                        "type": "mrkdwn",
                        "text": "*Pipeline Tag:*\n{}".format(crossmatch_source.pipelinetag)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*User Tag:*\n{}".format(crossmatch_source.usertag)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*Checked by:*\n{}".format(crossmatch_source.checkedby)
                    },
                    ]
                },
                {
                    "type": "actions",
                    "elements": [
                        {
                            "type": "button",
                            "text": {
                                "type": "plain_text",
                                "text": "SIMBAD",
                            },
                            "url": simbad_query,
                            "style": "primary",
                        },
                        {
                            "type": "button",
                            "text": {
                            "type": "plain_text",
                            "text": "NED",
                            },
                            "url": ned_query,
                            "style": "primary",
                        },
                    ]
                }
            ]
            detail_response = client.chat_postMessage(
                channel=settings.SLACK_CHANNEL_ID,
                blocks=blocks,
                thread_ts=ts,
                # reply_broadcast=True
            )

        return redirect(urltogo)

    except:
        pass

    try:
        slack_sent = request.GET["slack_sent"]
        if slack_sent:
            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                        'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                    'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                    'query':True, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page, 'slack_sent':True},)
    except:
        pass

    return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                        'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                    'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                    'query':True, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)


def image_crossmatches_detail(request,pk,cross_id):
    object_from_query={"crossmatches":Crossmatches,}
    title ={"crossmatches":"Crossmatches",}
    html ={"crossmatches":"crossmatches",}

    #Get the query id list to sort out next and previous
    try:
        query_ids = request.session['query_ids']
        total = len(query_ids)
        this_index = query_ids.index(int(cross_id))
        if this_index == 0:
            prev_id = -1
            next_id = query_ids[1]
        elif this_index == total-1:
            prev_id = query_ids[this_index-1]
            next_id = -1
        else:
            prev_id = query_ids[this_index-1]
            next_id = query_ids[this_index+1]
    except:
        #likely a user has clicked a link pasted somewhere
        next_id=-1
        prev_id=-1
        total = 0
        this_index = -1

    allsources = Crossmatches.objects.all()

    #See if we are turning on hotkeys
    try:
        hotkeys_on = request.GET["hotkeys_on"]
        request.session["hotkeys"]="true"
        # return render(request, 'crossmatch_detail_query.html', {"crossmatch_sources":allsources,"image":image, "querytype":querytype, "title":title_to_use, "total":total, "html":html[querytype],
        #     "subquery_type":subquery_type, "transient_filter":transient_filter})
    except:
        pass

    try:
        hotkeys_off = request.GET["hotkeys_off"]
        request.session["hotkeys"]="false"
        # return render(request, 'crossmatch_quickview.html', {"crossmatch_sources":allsources,"image":image, "querytype":querytype, "title":title_to_use, "total":total, "html":html[querytype],
        #     "subquery_type":subquery_type, "transient_filter":transient_filter})
    except:
        pass

    #Check if a navigation hotkey has been used
    try:
        go_to = request.GET['go_to']
        if go_to == "next":
            go_to_id = next_id
        elif go_to == "previous":
            go_to_id = prev_id
        print(go_to_id)
        return redirect('image_crossmatches_detail', pk=pk, cross_id=go_to_id)
    except:
        go_to = False

    crossmatch_source = Crossmatches.objects.get(id=cross_id)
    detail_table = CrossmatchDetailTable(allsources.filter(id=cross_id))
    ratio_table = CrossmatchDetailFluxTable(allsources.filter(id=cross_id))
    image = Image.objects.get(pk=pk)
    #Also as large ratio we can fetch the nearest sources table
    nearest_sources = crossmatch_source.nearest_sources
    nearest_sources = nearest_sources.split(",")
    nearest_sources = Crossmatches.objects.all().filter(image_id=image.id, master_name__in=nearest_sources)
    nearest_sources_table = NearestSourceDetailFluxTable(nearest_sources)
    radius = 5./3600.
    ra=float(crossmatch_source.ra)
    dec=float(crossmatch_source.dec)
    possible_assoc = allsources.extra(where=["q3c_radial_query(ra, dec, {:.6f}, {:.6f}, {:.6f})".format(ra,dec,radius)]).exclude(id=cross_id)
    possible_assoc_table = AssocFluxTable(possible_assoc)
    title_to_use = "Crossmatches"
    url_to_use = "crossmatches"
    simbad_query="http://simbad.u-strasbg.fr/simbad/sim-coo?CooEpoch=2000&Coord={}d{}d&Radius.unit=arcmin&CooEqui=2000&CooFrame=FK5&Radius=10".format(crossmatch_source.ra, crossmatch_source.dec)
    ned_query="https://ned.ipac.caltech.edu/conesearch?search_type=Near%20Position%20Search&coordinates={}d%20{}d&radius=2.0&in_csys=Equatorial&in_equinox=J2000.0&out_csys=Equatorial&out_equinox=J2000.0&hconst=67.8&omegam=0.308&omegav=0.692&wmap=4&corr_z=1".format(crossmatch_source.ra, crossmatch_source.dec)
    follow_up_page = "http://ada.physics.usyd.edu.au:8015/docs/ratio_query/{}_{}_followup.html".format(crossmatch_source.image_id, crossmatch_source.id)
    admin_only = ["transient", "maybe"]

    #Check if an assign button has been used.
    try:
        assign = request.GET['assign']
        usertag = request.GET['usertag']
        userreason = request.GET['userreason']
        hotkey = request.GET['hotkey']

        user = request.user
        username = user.get_username()

        #Have to be logged in to submit a category and to own the image
        if not user.is_authenticated:
            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":False, "updated":False, "conflict":False,
                                                            'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                        'query':False, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)

        elif username!=image.claimed_by and not user.is_staff:
            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                            'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":False, "updated":False, "conflict":False,
                                                        'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                    'query':False, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)

        elif not user.is_staff and usertag in admin_only:
            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                            'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":False, "updated":False, "conflict":False,
                                                        'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                    'query':False, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)

        else:

            if crossmatch_source.checkedby.lower()=="n/a":
                crossmatch_source.checkedby=username
                crossmatch_source.usertag=usertag
                crossmatch_source.userreason=userreason
                crossmatch_source.save()
                saved=True
                updated=False
                conflict=False

                #update image checked:
                transient_sources = Crossmatches.objects.all().filter(image_id=image.id).filter(pipelinetag="Good").filter(ratio__gte=2.0).exclude(checkedby="N/A")
                total_checked = len(list(transient_sources.values_list('id', flat=True)))
                image.number_candidates_checked = total_checked
                image.save()
                if hotkey=="true":
                    if next_id != -1:
                        return redirect(reverse('image_crossmatches_detail', kwargs={"pk":image.pk,"cross_id":next_id}) + '?hotkeyassign=true&update=false')
                    else:
                        return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                            'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                                        'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                                        'query':False, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)
            elif crossmatch_source.checkedby == username:
                crossmatch_source.checkedby=username
                crossmatch_source.usertag=usertag
                crossmatch_source.userreason=userreason
                crossmatch_source.save()
                updated=True
                saved=False
                conflict = False
                #update image checked:
                transient_sources = Crossmatches.objects.all().filter(image_id=image.id).filter(pipelinetag="Good").filter(ratio__gte=2.0).exclude(checkedby="N/A")
                total_checked = len(list(transient_sources.values_list('id', flat=True)))
                image.number_candidates_checked = total_checked
                image.save()
                if hotkey=="true":
                    if next_id != -1:
                        return redirect(reverse('image_crossmatches_detail', kwargs={"pk":image.pk, "cross_id":next_id}) + '?hotkeyassign=true&update=true')
                    else:
                        return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                            'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                                        'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                                        'query':False, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)
            else:
                if user.is_staff:
                    if usertag not in admin_only:
                        crossmatch_source.checkedby=username
                    crossmatch_source.usertag=usertag
                    crossmatch_source.userreason=userreason
                    crossmatch_source.save()
                    # image = Image.objects.get(pk=pk)
                    #update image checked:
                    transient_sources = Crossmatches.objects.all().filter(image_id=image.id).filter(pipelinetag="Good").filter(ratio__gte=2.0).exclude(checkedby="N/A")
                    total_checked = len(list(transient_sources.values_list('id', flat=True)))
                    image.number_candidates_checked = total_checked
                    image.save()
                    updated=True
                    saved=False
                    conflict = False
                    if hotkey=="true":
                        if next_id!=-1:
                            return redirect(reverse('image_crossmatches_detail', kwargs={"pk":image.pk, "cross_id":next_id}) + '?hotkeyassign=true&update=true')
                        else:
                            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                                'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                                            'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                                            'query':False, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)
                else:
                    saved=False
                    updated=False
                    conflict=True

    except:
        assign = False
        saved=False
        updated=False
        conflict=False
        hotkey=False

    #See if we are coming from a hotkey assign
    try:
        hotkeyassign = request.GET['hotkeyassign']
        update = request.GET["update"]
        if hotkeyassign == "true":
            if update == "true":
                hotkey_updated=True
            else:
                hotkey_updated=False
            prev_crossmatch_source = Crossmatches.objects.get(id=prev_id)
            hotkey_usertag = prev_crossmatch_source.usertag
            hotkey_userreason = prev_crossmatch_source.userreason
            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                                'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                            'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                            "hotkey_assign":True, "hotkey_prev_id":prev_id, "hotkey_usertag":hotkey_usertag, "hotkey_userreason":hotkey_userreason, "hotkey_updated":hotkey_updated,
                                                            'query':False, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)
    except:
        pass

    try:
        send_slack = request.GET["slack"]
        ploturl = request.GET["ploturl"]
        ploturl = "/static/"+ploturl
        settings_dir = os.path.dirname(__file__)
        PROJECT_ROOT = os.path.abspath(os.path.dirname(settings_dir))
        # print(PROJECT_ROOT)
        ploturl = PROJECT_ROOT + ploturl
        if send_slack == "true":
            username = request.user.get_username()
            client = slack.WebClient(token=settings.SLACK_TOKEN)
            url = request.build_absolute_uri().split("?")[0]
            urltogo = url + "?slack_sent=true"
            # direct users to source via image url, not query
            if "query" in url:
                url = url.replace("/query/", "/image/{}/crossmatches/".format(crossmatch_source.image_id))
            response = client.files_upload(
                channels=settings.SLACK_CHANNEL_ID,
                initial_comment="User *{}* would like to share the following source:\n\nLink: {}\n\nSee the thread for more information!".format(username, url),
                file=ploturl,
                filename=ploturl.split("/")[-1],
                title=crossmatch_source.master_name
                )
            ts = response['file']['shares']['private'][settings.SLACK_CHANNEL_ID][0]['ts']

            blocks = [
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": "Details of {}:".format(crossmatch_source.master_name),
                    }
                },
                {
                    "type": "section",
                    "fields": [
                    {
                        "type": "mrkdwn",
                        "text": "*RA, Dec:*\n{:.03f}, {:.03f}".format(crossmatch_source.ra, crossmatch_source.dec)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*Gal l, b:*\n{:.03f}, {:.03f}".format(crossmatch_source.gal_l, crossmatch_source.gal_b)
                    },
                    ]
                },
                {
                    "type": "divider"
                },
                {
                    "type": "section",
                    "fields": [
                    {
                        "type": "mrkdwn",
                        "text": "*Type:*\n{} ({})".format(crossmatch_source.transient_type, crossmatch_source.survey.upper())
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*Ratio:*\n{:.02f}".format(crossmatch_source.ratio)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*ASKAP Flux (scaled):*\n{:.02f} mJy".format(float(crossmatch_source.ratio_askap_flux)*1.e3)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*{} Flux:*\n{:.02f} mJy".format(crossmatch_source.survey.upper(), float(crossmatch_source.ratio_catalog_flux)*1.e3)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*Vs:*\n{:.02f}".format(crossmatch_source.vs)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*m:*\n{:.02f}".format(crossmatch_source.m)
                    },
                    ]
                },
                {
                    "type": "divider"
                },
                {
                    "type": "section",
                    "fields": [
                    {
                        "type": "mrkdwn",
                        "text": "*Pipeline Tag:*\n{}".format(crossmatch_source.pipelinetag)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*User Tag:*\n{}".format(crossmatch_source.usertag)
                    },
                    {
                        "type": "mrkdwn",
                        "text": "*Checked by:*\n{}".format(crossmatch_source.checkedby)
                    },
                    ]
                },
                {
                    "type": "actions",
                    "elements": [
                        {
                            "type": "button",
                            "text": {
                                "type": "plain_text",
                                "text": "SIMBAD",
                            },
                            "url": simbad_query,
                            "style": "primary",
                        },
                        {
                            "type": "button",
                            "text": {
                            "type": "plain_text",
                            "text": "NED",
                            },
                            "url": ned_query,
                            "style": "primary",
                        },
                    ]
                }
            ]
            detail_response = client.chat_postMessage(
                channel=settings.SLACK_CHANNEL_ID,
                blocks=blocks,
                thread_ts=ts,
                # reply_broadcast=True
            )

        return redirect(urltogo)

    except:
        pass

    try:
        slack_sent = request.GET["slack_sent"]
        if slack_sent:
            return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                        'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                    'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                    'query':False, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page, 'slack_sent':True},)
    except:
        pass

    return render(request, 'crossmatch_detail.html', {'crossmatch_source':crossmatch_source, 'image':image, 'type':"crossmatches",
                                                        'title':title_to_use, 'type_url':url_to_use, 'next_id':next_id, 'prev_id':prev_id, 'this_index':this_index+1, 'total':total, "saved":saved, "updated":updated, "conflict":conflict,
                                                    'simbad_query':simbad_query, 'ned_query':ned_query, 'detail_table':detail_table, "ratio_table":ratio_table, 'nearest_sources_table':nearest_sources_table,
                                                    'query':False, 'possible_assoc_table':possible_assoc_table, 'follow_up_page':follow_up_page},)