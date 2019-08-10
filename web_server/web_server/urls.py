"""web_server URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url, include
from django.contrib import admin
from django.contrib.auth import views as auth_views

from images import views
from accounts import views as accounts_views
from django.conf import settings
from django.conf.urls.static import static

urlpatterns = [
    url(r'^$', views.home, name='home'),
    url(r'^signup/$', accounts_views.signup, name='signup'),
    url(r'^login/$', auth_views.LoginView.as_view(template_name='login.html'), name='login'),
    url(r'^logout/$', auth_views.LogoutView.as_view(), name='logout'),
    url(r'^search/$', views.query_queries, name='search'),
    url(r'^comments/', include('django_comments.urls')),
    url(r'^query/$', views.transient_query, name='query'),
    url(r'^query/quickview/$', views.crossmatch_quickview_query, name='query_quickview'),
    url(r'^query/view_source/(?P<cross_id>\d+)/$', views.crossmatch_detail_query, name='crossmatch_detail_query'),
    url(r'^search/results/transient_type=(?P<transient_type>[\w\-]+)&user_tag=(?P<user_tag>[\w\-]+)&user=(?P<user>[\w\-]+)$', views.search_results, name='search_results'),
    url(r'^image/(?P<pk>\d+)/$', views.image_detail, name='image_detail'),
    url(r'^image/(?P<pk>\d+)/claim/$', views.image_detail_claim, name='image_detail_claim'),
    url(r'^image/(?P<pk>\d+)/claim_reset/$', views.image_detail_reset, name='image_detail_reset'),
    url(r'^image/(?P<pk>\d+)/transients/(?P<transient_filter>[\w\-]+)/$', views.transients, name='transients'),
    url(r'^image/(?P<pk>\d+)/(?P<querytype>[\w\-]+)/(?P<transient_filter>[\w\-]+)/quickview/$', views.crossmatch_quickview, name='crossmatch_quickview'),
    url(r'^image/(?P<pk>\d+)/(?P<querytype>[\w\-]+)/view_source/(?P<cross_id>\d+)/$', views.crossmatch_detail, name='crossmatch_detail'),
    url(r'^settings/password/$', auth_views.PasswordChangeView.as_view(template_name='password_change.html'),
        name='password_change'),
    url(r'^settings/password/done/$', auth_views.PasswordChangeDoneView.as_view(template_name='password_change_done.html'),
        name='password_change_done'),
    url(r'^admin/', admin.site.urls)
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
