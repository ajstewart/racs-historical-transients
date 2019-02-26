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
from django.conf.urls import url
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
    url(r'^image/(?P<pk>\d+)/$', views.image_detail, name='image_detail'),
    url(r'^image/(?P<pk>\d+)/sumssnomatch$', views.sumssnomatch, name='sumssnomatch'),
    url(r'^image/(?P<pk>\d+)/largeratio$', views.largeratio, name='largeratio'),
    url(r'^image/(?P<pk>\d+)/askapnotseen$', views.askapnotseen, name='askapnotseen'),
    url(r'^image/(?P<pk>\d+)/goodmatch$', views.goodmatch, name='goodmatch'),
    url(r'^image/(?P<pk>\d+)/(?P<querytype>[\w\-]+)/(?P<cross_id>\d+)$', views.crossmatch_detail, name='crossmatch_detail'),
    url(r'^image/(?P<pk>\d+)/(?P<querytype>[\w\-]+)/(?P<cross_id>\d+)/confirm$', views.crossmatch_commit, name='crossmatch_commit'),
    url(r'^settings/password/$', auth_views.PasswordChangeView.as_view(template_name='password_change.html'),
        name='password_change'),
    url(r'^settings/password/done/$', auth_views.PasswordChangeDoneView.as_view(template_name='password_change_done.html'),
        name='password_change_done'),
    url(r'^admin/', admin.site.urls)
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
