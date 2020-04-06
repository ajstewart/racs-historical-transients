from django import template
from django.utils.safestring import mark_safe
import math
from astropy import units as u
from astropy.coordinates import SkyCoord

register = template.Library()


SECONDS_IN_MINUTE = 60
SECONDS_IN_HOUR = 60**2
SECONDS_IN_DAY = 24 * SECONDS_IN_HOUR
RADIANS_IN_CIRCLE = 2 * math.pi
ASEC_IN_AMIN = 60
ASEC_IN_DEGREE = 60**2

@register.filter
def deg_to_gal(radec):
    """
    convert deg to galactic coordinates
    """
    c = SkyCoord(radec.replace("_", " "))

    return c.galactic.to_string()

@register.filter
def deg_to_hms(deg):
    """
    convert degrees to sexagesimal hours, minutes, seconds
    """
    rad = float(deg) * (math.pi/180)
    rad %= RADIANS_IN_CIRCLE
    seconds = SECONDS_IN_DAY * rad / RADIANS_IN_CIRCLE
    hours, seconds = divmod(seconds, SECONDS_IN_HOUR)
    minutes, seconds = divmod(seconds, SECONDS_IN_MINUTE)
    return "{:02.0f}:{:02.0f}:{:04.1f}".format(hours, minutes, seconds)

@register.filter
def deg_to_dms(deg):
    """
    convert degrees to sexagesimal  degrees, arecmin and arcsec
    """
    rad = float(deg) * (math.pi/180)
    sign = "+" if rad >= 0 else "-"
    rad = abs(rad)
    seconds = math.degrees(rad) * ASEC_IN_DEGREE
    degrees, seconds = divmod(seconds, ASEC_IN_DEGREE)
    minutes, seconds = divmod(seconds, ASEC_IN_AMIN)
    return "{}{:02.0f}:{:02.0f}:{:02.0f}".format(sign, degrees, minutes, seconds)

@register.filter
def jy_to_mjy(flux):
    """
    convert jansky to millijansky
    """
    return float(flux) * 1.e3

@register.filter(name='subtract')
def subtract(value, arg):
    return value - arg
