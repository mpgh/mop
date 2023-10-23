from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import astropy.units as u
from astropy.coordinates import SkyCoord
from mop.brokers import vizier_utils

def query_gsc(target, radius=Angle(0.0083, "deg"), row_limit=-1):
    """Function to perform a Vizier search of the Guide Star Catalogue, I/353"""

    columns_list = ['GSC2', 'RA_ICRS', 'DE_ICRS',
                         'Gmag', 'e_Gmag',
                         'Bmag', 'e_Bmag',
                         'Vmag', 'e_Vmag',
                         'Jmag', 'e_Jmag',
                         'Hmag', 'e_Hmag',
                         'Ksmag', 'e_Ksmag']

    result = vizier_utils.query_vizier(target, 'I/353', columns_list, radius=radius, row_limit=row_limit)

    return result