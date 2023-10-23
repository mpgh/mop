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
                         'Ksmag', 'e_Ksmag',
                         'W1mag', 'e_W1mag',
                         '+_r', 'pmRA', 'pmDE', 'plx']

    result = vizier_utils.query_vizier(target, 'I/353', columns_list, radius=radius, row_limit=row_limit)

    return result

def estimate_Ksmag(gsc_table):
    """
    Function estimate a magnitude in Ks-band from data in other filters, if no catalog value is available.
    Parameters:
        gsc_entry   vizier table   Table returned from GSC Vizier query
    Returns
        gsc_entry   vizier table    with missing Ksmag, e_Ksmag entries filled
    """


def Ksmag_from_JH(gsc_table):

    idx = np.where(gsc_table['Ksmag'])  # What does Vizier return if there are no data?  Is this a masked array?
    Ksmag = gsc_table['Hmag'].value + 0.2 * (gscFTp['Jmag'].value - gscFTp['Hmag'].value)