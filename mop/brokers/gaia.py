from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target,TargetExtra
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import astropy.units as u
from astropy.coordinates import SkyCoord

#for a given mag computes new error-bar
#from Gaia DR2 papers, degraded by x10 (N=100 ccds), in log
def estimateGaiaError(mag) :

    a1=0.2
    b1= -5.3#-5.2
    log_err1 = a1*mag + b1
    a2=0.2625
    b2= -6.3625#-6.2625
    log_err2 = a2*mag + b2

    if (mag<13.5): expectedStdAtBaselineMag = 10**(a1*13.5+b1)
    if (mag>=13.5 and mag<17) : expectedStdAtBaselineMag = 10**log_err1
    if (mag>=17) : expectedStdAtBaselineMag = 10**log_err2
    #this works until 21 mag.

    return expectedStdAtBaselineMag

def update_gaia_errors(target):

    datasets = ReducedDatum.objects.filter(target=target)




    for i in datasets:

        if (i.data_type == 'photometry') & ("error" not in i.value.keys())  &  ('Gaia' in i.source_name):
           
            magnitude = i.value['magnitude']
            error = estimateGaiaError(magnitude)
            i.value['error'] = error 
            
            i.save()

def query_gaia_dr3(target, radius=Angle(0.004, "deg")):
    """Function to query the Gaia DR3 catalog for information on a target and stars nearby"""

    Vizier.ROW_LIMIT = -1
    coord = SkyCoord(ra=target.ra, dec=target.dec, unit=(u.deg, u.deg), frame='icrs')
    result = Vizier.query_region(coord, radius=radius, catalog='I/355/gaiadr3')

    return result

def fetch_gaia_dr3_entry(target):
    """Function to retrieve the Gaia photometry for a target and store it in the Target's ExtraParameters"""

    # Search the Gaia DR3 catalog:
    results = query_gaia_dr3(target)
    print(results[0])

    if len(results) > 0:

        extra_params = target.extra_fields
        print('INIT: ',extra_params)

        #extra_params = {}
        fields = {
            'Gmag': 'Gmag',
            'e_Gmag': 'Gmag_error',
            'RPmag': 'RPmag',
            'e_RPmag': 'RPmag_error',
            'BPmag': 'BPmag',
            'e_BPmag': 'BPmag_error',
            'BP-RP': 'BP-RP',
            'E(BP-RP)': 'Reddening(BP-RP)',
            'AG': 'Extinction_G',
            'Dist':'Distance',
            'Teff': 'Teff',
            'logg': 'logg',
            '[Fe/H]': '[Fe/H]'}
        for cat_field, mop_field in fields.items():
            try:
                if results[0][0][cat_field]:
                    extra_params[mop_field] = results[0][0][cat_field]
                    print(cat_field, results[0][0][cat_field])
            except KeyError:
                pass
        print('EXTRAS: ',extra_params)
        target.save(extras = extra_params)
