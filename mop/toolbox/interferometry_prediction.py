from tom_dataproducts.models import ReducedDatum
from astropy.coordinates import SkyCoord
from astropy import units as u
from mop.brokers import gaia
from mop.toolbox import utilities
import numpy as np
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from astropy.table import Table, Column
from astropy.time import Time, TimezoneInfo
import logging
from mop import settings

logger = logging.getLogger(__name__)

def unmask_column(mColumn):
    """Gaia catalog searches can return MaskedColumns, which cannot easily be converted. This function
    extracts the data"""

    data = np.zeros(len(mColumn))
    mask = np.where(mColumn.mask == False)
    data[mask] = mColumn.data[mask]

    return data

def find_companion_stars(target, star_catalog):
    """Function to identify stars nearby to the target that have JHK brightneses high enough
    for inteferometry.  """

    # Identifies stars with valid photometry in all passbands
    try:
        mask = ~star_catalog[0]['Gmag'].data.mask & ~star_catalog[0]['BP-RP'].data.mask

        # For this subset of stars, calculates the angular separation between each star and the target
        # Note: this is the wrong formula
        star = SkyCoord(target.ra, target.dec, frame='icrs', unit=(u.deg, u.deg))
        neighbours = SkyCoord(star_catalog[0][mask]['RA_ICRS'].data.data,
                              star_catalog[0][mask]['DE_ICRS'].data.data,
                              frame='icrs', unit=(u.deg, u.deg))
        separations = star.separation(neighbours)

        # Sort into order of ascending distance from the target:
        index = np.argsort(separations)

        # Estimate photometric uncertainties on Gaia photometry
        BPmag_error = unmask_column(star_catalog[0][mask]['e_BPmag'][index])
        RPmag_error = unmask_column(star_catalog[0][mask]['e_RPmag'][index])
        BPRP_error = np.sqrt( BPmag_error*BPmag_error + RPmag_error*RPmag_error )

        # Returns lists of the star subset's photometry in order of ascending distance from
        # the target
        column_list = [
            Column(name='Gaia_Source_ID', data=star_catalog[0][mask]['Source'][index]),
            Column(name='Gmag', data=unmask_column(star_catalog[0][mask]['Gmag'][index])),
            Column(name='Gmag_error', data=unmask_column(star_catalog[0][mask]['e_Gmag'][index])),
            Column(name='BPmag', data=unmask_column(star_catalog[0][mask]['BPmag'][index])),
            Column(name='BPmag_error', data=BPmag_error),
            Column(name='RPmag', data=unmask_column(star_catalog[0][mask]['RPmag'][index])),
            Column(name='RPmag_error', data=RPmag_error),
            Column(name='BP-RP', data=unmask_column(star_catalog[0][mask]['BP-RP'][index])),
            Column(name='BP-RP_error', data=BPRP_error),
            Column(name='Reddening(BP-RP)', data=unmask_column(star_catalog[0][mask]['E_BP-RP_'][index])),
            Column(name='Extinction_G', data=unmask_column(star_catalog[0][mask]['AG'][index])),
            Column(name='Distance', data=unmask_column(star_catalog[0][mask]['Dist'][index])),
            Column(name='Teff', data=unmask_column(star_catalog[0][mask]['Teff'][index])),
            Column(name='logg', data=unmask_column(star_catalog[0][mask]['logg'][index])),
            Column(name='[Fe/H]', data=unmask_column(star_catalog[0][mask]['__Fe_H_'][index])),
            Column(name='RUWE', data=unmask_column(star_catalog[0][mask]['RUWE'][index])),
            Column(name='Separation', data=separations[index])
        ]
    except IndexError:
        column_list = [
            Column(name='Gaia_Source_ID', data=np.array([])),
            Column(name='Gmag', data=np.array([])),
            Column(name='Gmag_error', data=unp.array([])),
            Column(name='BPmag', data=np.array([])),
            Column(name='BPmag_error', data=np.array([])),
            Column(name='RPmag', data=np.array([])),
            Column(name='RPmag_error', data=np.array([])),
            Column(name='BP-RP', data=np.array([])),
            Column(name='BP-RP_error', data=np.array([])),
            Column(name='Reddening(BP-RP)', data=np.array([])),
            Column(name='Extinction_G', data=np.array([])),
            Column(name='Distance', data=np.array([])),
            Column(name='Teff', data=np.array([])),
            Column(name='logg', data=np.array([])),
            Column(name='[Fe/H]', data=np.array([])),
            Column(name='RUWE', data=np.array([])),
            Column(name='Separation', data=np.array([]))
        ]
        
    return Table(column_list)

def convert_Gmag_to_JHK(Gmag,BpRp):
    """Function to convert Gaia photometry in G band, together with colour photometry in BP-RP to
    J, H, Ks filters, using the conversion from
    https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html#Ch5.T8
    """
    J = Gmag - (-0.01798 + 1.389 * BpRp - 0.09338 * BpRp ** 2)
    H = Gmag - (-0.1048 + 2.011 * BpRp - 0.1758 * BpRp ** 2)
    K = Gmag - (-0.0981 + 2.089 * BpRp - 0.1579 * BpRp ** 2)

    # If the input magnitudes are Table arrays, ensure the output
    # tables column names are accurately renamed
    if type(Gmag) == type(Column()):
        J = Column(name='Jmag', data=J.data)
        H = Column(name='Hmag', data=H.data)
        K = Column(name='Kmag', data=K.data)
    return J, H, K

def estimate_target_Gaia_phot_uncertainties(Gmag, u0, u0_error):
    """Function to estimate the photometric uncertainties for the Gaia photometry of a given star"""

    # Simulate a distribution of event baseline magnitudes and colours
    g_sample = np.random.normal(Gmag, 0.1, 10000)

    # Simulate a distribution of event impact parameters, and calculate the corresponding
    # peak magnification
    u0_sample = abs(np.random.normal(u0, u0_error, 10000))
    magnification = peak_magnification(u0_sample)

    # Calculate the corresponding distribution of G-band peak event magnitudes,
    # and estimate the uncertainty from the width of the distribution
    g_sample -= 2.5 * np.log10(magnification)
    Gmag_error = g_sample.std()

    return Gmag_error

def interferometry_decision(G_lens, BPRP_lens, K_neighbours):

    mode = 'No'
    guide = 0

    if len(K_neighbours) > 0:
        g = np.random.normal(G_lens, 0.1, 10000)
        bprp = np.random.normal(BPRP_lens, 0.1, 10000)
        (J_lens, K_lens, H_lens) = convert_Gmag_to_JHK(g, bprp)

        percentiles = np.percentile(K_lens,[16,50,84])
        median = percentiles[1]
        std = (percentiles[2]-percentiles[0])/2
        if (median+2*std<10) & (std<1):

            mode = 'Single Field'

        else:

            if np.min(K_neighbours)<11:

                mode = 'Dual Field Wide'
                guide = np.argmin(K_neighbours)
    return mode,guide        
        
            
def GAIA_toJHK(G,BpRp):
    
    ### https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html#Ch5.T8
    
    try:
    
        J = G-(-0.01798+1.389*BpRp-0.09338*BpRp**2)+np.random.normal(0,0.04762,len(BpRp))
        H = G-(-0.1048+2.011*BpRp-0.1758*BpRp**2)+np.random.normal(0,0.07805,len(BpRp))
        K = G-(-0.0981+2.089*BpRp-0.1579*BpRp**2)+np.random.normal(0,0.08553,len(BpRp))
    
    except:
        
        J = G-(-0.01798+1.389*BpRp-0.09338*BpRp**2)
        H = G-(-0.1048+2.011*BpRp-0.1758*BpRp**2)
        K = G-(-0.0981+2.089*BpRp-0.1579*BpRp**2)
    
    
    mask = ~np.isnan(J) & ~np.isinf(J) &  ~np.isnan(H) & ~np.isinf(H) & ~np.isnan(K) & ~np.isinf(K)            
    
    return J[mask],H[mask],K[mask]
    

def peak_magnification(uuu):

    A=(uuu**2+2)/(uuu*(uuu**2+4)**0.5)
    
    return A
        
def interfero_plot():
    ### Example with Gaia19bld

    #RA = 189.38565
    #DEC = -66.11136
    #U0 = 0.0193
    #EU0 = 0.0001


    ### Example with Gaia22ahy

    #RA = 225.25274
    #DEC = -54.40000
    #U0 = 0.55601
    #EU0 = 0.2451797

    ### Example with Kojima-1

    RA = 76.925
    DEC = 24.79888
    U0 = 0.08858
    EU0 = 0.0003


    # Queries Gaia DR3 catalog for the target object and near neighbours
    Vizier.ROW_LIMIT = -1

    coord = SkyCoord(ra=RA, dec=DEC, unit=(u.degree, u.degree), frame='icrs')
    result = Vizier.query_region(coord, radius=Angle(0.004, "deg"), catalog='I/355/gaiadr3')

    # Identifies stars with valid photometry in all passbands
    mask = ~result[0]['Gmag'].data.mask & ~result[0]['BP-RP'].data.mask

    # For this subset of stars, calculates the angular separation between each star and the target
    # Note: this is the wrong formula
    distance = np.sqrt((result[0][mask]['RA_ICRS'].data.data-RA)**2+(result[0][mask]['DE_ICRS'].data.data-DEC)**2)

    # Returns lists of the star subset's photometry in order of ascending distance from
    # the target
    G = result[0][mask]['Gmag'][distance.argmin()]
    BP_RP = result[0][mask]['BP-RP'][distance.argmin()]

    # Generates a set of simulated events?
    g = np.random.normal(G,0.1,10000)
    bp_rp = np.random.normal(BP_RP,0.1,10000)
    uo = np.random.normal(U0,EU0,10000)
    magnification = peak_magnification(uo)
    g -= 2.5*np.log10(magnification)

    # Converts the set of stars from Gaia photometric bands to J,H,K
    J,H,K = GAIA_toJHK(g,bp_rp)

    plt.hist(J,25,alpha=0.5,label='J')
    plt.hist(H,25,alpha=0.5,label='H')
    plt.hist(K,25,alpha=0.5,label='Ks')
    plt.legend()
    plt.xlabel('Mag')
    plt.show()


    JJ = []
    HH = []
    KK = []

    for close_stars_index in distance.argsort()[1:]:


        jjj,hhh,kkk = GAIA_toJHK(result[0][mask]['Gmag'][close_stars_index],result[0][mask]['BP-RP'][close_stars_index])

        JJ.append(jjj)
        HH.append(hhh)
        KK.append(kkk)

        if kkk<11:

            break


    mode,guide = interferometry_decision(K,np.array(KK))

    print(mode,guide)

def evaluate_target_for_interferometry(target):
    """Function to calculate the necessary parameters in order to determine whether this target
    is suitable for interferometry"""
    logger.info('INTERFERO: Evaluating event ' + target.name)

    # Extract the necessary target parameters and sanity check.  For reasons I don't understand,
    # u0_error is sometimes stored as a tag rather than an extra parameter during testing.
    u0 = utilities.fetch_extra_param(target, 'u0')
    u0_error = utilities.fetch_extra_param(target, 'u0_error')
    if u0 == None or u0_error == None:
        logger.info('INTERFERO: Insufficient u0 info to evaluate target')
        return

    # Query the Gaia DR3 catalog for the target object and near neighbours within 14.4 arcsec
    # Search radius is set by the interferometry requirement to have bright neighbouring stars
    #         # that are accessible to the GRAVITY instrument
    neighbour_radius = Angle(20.0/3600.0, "deg")
    star_catalog = gaia.query_gaia_dr3(target, radius=neighbour_radius)
    neighbours = find_companion_stars(target, star_catalog)
    logger.info('INTERFERO: Identified ' + str(len(neighbours)) + ' stars in the neighbourhood of '+target.name)

    # The neighbours table is order in ascending separation from the target coordinates,
    # so the target should be the first entry.  Extract the target's Gaia photometry and estimate
    # the uncertainties
    G_lens = neighbours['Gmag'][0]
    BPRP_lens = neighbours['BP-RP'][0]
    G_lens_error = estimate_target_Gaia_phot_uncertainties(G_lens, u0, u0_error)
    logger.info('INTERFERO: Calculated uncertainties Gmag=' + str(G_lens)
                + '+/-' + str(G_lens_error) + 'mag')

    # Predict the NIR photometry of all stars in the region
    (J, H, K) = convert_Gmag_to_JHK(neighbours['Gmag'], neighbours['BP-RP'])
    logger.info('INTERFERO: Computed JHK photometry for neighbouring stars')

    # Evaluate whether this target is suitable for inteferometry
    (mode, guide) = interferometry_decision(G_lens, BPRP_lens, np.array(K.data)[1:])
    logger.info('INTERFERO: Evaluation for interferometry for ' + target.name + ': ' + str(mode) + ' guide=' + str(guide))

    # Store the results
    extras = {
            'Gaia_Source_ID': neighbours['Gaia_Source_ID'][0],
            'Gmag': G_lens, 'Gmag_error': G_lens_error,
            'RPmag': neighbours['RPmag'][0], 'RPmag_error': neighbours['RPmag_error'][0],
            'BPmag': neighbours['BPmag'][0], 'BPmag_error': neighbours['BPmag_error'][0],
            'BP-RP': BPRP_lens, 'BP-RP_error': neighbours['BP-RP_error'][0],
            'Reddening(BP-RP)': neighbours['Reddening(BP-RP)'][0],
            'Extinction_G': neighbours['Extinction_G'][0],
            'Distance': neighbours['Distance'][0],
            'Teff': neighbours['Teff'][0],
            'logg': neighbours['logg'][0],
            '[Fe/H]': neighbours['[Fe/H]'][0],
            'RUWE': neighbours['RUWE'][0],
            'Interferometry_mode': mode,
            'Interferometry_guide_star': guide,
            'Mag_peak_J': J[0],
            'Mag_peak_H': H[0],
            'Mag_peak_K': K[0],
              }
    target.save(extras=extras)

    # Repackage data into a convenient form for storage
    datum = {
        'Gaia_Source_ID': [str(x) for x in neighbours['Gaia_Source_ID']],
        'Gmag': [x for x in neighbours['Gmag']],
        'Gmag_error': [x for x in neighbours['Gmag_error']],
        'BPmag' : [x for x in neighbours['BPmag']],
        'BPmag_error' : [x for x in neighbours['BPmag_error']],
        'RPmag' : [x for x in neighbours['RPmag']],
        'RPmag_error' : [x for x in neighbours['RPmag_error']],
        'BP-RP': [x for x in neighbours['BP-RP']],
        'BP-RP_error': [x for x in neighbours['BP-RP_error']],
        'Jmag': [x for x in J],
        'Hmag': [x for x in H],
        'Kmag': [x for x in K],
        'Reddening(BP-RP)': [x for x in neighbours['Reddening(BP-RP)']],
        'Extinction_G': [x for x in neighbours['Extinction_G']],
        'Distance': [x for x in neighbours['Distance']],
        'Teff': [x for x in neighbours['Teff']],
        'logg': [x for x in neighbours['logg']],
        '[Fe/H]': [x for x in neighbours['[Fe/H]']],
        'RUWE': [x for x in neighbours['RUWE']],
        'Separation': [x for x in neighbours['Separation']],
    }

    # To avoid accumulating entries, search for any existing tabular
    # data for this object and remove it from the DB:
    qs = ReducedDatum.objects.filter(target=target)
    for rd in qs:
        if rd.data_type == 'tabular' and rd.source_name == 'Interferometry_predictor':
            rd.delete()

    # Now store the tabular results
    tnow = Time.now()
    rd, created = ReducedDatum.objects.get_or_create(
        timestamp=tnow.to_datetime(timezone=TimezoneInfo()),
        value=datum,
        source_name='Interferometry_predictor',
        source_location=target.name,
        data_type='tabular',
        target=target)
    logger.info('INTERFERO: Stored neighbouring star data in MOP')
