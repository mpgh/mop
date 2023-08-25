import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.gaia import Gaia
import numpy as np
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from astropy.table import Table, Column

def find_companion_stars(target, star_catalog):
    """Function to identify stars nearby to the target that have JHK brightneses high enough
    for inteferometry.  """

    # Identifies stars with valid photometry in all passbands
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

    # Returns lists of the star subset's photometry in order of ascending distance from
    # the target
    column_list = [
        Column(name='Source', data=star_catalog[0][mask]['Source'][index]),
        Column(name='Gmag', data=star_catalog[0][mask]['Gmag'][index]),
        Column(name='BP-RP', data=star_catalog[0][mask]['BP-RP'][index]),
        Column(name='Separation', data=separations[index])
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
