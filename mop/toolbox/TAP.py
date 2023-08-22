from django.core.management.base import BaseCommand
from tom_observations.models import ObservationRecord
from tom_dataproducts.models import ReducedDatum

from tom_targets.models import Target
from astropy import units as u
from astropy.coordinates import Angle
from astropy.time import Time
import datetime
import numpy as np
from pyLIMA import event
from pyLIMA import telescopes

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


ZP = 27.4 #pyLIMA convention

def TAP_anomaly():

    pass

def TAP_observing_mode(planet_priority, planet_priority_error,
                           long_priority, long_priority_error,
                           mag_now, mag_baseline):

    if (planet_priority>10) & (planet_priority/planet_priority_error>3) & (mag_baseline-mag_now>2) & (mag_now<19): #mag cut for high blended events
        return 'priority_stellar_event'

    elif (long_priority > 50) & (mag_now < 19) & (mag_baseline < 19):
        return 'priority_long_event'

    elif (long_priority > 10) & (mag_now < 19) & (mag_baseline < 19):
        return 'regular_long_event'

    else:
        return None

def calculate_exptime_floyds(magin):
    """
    This function calculates the required exposure time
    for a given iband magnitude for the floyds spectra
    """
    exposure_time = 3600 #s

    if magin<11:
       exposure_time = 1800 #s

    return exposure_time


def calculate_exptime_omega_sdss_i(magin):
    """
    This function calculates the required exposure time
    for a given iband magnitude (e.g. OGLE I which also
    roughly matches SDSS i) based on a fit to the empiric
    RMS diagram of DANDIA light curves from 2016. The
    output is in seconds.
    """
    if magin < 14.7:
        mag = 14.7
    else:
        mag = magin
    lrms = 0.14075464 * mag * mag - 4.00137342 * mag + 24.17513298
    snr = 1.0 / np.exp(lrms)


    # target 4% -> snr 25
    exptime = np.round((25. / snr)**2 * 300.,1)

    #Scaling for bright events
    if magin<14.7:

        exptime *= 10**((magin-mag)/2.5)

    #no need to more 5 min exposure time, since we have different apertures, but more than 5 s at least

    exptime = float(np.max((5,np.min((exptime,300)))))
    return  exptime



def event_in_the_Bulge(ra,dec):

    Bulge_limits = [[255,275],[-36,-22]]

    if (ra>Bulge_limits[0][0]) & (ra<Bulge_limits[0][1]) & (dec>Bulge_limits[1][0]) & (dec<Bulge_limits[1][1]):
        in_the_Bulge = True
    else:

        in_the_Bulge = False

    return in_the_Bulge

def psi_derivatives_squared(t,te,u0,t0):
    """if you prefer to have the derivatives for a simple
       error propagation without correlation
    """
    x0 = u0**2
    x1 = te**(-2)
    x2 = (t - t0)**2
    x3 = x1*x2
    x4 = x0 + x3
    x5 = x2/te**3
    x6 = x4*x5
    x7 = x4 + 4.0
    x8 = x5*x7
    x9 = (x4*x7)**0.5
    x10 = 1/(x4*x7)
    x11 = x10/x9
    x12 = x10*x9
    x13 = 0.125/(0.5*x0 + 0.5*x3 + 0.5*x9 + 1)**2
    x14 = u0*x4
    x15 = u0*x7
    x16 = x1*(-2*t + 2*t0)
    x17 = (1/2)*x16
    x18 = x17*x4
    x19 = x17*x7

    c0 = 16.0*(x11*(x6 + x8) - x13*(-x12*(-x6 - x8) + 2*x5))**2
    c1 = 16.0*(x11*(-x14 - x15) - x13*(-2*u0 - x12*(x14 + x15)))**2
    c2 = 16.0*(x11*(-x18 - x19) - x13*(-x12*(x18 + x19) - x16))**2
    #i.e. for te, u0, to
    return [c0, c1, c2 ]

def TAP_planet_priority_error(time_now,t0_pspl,u0_pspl,tE_pspl,covariance):
    """
    This function calculates the priority for ranking
    microlensing events based on the planet probability psi
    as defined by Dominik 2009 and estimates the cost of
    observations based on an empiric RMS estimate
    obtained from a DANDIA reduction of K2C9 Sinistro
    observations from 2016. It expects the overhead to be
    60 seconds and also return the current Paczynski
    light curve magnification.
    """

    usqr = u0_pspl**2 + ((time_now - t0_pspl) / tE_pspl)**2

    dpsipdu = -8*(usqr+2)/(usqr*(usqr+4)**1.5)
    dpsipdu += 4*(usqr**0.5*(usqr+4)**0.5+usqr+2)/(usqr+2+(usqr+4)**0.5)**2*1/(usqr+4)**0.5




    dUdto = -(time_now - t0_pspl) / (tE_pspl ** 2 *usqr**0.5)
    dUduo = u0_pspl/ usqr**0.5
    dUdtE = -(time_now - t0_pspl) ** 2 / (tE_pspl ** 3 * usqr**0.5)

    Jacob = np.zeros(len(covariance))
    Jacob[0] = dpsipdu*dUdto
    Jacob[1] = dpsipdu*dUduo
    Jacob[2] = dpsipdu*dUdtE


    error_psip = np.dot(Jacob.T,np.dot(covariance,Jacob))**0.5

    return error_psip #/ (calculate_exptime_omega_sdss_i(mag) + 60.)

def TAP_planet_priority(time_now,t0_pspl,u0_pspl,tE_pspl):
    """
    This function calculates the priority for ranking
    microlensing events based on the planet probability psi
    as defined by Dominik 2009 and estimates the cost of
    observations based on an empiric RMS estimate
    obtained from a DANDIA reduction of K2C9 Sinistro
    observations from 2016. It expects the overhead to be
    60 seconds and also return the current Paczynski
    light curve magnification.
    """
    usqr = u0_pspl**2 + ((time_now - t0_pspl) / tE_pspl)**2
    pspl_deno = (usqr * (usqr + 4.))**0.5
    if pspl_deno < 1e-10:
        pspl_deno = 10000.
    psip = 4.0 / (pspl_deno) - 2.0 / (usqr + 2.0 + pspl_deno)

    return psip




def TAP_regular_mode(in_the_Bulge,survey_cadence,sdssi_baseline,tE_fit):

    cadence =   20/tE_fit  #visits per day
    cadence = np.max(0.5,np.min(4,cadence)) # 0.5<cadence<4

    #Inside the Bulge?
    if not in_the_Bulge:
        return cadence

    if survey_cadence<1:

       if sdssi_baseline<16.5:

          return cadence

    return None


def TAP_priority_mode():

    cadence = 24 #pts/day
    duration = 3
    return duration,cadence



def TAP_telescope_class(sdss_i_mag):

    telescope_class = '2m'

    #change telescope class limit to 18.5 to save 2 m time
    if sdss_i_mag<18.0:

        telescope_class = '1m'

    # RAS: Switched off for 2023+ KP since no time allocated
    # on 0.4m network
    use_04_network = False
    if use_04_network:
        if sdss_i_mag<14:
            telescope_class = '0.4m'

    return telescope_class

#def TAP_mag_now(target):

#   fs = 10**((ZP-target.extra_fields['Source_magnitude'])/2.5)
#
#   try:
#       fb = 10**((ZP-target.extra_fields['Blend_magnitude'])/2.5)

#   except:
#      fs = 10**((ZP-target.extra_fields['Baseline_magnitude'])/2.5)
#      fb = 0
#   fit_parameters = [target.extra_fields['t0'],target.extra_fields['u0'],target.extra_fields['tE'],
#                     target.extra_fields['piEN'],target.extra_fields['piEE'],
#                     fs,
#                     fb]
#
#   current_event = event.Event()
#   current_event.name = 'MOP_to_fit'

#   current_event.ra = target.ra
#   current_event.dec = target.dec

#   time_now = Time(datetime.datetime.now()).jd
#   fake_lightcurve = np.c_[time_now,14,0.01]
#   telescope = telescopes.Telescope(name='Fake', camera_filter='I',
#                                            light_curve_magnitude= fake_lightcurve,
#                                            clean_the_lightcurve='No')
#   current_event.telescopes.append(telescope)
#   t0_par = fit_parameters[0]

#   Model_parallax = microlmodels.create_model('PSPL', current_event, parallax=['Full', t0_par],blend_flux_ratio=False)
#   Model_parallax.define_model_parameters()
#   pyLIMA_parameters = Model_parallax.compute_pyLIMA_parameters(fit_parameters)
#   ml_model, f_source, f_blending = Model_parallax.compute_the_microlensing_model(telescope, pyLIMA_parameters)
#
#   mag_now = ZP-2.5*np.log10(ml_model)
#   return mag_now

def TAP_mag_now(target):
    lightcurve = ReducedDatum.objects.filter(target=target,data_type='lc_model')
    time_now = Time(datetime.datetime.now()).jd

    if len(lightcurve) > 0:
        closest_mag = np.argmin(np.abs(lightcurve[0].value['lc_model_time']-time_now))

        mag_now =  lightcurve[0].value['lc_model_magnitude'][closest_mag]
    else:
        mag_now = None

    return mag_now
   
def load_KMTNet_fields():
    """
    Returns a numpy array with the vertices
    describing the KMTNet zone polygon.
    """
    fields = np.array([[264.00, -37.40],
                        [270.50, -37.40],
                        [270.50, -33.75],
                        [272.50, -33.75],
                        [272.50, -32.00],
                        [275.50, -32.00],
                        [275.50, -25.30],
                        [275.60, -25.30],
                        [275.60, -21.90],
                        [272.00, -21.90],
                        [272.00, -23.00],
                        [270.40, -23.00],
                        [270.40, -20.50],
                        [264.50, -20.50],
                        [264.50, -22.70],
                        [262.00, -22.70],
                        [262.00, -26.25],
                        [260.50, -26.25],
                        [260.50, -31.40],
                        [262.00, -31.40],
                        [262.00, -36.00],
                        [264.00, -36.00]])
    return fields
    
def event_not_in_OMEGA_II(ra, dec, KMTNet_fields):

    """
    This function checks if the event is within the KMTNet fields.
    If it is, the event has to be rejected and not followed by OMEGA-II.

    :param ra: Right Ascension of the event.
    :param dec: Declination of the event.
    :param KMTNet_fields: A numpy array that contains a series of
                          points describing roughly
    :return: Boolean value if the event is not ok for OMEGA II.
    """

    exclusion_zone = Polygon(zip(KMTNet_fields[:,0], KMTNet_fields[:,1]))

    not_in_OMEGA_II = False

    point = Point(ra, dec)
    if (exclusion_zone.contains(point)):
        not_in_OMEGA_II = True

    return not_in_OMEGA_II

def TAP_long_event_priority(t_now, t_last, t_E):
    """
    This function calculates priority of a long microlensing event.
    If the event has a timescale equal to t_E_base, and was not
    observed for t_obs_gap days, the priority should be close to 10.
    If the event was not observed for more than 10 days,
    the priority is boosted by 10.
    If the fit was bad (indicated by t_E being a NaN), than
    the priority drops to 0.
    Inspired by Horn et al. 2009.

    :param t_now: current JD
    :param t_last: JD of the latest datapoint in any lightcurve
                    assigned to the dataset
    :param t_E: Einstein timescale of the event

    :return: priority of the long microlensing event
    """

    if (np.isnan(t_E)):
        # t_E may hit the fit bounds. This means, the fit was
        # wrong. In this case the event should not be observed
        # and wait for a better estimate of the values.
        psi = 0.
    else:
        t_E_base = 250.
        t_obs_gap = 2.

        # Events with t_E of t_E_base days will have priority 10.
        psi = t_E / t_E_base

        # If the event was not observed for t_obs_gap days the priority
        # should rise.
        delta_t = t_now - t_last
        x = (delta_t / t_obs_gap) - t_obs_gap
        psi *= 10. * 1. / (1. + np.exp(-x))

        # The event priority should be boosted if it was not
        # observed for more than 10 days and it is a long event.
        if (delta_t > 10. and t_E / t_E_base > 0.8 ):
            psi *= 10.

    return psi

def TAP_long_event_priority_error(t_E, covariance):
    """
    This function calculates error of the priority
    of a long microlensing event.

    :param t_E: Einstein timescale of the event
    :param covariance: covariance matrix

    :return: error of the priority of the long microlensing event.
    """

    t_E_base = 250.

    if(np.isnan(t_E)):
        err_psi = 0.
    else:
        # Simple derevative
        err_t_E = covariance[2,2]
        err_psi = err_t_E / t_E_base

    return err_psi

def TAP_time_last_datapoint(target):
     """
     Returns time of the latest datapoint in the lightcurve.
     """
     datasets = ReducedDatum.objects.filter(target=target)
     time = [Time(i.timestamp).jd for i in datasets if i.data_type == 'photometry']
     sorted_time = np.sort(time)
     t_last = sorted_time[-1]

     return t_last

def categorize_event_timescale(target, threshold=75.0):
    """
    This function is designed to classify an event based on the best
    currently available estimate of its Einstein crossing time.

    Events with a tE >= threshold [days] are considered to be
    long-timescale events and black hole candidates.
    """

    category = 'Microlensing stellar/planet'

    if target.extra_fields['tE'] >= threshold:
        category = 'Microlensing long-tE'

    extras = {'Category': category}
    target.save(extras=extras)

    return category