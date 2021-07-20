from astroplan import AltitudeConstraint, AirmassConstraint, AtNightConstraint, is_observable
from astroplan import moon_illumination as moon_illum
from astroplan import moon_phase_angle as moon_phase_angle
from astropy.coordinates import AltAz, SkyCoord
from astropy.coordinates import get_moon as get_moon
from astropy.time import Time
import astropy.units as u
import numpy as np
import warnings

from toolbox.LCO_obs_locs import choose_loc


def timeobj(date):
    obs_night = Time(date, format='iso', scale='utc')
    return obs_night


def calculate_visibility(ra, dec, obs_night, observatory, max_airmass=2.0):
    """
    Visibility calculator for a single target.  Constraints include airmass, altitude,
    and AtNight, among others.
    """
    try:
        coords = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')

        obs_loc = choose_loc(observatory)
        if obs_loc is None:
            raise Exception('Please input a valid LCO observatory in string format.')
        else:
            pass

        obs_begin = obs_loc.twilight_evening_astronomical(obs_night, which='nearest')
        obs_end = obs_loc.twilight_morning_astronomical(obs_night, which='next')
        observing_range = [obs_begin, obs_end]
        constraints = [AirmassConstraint(max_airmass), AltitudeConstraint(20*u.deg, 85*u.deg),
                       AtNightConstraint.twilight_astronomical()]
        ever_observable = is_observable(constraints, obs_loc, coords, time_range=observing_range)
        if ever_observable:
            return True
        else:
            return False
    except ValueError:
        print('Your dates were not inputted correctly.')

def all_night_moon_sep(ra, dec, obs_night, observatory, sample_size=25):
    """
    Determines the min and max separations of the target object and the moon over a full
    observing night at the desired observatory. If it registers <15 degree separation at
    minimum, it prints a warning that the target object is too close to the moon.
    If it registers a <15 degree separation at maximum, the observation request is rejected.
    """
    try:
        coords = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')

        obs_loc = choose_loc(observatory)
        # I can either leave this in here, or I can require the input to be choose_loc('Whatever')
        if obs_loc is None:
            raise Exception('Please input a valid LCO observatory in string format.')
        else:
            pass

        obs_begin = obs_loc.twilight_evening_astronomical(obs_night, which='nearest')
        obs_end = obs_loc.twilight_morning_astronomical(obs_night, which='next')
        midnight = obs_loc.midnight(obs_night, which='nearest')
        lower_lim = (obs_begin - midnight).to(u.h)
        upper_lim = (obs_end - midnight).to(u.h)

        delta_midnight = np.linspace(lower_lim.value, upper_lim.value, sample_size)*u.hour
        frame_observing_night = AltAz(obstime=midnight+delta_midnight, location=obs_loc.location)
        targetaltaz_obsnight = coords.transform_to(frame_observing_night)
        moonaltaz_obsnight = get_moon(
            time=midnight+delta_midnight,
            location=obs_loc.location).transform_to(frame_observing_night)

        moon_frac = moon_illum(time=midnight+delta_midnight) * 100
        avg_moonill = np.mean(moon_frac)
        mphase = moon_phase_angle(time=midnight+delta_midnight).to(u.deg)
        avg_mphase = np.mean(mphase)
        # does not go to negatives: simply moves between 0 and 180 degrees, 0 being full moon and 180 being new moon

        sep_array = [y.separation(x) for x, y in zip(targetaltaz_obsnight, moonaltaz_obsnight)]
        sep_array_deg = [x.degree for x in sep_array]
        avg_sep = np.mean(sep_array_deg)

        return sep_array_deg, avg_sep, avg_moonill, avg_mphase

        #if max(sep_array_deg) < 15:
            #raise Exception('Object is too close to the moon on this date.')
        #elif min(sep_array_deg) >= 15:
            #print('Average moon separation is {0:.1f} degrees'.format(avg_sep))
            #print('{0:.1f} percent of the moon is illuminated'.format(avg_moonill))
            #print('The average moon phase angle is {0:.1f}'.format(avg_mphase))
        #else:
            #warnings.warn('Warning: Object is very close to the moon on this date.')
            #print('Average separation is {0:.1f} degrees'.format(avg_sep))
    except (ValueError, AttributeError):
        print('Your dates were not inputted correctly.')
