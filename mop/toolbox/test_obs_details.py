from astroplan import AltitudeConstraint, AirmassConstraint, AtNightConstraint, is_observable
from astroplan import moon_illumination as moon_illum
from astroplan import moon_phase_angle as moon_phase_angle
from astropy.coordinates import SkyCoord, AltAz
from astropy.coordinates import get_moon as get_moon
from astropy.time import Time
import astropy.units as u
import numpy as np
from unittest.mock import patch

from django.test import TestCase

from mop.toolbox.obs_details import all_night_moon_sep, calculate_visibility
from mop.toolbox.LCO_obs_locs import choose_loc


OGG = choose_loc('OGG')
v_test_target = ['Sirius', 100.7362500*u.deg, -16.6459444*u.deg]
v_date = Time("2019-12-25 00:00:00", scale='utc')
v_coords = SkyCoord(v_test_target[1], v_test_target[2], frame='icrs')
v_obs_begin = OGG.twilight_evening_astronomical(v_date, which='nearest')
v_obs_end = OGG.twilight_morning_astronomical(v_date, which='next')
v_observing_range = [v_obs_begin, v_obs_end]
constraints = [AirmassConstraint(2.0), AltitudeConstraint(20*u.deg, 85*u.deg),
               AtNightConstraint.twilight_astronomical()]
ever_observable = is_observable(constraints, OGG, v_coords, time_range=v_observing_range)

v_fail_start = Time("2019-12-24 10:00:00", scale='utc')
v_fail_end = Time("2019-12-25 10:00:00", scale='utc')


class TestVisibilityCalc(TestCase):

    def test_timeobj(self):
        self.assertEqual(v_date.scale, 'utc')
        self.assertEqual(v_date.value, '2019-12-25 00:00:00.000')

    def test_coords(self):
        self.assertEqual(v_coords.ra, v_test_target[1])
        self.assertEqual(v_coords.dec, v_test_target[2])

    def test_dates(self):
        tv_obs_begin = v_obs_begin.to_value(format='iso')
        tv_obs_end = v_obs_end.to_value(format='iso')
        self.assertEqual(tv_obs_begin, '2019-12-25 05:10:04.696')
        self.assertEqual(tv_obs_end, '2019-12-25 15:39:42.359')

    def test_ever_obs(self):
        self.assertTrue(ever_observable)

    def test_raise_exceptions(self):
        with self.subTest('Test that an invalid object returns an exception.'):
            vis = calculate_visibility(37.954, 89.264, v_fail_start, v_fail_end, 'OGG')
            self.assertFalse(vis)

        with self.subTest('Test that an exception is raised if the location is incorrect.'):
            with patch('mop.toolbox.obs_details.calculate_visibility') as mock_calculate_visibility:
                mock_calculate_visibility.side_effect = Exception()
                with self.assertRaisesRegex(Exception, 'Please input a valid LCO observatory in string format.'):
                    calculate_visibility(37.954, 89.264, v_fail_start, v_fail_end, OGG)


m_test_target = ['M31', 10.6847083*u.deg, 41.2687500*u.deg]
m_date = Time("2021-05-04 00:00:00", scale='utc')
m_coords = SkyCoord(m_test_target[1], m_test_target[2], frame='icrs')
m_obs_begin = OGG.twilight_evening_astronomical(m_date, which='nearest')
m_obs_end = OGG.twilight_morning_astronomical(m_date, which='next')
midnight = OGG.midnight(m_date, which='nearest')
lower_lim = (m_obs_begin - midnight).to(u.h)
upper_lim = (m_obs_end - midnight).to(u.h)

delta_midnight = np.linspace(lower_lim.value, upper_lim.value, 25)*u.hour
frame_observing_night = AltAz(obstime=midnight+delta_midnight, location=OGG.location)
targetaltaz_obsnight = m_coords.transform_to(frame_observing_night)
moonaltaz_obsnight = get_moon(time=midnight+delta_midnight, location=OGG.location).transform_to(frame_observing_night)

moon_frac = moon_illum(time=midnight+delta_midnight) * 100
avg_moonill = np.mean(moon_frac)
mphase = moon_phase_angle(time=midnight+delta_midnight).to(u.deg)
avg_mphase = np.mean(mphase)

sep_array = [y.separation(x) for x, y in zip(targetaltaz_obsnight, moonaltaz_obsnight)]
sep_array_deg = [x.degree for x in sep_array]
avg_sep = np.mean(sep_array_deg)

m_fail_start = Time("2021-05-03 10:00:00", scale='utc')
m_fail_end = Time("2021-05-04 10:00:00", scale='utc')


class MoonSepCalc(TestCase):

    def test_coords(self):
        self.assertEqual(m_coords.ra, m_test_target[1])
        self.assertEqual(m_coords.dec, m_test_target[2])

    def test_dates(self):
        tm_obs_begin = m_obs_begin.to_value(format='iso')
        tm_obs_end = m_obs_end.to_value(format='iso')
        tm_midnight = midnight.to_value(format='iso')
        self.assertEqual(tm_obs_begin, '2021-05-04 06:09:24.258')
        self.assertEqual(tm_obs_end, '2021-05-04 14:33:56.825')
        self.assertEqual(tm_midnight, '2021-05-04 10:21:40.084')

    def test_good_sep(self):
        f_avg_sep = round(avg_sep, 1)
        f_avg_moonill = round(avg_moonill, 1)
        f_avg_mphase = round(avg_mphase.value, 1)

        with self.subTest('Test the average separation.'):
            self.assertEqual(f_avg_sep, 73.4)

        with self.subTest('Test the average moon illumination.'):
            self.assertEqual(f_avg_moonill, 43.8)

        with self.subTest('Test the average moon phase angle.'):
            self.assertEqual(f_avg_mphase, 97.2)

    def test_raise_exceptions(self):
        with self.subTest('Test that an exception is raised if the location is incorrect.'):
            with patch('mop.toolbox.obs_details.all_night_moon_sep') as mock_all_night_moon_sep:
                mock_all_night_moon_sep.side_effect = Exception()
                with self.assertRaisesRegex(Exception, 'Please input a valid LCO observatory in string format.'):
                    all_night_moon_sep(323.2651667, -19.8063111, m_fail_start, m_fail_end, OGG)
