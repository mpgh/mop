from django.test import TestCase
from tom_targets.tests.factories import SiderealTargetFactory
from tom_dataproducts.models import ReducedDatum
from .test_fittools import generate_test_ReducedDatums
from mop.toolbox import TAP
from astropy.time import Time, TimeDelta
from astropy import units as u
import numpy as np

class TestObservingMode(TestCase):
    def setUp(self):
        self.params = [{
                            'planet_priority': 7.5e-5,
                            'planet_priority_error': 0.026,
                            'long_priority': 0.01,
                            'long_priority_error': 1341260.75767964,
                            'mag_now': 14.292,
                            'mag_baseline': 15.0,
                            'obs_mode': None
                        },
                        {
                            'planet_priority': 15.0,
                            'planet_priority_error': 1.0,
                            'long_priority': 0.01,
                            'long_priority_error': 1e5,
                            'mag_now': 15.5,
                            'mag_baseline': 18.0,
                            'obs_mode': 'priority_stellar_event'
                        },
                        {
                            'planet_priority': 0.01,
                            'planet_priority_error': 10.0,
                            'long_priority': 75.0,
                            'long_priority_error': 5.0,
                            'mag_now': 17.0,
                            'mag_baseline': 18.5,
                            'obs_mode': 'priority_long_event'
                        },
                        {
                            'planet_priority': 0.01,
                            'planet_priority_error': 10.0,
                            'long_priority': 15.0,
                            'long_priority_error': 5.0,
                            'mag_now': 17.0,
                            'mag_baseline': 18.5,
                            'obs_mode': 'regular_long_event'
                        },
                        ]

    def test_TAP_observing_mode(self):
        for event in self.params:
            obs_mode = TAP.TAP_observing_mode(
                event['planet_priority'],
                event['planet_priority_error'],
                event['long_priority'],
                event['long_priority_error'],
                event['mag_now'],
                event['mag_baseline']
            )
            self.assertEqual(obs_mode, event['obs_mode'])

class TestEventLocation(TestCase):
    def setUp(self):
        self.kmtnet_fields = TAP.load_KMTNet_fields()
        self.params = [
            {'test': (268.0, -28.60972),
             'expect': True},
            {'test': (150.0, -28.60972),
             'expect': False},
            {'test': (268.0, 15.0),
             'expect': False},
        ]

    def test_event_not_in_OMEGA_II(self):
        for config in self.params:
            status = TAP.event_not_in_OMEGA_II(config['test'][0], config['test'][1], self.kmtnet_fields)
            assert (status == config['expect'])

class TestLightcurveData(TestCase):
    def setUp(self):
        st1 = SiderealTargetFactory.create()
        tel_configs = {
            'I': [2000, 18.0, 0.01,
                  Time('2023-08-01T00:00:00.0', format='isot'),
                  TimeDelta(0.25*u.day)],
            'G': [100, 17.5, 0.005,
                  Time('2023-08-01T00:00:00.0', format='isot'),
                  TimeDelta(1.0*u.day)]
        }

        self.params = {
            'target': st1,
            'tel_configs': tel_configs
        }

        # Use this configuration to generate some photometry datapoints for testing
        photometry = generate_test_ReducedDatums(self.params['target'], self.params['tel_configs'])

    def test_TAP_time_last_datapoint(self):
        # Retrieve the photometry for the test target.  Note that we retrieve it here rather than
        # pass the photometry array into the test because the DB applies a TimeZone correction
        # which adjusts the datapoints.  So the timestamps of the ReducedDatums are slightly different.
        photometry = ReducedDatum.objects.filter(target=self.params['target'])

        # Review timestamps of generated data
        ts = [Time(rd.timestamp, format='datetime').jd for rd in photometry if rd.data_type == 'photometry']
        ts = np.array(ts)
        idx = np.argsort(ts)
        expected_t_last = ts.max()

        # Calculate the timestamp of the latest datapoint
        t_last = TAP.TAP_time_last_datapoint(self.params['target'])

        # Expect a floating point number in Julian Date
        assert(t_last > 2450000.0)

        # Test the correct most-recent timestamp is returned
        assert(t_last == expected_t_last)