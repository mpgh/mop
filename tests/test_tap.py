from django.test import TestCase

from mop.toolbox import TAP

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