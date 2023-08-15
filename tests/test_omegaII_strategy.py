from django.test import TestCase
from tom_targets.models import Target
from tom_targets.tests.factories import SiderealTargetFactory
from astropy.time import Time, TimeDelta
import os

from mop.toolbox import omegaII_strategy

class TestObsConfig(TestCase):
    def setUp(self):
        time_now = Time.now()
        tE = TimeDelta(30.0, format='days')
        self.params = {
                        'target': SiderealTargetFactory.create(),
                        'event_in_hcz': False,
                        'event_category': 'Stellar',
                        'current_mag': 16.0,
                        'tE': tE,
                        't0': Time((time_now + 0.1*tE), format='isot')
                    }

    def test_get_default_obs_config(self):
        config = omegaII_strategy.get_default_obs_config(self.params['target'])
        self.assertTrue(config['observation_mode']=='NORMAL')
        self.assertTrue(config['instrument_type']=='1M0-SCICAM-SINISTRO')
        self.assertTrue(config['proposal']==os.getenv('LCO_PROPOSAL_ID'))
        self.assertTrue(config['facility']=='LCO')
        self.assertTrue(config['max_airmass']==2.0)
        self.assertTrue(config['target_id']==target.id)

    def test_determine_obs_config(self):
        default_config = omegaII_strategy.get_default_obs_config(self.params['target'])
        configs = omegaII_strategy.determine_obs_config(self.params['target'],
                                                        self.params['event_in_hcz'],
                                                        self.params['event_category'],
                                                        self.params['current_mag'],
                                                        self.params['t0'],
                                                        self.params['tE'])

        # Test that this target results in the expected two configurations:
        self.assertTrue(len(configs) == 2)

        # Test that the contents of the list of configurations is a set
        # of dictionaries with the following keys populated:
        expected_keys = ['filter', 'ipp_value','start','end','name',
                         'period', 'jitter','exposure_time','exposure_count']
        for conf in configs:
            for key in expected_keys:
                self.assertTrue(key in conf.keys())

