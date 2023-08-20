import os
from django.test import TestCase
from datetime import datetime, timedelta
from tom_observations.facilities import lco
from tom_targets.tests.factories import SiderealTargetFactory

from mop.toolbox import obs_control

class TestObsConfig(TestCase):
    def setUp(self):
        # Input dictionary format based on example in TOM Toolkit tom_observations/tests/facilities/test_lco.py
        self.st = SiderealTargetFactory.create()
        self.params = [ {
                            'observation_mode': 'NORMAL', *
                            'instrument_type': '1M0-SCICAM-SINISTRO',
                            'proposal': os.getenv('LCO_PROPOSAL_ID'),
                            'facility': 'LCO',
                            'max_airmass': 2.0,
                            'min_lunar_distance': 15.0,
                            'target': self.st
                            'filters': ['gp','ip'],
                            'ipp_value': 0.8,
                            'name': 'TEST_reg_phot_ip',
                            'period': 24.0,
                            'jitter': 24.0,
                            'exposure_times': [30.0, 30.0],
                            'exposure_counts': [1,1],
                            'start': datetime.utcnow(),
                            'end': (datetime.utcnow() + timedelta(days=2.0)),
                         }
                        ]
    def test_build_lco_imaging_request(self):
        obs_request = obs_control.build_lco_imaging_request(config)
        assert(type(obs_request) == type(lco.LCOImagingObservationForm()))
        assert(obs_request.is_valid())

class TestVisibility(TestCase):
    def setUp(self):
        self.st = SiderealTargetFactory.create()
        self.st.ra = 241.2806       # 16:05:07.3
        self.st.dec = -56.4367      # -56:26:12.26
        self.params = [{
                        'target': self.st,
                        'timenow': Time('2023-08-20').decimalyear,
                        'result': True
                       },{
                        'target': self.st,
                        'timenow': Time('2023-01-01').decimalyear,
                        'result': False
                       }]

    def test_check_visibility(self):
        for test_case in self.params:
            visible = obs_control.check_visibility(test_case['target'],
                                                   test_case['timenow'])
            assert(visible == test_case['result'])
            