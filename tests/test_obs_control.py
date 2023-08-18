import os
from django.test import TestCase
from datetime import datetime, timedelta
from tom_observations.facilities import lco
from tom_targets.tests.factories import SiderealTargetFactory

from mop.toolbox import obs_control

class TestObsConfig(TestCase):
    def setUp(self):
        # Missing params group_id, submitter, target(name,ra,dec) -> target.id,
        # Changed names: observation_type -> observation_mode, name -> group_id
        self.params = [
                        {
                            'observation_mode': 'NORMAL',
                            'operator': 'SINGLE',
                            'telescope_class': '1m0',
                            'instrument_type': '1M0-SCICAM-SINISTRO',
                            'proposal': os.getenv('LCO_PROPOSAL_ID'),
                            'facility': 'LCO',
                            'max_airmass': 2.0,
                            'min_lunar_distance': 15.0,
                            'max_lunar_phase': 1.0,
                            'optimization_type': 'TIME',
                            'acceptability_threshold': 100,
                            'configuration_repeats': 1,
                            'target': SiderealTargetFactory.create(),
                            'filters': ['gp','ip'],
                            'ipp_value': 0.8,
                            'name': 'TEST_reg_phot_ip',
                            'period': 24.0,
                            'jitter': 24.0,
                            'exposure_times': [30.0, 30.0],
                            'exposure_counts': [1,1],
                            'start': datetime.utcnow().isoformat(),
                            'end': (datetime.utcnow() + timedelta(days=2.0)).isoformat(),
                         }
                        ]
    def test_build_lco_imaging_request(self):
        for config in self.params:
            obs_request = obs_control.build_lco_imaging_request(config)
            assert(type(obs_request) == type(lco.LCOImagingObservationForm()))
            assert(obs_request.is_valid())
