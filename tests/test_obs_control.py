import os
from django.test import TestCase
from unittest import skip
from datetime import datetime, timedelta
from tom_observations.facilities import lco
from tom_targets.tests.factories import SiderealTargetFactory
from astropy.time import Time

from mop.toolbox import obs_control, lco_observations

class TestObsConfig(TestCase):
    def setUp(self):
        # Input dictionary format based on example in TOM Toolkit tom_observations/tests/facilities/test_lco.py
        time_start = datetime.utcnow() + timedelta(days=1.0)
        time_end = time_start + timedelta(days=1.0)
        self.st = SiderealTargetFactory.create()
        self.st.name = 'Gaia23bvd'
        self.st.ra = 285.8068
        self.st.dec = -30.4068
        self.params = [
                        (
                            [{
                                "group_id": "TEST1",
                                "observation_type": "NORMAL",
                                "telescope_class": "1m0",
                                "instrument_type": "1M0-SCICAM-SINISTRO",
                                "target": self.st,
                                "target_type": "ICRS",
                                "max_airmass": 1.6,
                                "min_lunar_distance": 30.0,
                                "max_lunar_phase": 1.0,
                                "tel_code": "ogg-clma-2m0a",
                                "exposure_counts": [1],
                                "exposure_times": [30.0],
                                "filters": ["Bessell-R"],
                                "ipp": 1.05,
                                "tstart": time_start,
                                "tend": time_end
                            }],
                            [{
                                "submitter": os.getenv("LCO_USERNAME"),
                                "name": "TEST1",
                                "observation_type": "NORMAL",
                                "operator": "SINGLE",
                                "ipp_value": 1.05,
                                "proposal": os.getenv('LCO_PROPOSAL_ID'),
                                "requests": [
                                    {"location": {
                                        'telescope_class': '1m0',
                                    },
                                        "configurations": [{
                                            "type": "EXPOSE",
                                            "instrument_type": "1M0-SCICAM-SINISTRO",
                                            "instrument_configs": [{
                                                "exposure_time": 30.0,
                                                "exposure_count": 1,
                                                "mode": "full_frame",
                                                "rotator_mode": "",
                                                "optical_elements": {
                                                    "filter": "r"
                                                },
                                                "extra_params": {
                                                    "defocus": 0,
                                                    "offset_ra": 0,
                                                    "offset_dec": 0
                                                }
                                            }],
                                            "acquisition_config": {
                                                "mode": "OFF",
                                                "extra_params": {}
                                            },
                                            "guiding_config": {
                                                "mode": "ON",
                                                "optional": "true",
                                                "extra_params": {}
                                            },
                                            "constraints": {
                                                "max_airmass": 1.6,
                                                "min_lunar_distance": 30.0,
                                                "max_lunar_phase": 1.0
                                            },
                                            "target": {
                                                "type": "ICRS",
                                                "name": self.st.name,
                                                "ra": self.st.ra,
                                                "dec": self.st.dec,
                                                "proper_motion_ra": 0,
                                                "proper_motion_dec": 0,
                                                "parallax": 0,
                                                "epoch": 2000,
                                                "extra_params": {}
                                            },
                                        }],
                                        "windows": [{
                                            "start": time_start.strftime("%Y-%m-%d %H:%M:%S"),
                                            "end": time_end.strftime("%Y-%m-%d %H:%M:%S")
                                        }],
                                        "observation_note": "",
                                        "acceptability_threshold": 90,
                                        "configuration_repeats": 1,
                                        "optimization_type": "TIME",
                                        "extra_params": {}
                                    }
                                ]
                            }]
                        )
                        ]
    def test_build_lco_imaging_request(self):
        for config in self.params:
            obs_list = obs_control.build_lco_imaging_request(config[0])
            for obs in obs_list:
                assert(type(obs) == type(lco_observations.LasCumbresObservation()))
                for expected_obs in config[1]:
                    for key, value in expected_obs.items():
                        if not key=='requests':
                            assert (obs.request[key] == expected_obs[key])
                        else:
                            for key2,value2 in expected_obs['requests'][0].items():
                                assert (obs.request[key][0][key2] == value2)

    @skip("Uncomment this to submit live test observations")
    def test_submit_observations(self):
        for config in self.params:
            obs = obs_control.build_lco_imaging_request(config[0])
            obs_control.submit_lco_obs_request(obs, self.st)


class TestVisibility(TestCase):
    def setUp(self):
        self.st1 = SiderealTargetFactory.create()
        self.st1.ra = 241.2806       # 16:05:07.3
        self.st1.dec = -56.4367      # -56:26:12.26
        self.st2 = SiderealTargetFactory.create()
        self.st2.ra = 255.0          # 17:00:00.0
        self.st2.dec = 10.0          # 10:00:00.0
        self.params = [{
                        'target': self.st1,
                        'timenow': Time('2023-08-20').decimalyear,
                        'result': True
                       },{
                        'target': self.st2,
                        'timenow': Time('2023-01-01').decimalyear,
                        'result': False
                       }]

    def test_check_visibility(self):
        for test_case in self.params:
            visible = obs_control.check_visibility(test_case['target'],
                                                   test_case['timenow'])
            assert(visible == test_case['result'])
            