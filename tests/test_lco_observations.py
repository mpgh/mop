import os
from django.test import TestCase
from datetime import datetime, timedelta
from tom_observations.facilities import lco
from tom_targets.tests.factories import SiderealTargetFactory
from astropy.time import Time

from mop.toolbox import lco_observations

class TestLCOObservation(TestCase):
    def setUp(self):
        self.params = [
                    ({
                        "group_id": "TEST1",
                        "submitter": "rstreet",
                        "proposal_id": "LCO2023A-001",
                        "observation_type": "NORMAL",
                        "telescope_class": "1m0",
                        "instrument_type": "1M0-SCICAM-SINISTRO",
                        "target_name": "M31",
                        "target_type": "ICRS",
                        "ra": 121.174329,
                        "dec": -21.573309,
                        "max_airmass": 1.6,
                        "min_lunar_distance": 30.0,
                        "max_lunar_phase": 1.0,
                        "tel_code": "ogg-clma-2m0a",
                        "exposure_counts": [1],
                        "exposure_times": [30.0],
                        "filters": ["Bessell-R"],
                        "ipp": 1.05,
                        "tstart": datetime.strptime("2017-03-22 14:26:08","%Y-%m-%d %H:%M:%S"),
                        "tend": datetime.strptime("2017-03-22 14:26:08","%Y-%m-%d %H:%M:%S")
                      },
                      {
                        "submitter": "rstreet",
                        "name": "TEST1",
                        "observation_type": "NORMAL",
                        "operator": "SINGLE",
                        "ipp_value": 1.05,
                        "proposal": "LCO2023A-001",
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
                                                                "name": "M31",
                                                                "ra": 121.174329,
                                                                "dec": -21.573309,
                                                                "proper_motion_ra": 0,
                                                                "proper_motion_dec": 0,
                                                                "parallax": 0,
                                                                "epoch": 2000,
                                                                "extra_params": {}
                                                            },
                                                }],
                                    "windows": [{
                                                "start": "2017-03-22 14:26:08",
                                                "end": "2017-03-22 14:26:08"
                                                }],
                                    "observation_note": "",
                                    "acceptability_threshold": 90,
                                    "configuration_repeats": 1,
                                    "optimization_type": "TIME",
                                    "extra_params": {}
                                  }
                                ]
                      }
                    )]
    def test_build_obs_request(self):
        for config in self.params:
            obs = lco_observations.LasCumbresObservation(config[0])
            obs.build_obs_request()
            for key, value in config[1].items():
                if not key=='requests':
                    assert (obs.request[key] == config[1][key])
                else:
                    for key2,value2 in config[1]['requests'][0].items():
                        assert (obs.request[key][0][key2] == value2)