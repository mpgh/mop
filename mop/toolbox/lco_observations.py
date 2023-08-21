from astropy import units as u
from astropy.coordinates import SkyCoord
import json
from os import path
import requests


class LasCumbresObservation():

    def __init__(self, params=None):
        self.group_id = None

        # Credentials
        self.submitter = None
        self.proposal_id = None

        # Facility & instrument
        self.telescope_class = None
        self.instrument_type = None

        # Target
        self.target_name = None
        self.target_type = 'ICRS'
        self.ra = None
        self.dec = None
        self.proper_motion_ra = None
        self.proper_motion_dec = None
        self.parallax = None
        self.epoch = 2000

        # Observation parameters
        self.exposure_times = []
        self.exposure_counts = []
        self.filters = []
        self.ipp = None
        self.obs_type = 'NORMAL'

        # Constraints
        self.max_airmass = None
        self.min_lunar_distance = None
        self.max_lunar_phase = None

        # Scheduling
        self.tstart = None
        self.tend = None

        # Assign parameters, if provided
        if params:
            for key, value in params.items():
                if hasattr(self, key):
                    if key in ['exposure_times', 'exposure_counts', 'filters']:
                        data = getattr(self, key)
                        for entry in value:
                            data.append(entry)
                        setattr(self, key, data)
                    else:
                        setattr(self, key, value)

    def build_target_dict(self):
        # Check target coordinates are in decimal degrees:
        if type(self.ra) == type(1.0):
            s = SkyCoord(self.ra, self.dec, frame='icrs', unit=(u.deg, u.deg))
        else:
            s = SkyCoord(self.ra, self.dec, frame='icrs', unit=(u.hourangle, u.deg))

        target =   {
                    'name': str(self.target_name),
                    'type': self.target_type,
                    'ra': s.ra.deg,
                    'dec': s.dec.deg,
                    'proper_motion_ra': 0,
                    'proper_motion_dec': 0,
                    'parallax': 0,
                    'epoch': 2000,
                    'extra_params': {}
                    }

        return target

    def build_constraints_dict(self):

        constraints = {}
        for key in ['max_airmass', 'min_lunar_distance', 'max_lunar_phase']:
            if getattr(self,key):
                constraints[key] = float(getattr(self, key))

        return constraints

    def build_location_dict(self):

        location = {
                    'telescope_class' : self.telescope_class,
 #                   'site':             str(self.site),
 #                   'enclosure':      str(self.enclosure)
                    }

        return location
    def build_instrument_configs(self, target, constraints):
        """Function to compose the instrument configuration dictionary for a
        set of exposures"""

        def parse_filter(f):
            filters = { 'SDSS-g': 'gp', 'SDSS-r': 'rp', 'SDSS-i': 'ip',
                       'Bessell-B': 'b', 'Bessell-V': 'v', 'Bessell-R': 'r',
                       'Cousins-Ic': 'i', 'Pan-STARRS-Z': 'zs'
                       }
            if f in filters.keys():
                return filters[f]
            else:
                raise ValueError('Unrecognized filter ('+f+') requested')

        config_list = []
        for i in range(0,len(self.exposure_times),1):
            config = {
                    'type': 'EXPOSE',
                    'instrument_type': self.instrument_type,
                    'instrument_configs': [
                        {
                            'exposure_count': int(self.exposure_counts[i]),
                            'exposure_time': float(self.exposure_times[i]),
                            'mode': 'full_frame',
                            'rotator_mode': '',
                            "extra_params": {
                                "defocus": 0,
                                "offset_ra": 0,
                                "offset_dec": 0
                            },
                            'optical_elements': {
                                'filter': parse_filter(self.filters[i])
                            }
                        }
                    ],
                    'acquisition_config': {
                        'mode': 'OFF',
                        'extra_params': {}
                    },
                    'guiding_config': {
                        'mode': 'ON',
                        'optional': 'true',
                        'extra_params': {}
                    },
                    'target': target,
                    'constraints': constraints,
                }
            config_list.append(config)

        return config_list

    def build_obs_request(self):

        request_group = {
                        "submitter": self.submitter,
                        "name": self.group_id,
                        "observation_type": self.obs_type,
                        "operator": "SINGLE",
                        "ipp_value": float(self.ipp),
                        "proposal": self.proposal_id,
                         }

        target = self.build_target_dict()
        location = self.build_location_dict()
        constraints = self.build_constraints_dict()
        inst_config_list = self.build_instrument_configs(target, constraints)
        windows = {'start': self.tstart.strftime("%Y-%m-%d %H:%M:%S"),
                   'end': self.tend.strftime("%Y-%m-%d %H:%M:%S")}

        request_group['requests'] = [{
                                        "location": location,
                                        "configurations": inst_config_list,
                                        "windows": [windows],
                                        "observation_note": "",
                                        "acceptability_threshold": 90,
                                        "configuration_repeats": 1,
                                        "optimization_type": "TIME",
                                        "extra_params": {}
                                    }]

        self.request = request_group

    def submit(self, credentials):
        """Function to communicate with various APIs of the LCO network.
        ur should be a user request in the form of a Python dictionary,
        while end_point is the URL string which
        should be concatenated to the observe portal path to complete the URL.
        Accepted methods are:
            POST
        """
        PORTAL_URL = 'https://observe.lco.global/api/requestgroups/'

        headers = {'Authorization': 'Token ' + credentials['lco_token']}

        response = requests.post(PORTAL_URL, headers=headers, json=self.request).json()

        return response
