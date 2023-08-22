from django.test import TestCase
from mop.brokers import ogle
from tom_targets.tests.factories import SiderealTargetFactory
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

class TestOgleBroker(TestCase):
    def setUp(self):
        self.params = {
                        'years': ['2022', '2023'],
                        'events': {'OGLE-2023-BLG-0001': ('17:53:44.28', '-28:40:04.6'),
                                   'OGLE-2023-BLG-0002': ('17:54:01.37', '-28:26:55.9')},
                        'photometry': np.array([
                                [2459795.61659, 21.078, 0.288, 5.21, 386.0],
                                [2459795.65233, 20.418, 0.179, 5.75, 401.0],
                                [2459795.68540, 20.710, 0.277, 6.47, 416.0],
                                [2459796.51954, 20.482, 0.282, 7.07, 481.0],
                                [2459796.62883, 20.621, 0.237, 5.91, 412.0],
                                [2459796.66275, 20.645, 0.223, 5.49, 420.0],
                                [2459797.53976, 20.256, 0.285, 7.47, 645.0],
                                [2459798.71855, 20.042, 0.290, 6.56, 1507.0],
                                [2459804.52573, 20.667, 0.270, 4.91, 691.0]
                        ])
        }

    def test_fetch_lens_model_parameters(self):
        broker = ogle.OGLEBroker()
        events = broker.fetch_lens_model_parameters(self.params['years'])
        assert(type(events) == type({}))
        assert(len(events) > 0)
        for key, value in events.items():
            assert(type(value) == type(()))
            assert(len(value) == 2)

    def test_ingest_ogle_events(self):
        broker = ogle.OGLEBroker()
        test_target = SiderealTargetFactory.create()
        target_list = broker.ingest_events(self.params['events'])
        assert(type(target_list) == type([]))
        for target in target_list:
            assert(type(target) == type(test_target))

    def test_read_ogle_lightcurve(self):
        test_target = self.generate_test_target()

        broker = ogle.OGLEBroker()
        photometry = broker.read_ogle_lightcurve(test_target)
        assert(type(photometry) == type(np.array([])))

    def test_ingest_ogle_photometry(self):
        test_target = self.generate_test_target()

        broker = ogle.OGLEBroker()
        status = broker.ingest_ogle_photometry(test_target, self.params['photometry'])
        assert(status == 'OK')

    def generate_test_target(self):

        test_target = SiderealTargetFactory.create()
        target_name = 'OGLE-2023-BLG-0001'
        s = SkyCoord(self.params['events'][target_name][0], self.params['events'][target_name][0],
                     unit=(u.hourangle, u.deg), frame='icrs')
        test_target.name = target_name
        test_target.ra = s.ra.deg
        test_target.dec = s.dec.deg

        return test_target
