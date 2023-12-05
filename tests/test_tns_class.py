from django.test import TestCase
from tom_targets.tests.factories import SiderealTargetFactory
from tom_dataproducts.models import ReducedDatum

from mop.brokers import tns

from os import getcwd, path
import numpy as np
from astropy.time import Time, TimezoneInfo
from astropy.coordinates import SkyCoord
from astropy import units as u

class TestClassSN(TestCase):
    """
    Class describing unittests for the target planet priority functions
    """
    def setUp(self):
        st1 = SiderealTargetFactory.create()
        st1.name = 'Gaia23dje'
        st1.ra = 41.41608
        st1.dec = -20.80996
        cwd = getcwd()
        lightcurve_file = path.join(cwd, 'tests/data/Gaia23dje.csv')
        # photometry = generate_test_ReducedDatums(st1, lightcurve_file, 'G')

        self.model_params = {'t': 2460284.3426,
                             't0': 2460273.62618,
                             'u0': 0.08516,
                             'te': 77.25305,
                             'Source_magnitude': 19.111,
                             'Blend_magnitude': 'nan',
                             'Baseline_magnitude': 19.111,
                             'chi2': 16.227,
                             'red_chi2': 1.803,
                             }
        self.params = {
            'target': st1,
            'lightcurve_file': lightcurve_file,
            # 'photometry': photometry,
            'Latest_data_HJD': 2460277.45285
        }

    def test_tns_response(self):

        parameters = {
            'ra' : self.params['target'].ra,
            'dec' : self.params['target'].dec,
            'radius' : 1.0,
            'units' : 'arcsec'
        }
        tns_classes = tns.Custom_TNS.fetch_tns_class(cls, parameters)
        print(tns_classes)
        self.assertEqual(tns_classes, ['SN Ic-BL'])