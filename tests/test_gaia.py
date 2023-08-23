from django.test import TestCase
from tom_targets.models import Target,TargetExtra
from unittest import skip
from tom_targets.tests.factories import SiderealTargetFactory
from mop.brokers import gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.utils.commons import TableList
from astropy.coordinates import Angle

class TestGaia(TestCase):
    def setUp(self):
        self.st = SiderealTargetFactory.create()
        self.st.name = 'Gaia23aiy'
        self.st.ra = 241.2805
        self.st.dec = -56.4367
        self.params = {'target': self.st,
                       'radius': Angle(0.004, "deg")}

    def test_query_gaia_dr3(self):
        results = gaia.query_gaia_dr3(self.params['target'], radius=self.params['radius'])
        assert(type(results) == type(TableList([])))
        assert(len(results) > 0)
        print(results)
        for star in results:
            print(star.columns)

    def test_fetch_gaia_photometry(self):
        target = self.params['target']

        target = gaia.fetch_gaia_photometry(target)
        updated_target = Target.objects.all()[0]
        print(updated_target.extra_fields)
        print('TARGET: ',target.extra_fields)

        expected_minimum_fields = ['Gmag',
                           'RPmag',
                           'BPmag',
                           'BP-RP']

        for field in expected_minimum_fields:
            print('EXPECT field: ',field, target.extra_fields.keys())
            assert(field in target.extra_fields.keys())
            assert(target.extra_fields[field])
