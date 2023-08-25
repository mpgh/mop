from django.test import TestCase
from unittest import skip
from tom_targets.tests.factories import SiderealTargetFactory
from mop.toolbox import interformetry_prediction
from astropy.table import Table, Column
from mop.brokers import gaia
from astropy.coordinates import Angle
from astroquery.utils.commons import TableList

class TestInterferometryFunctions(TestCase):
    def setUp(self):
        self.st = SiderealTargetFactory.create()
        self.st.ra = 274.2974
        self.st.dec = -22.3452
        extra_params = {'u0': 0.08858,
                        #EU0 = 0.0003,
                        'baseline_magnitude': 17.989,
                        }
        self.st.save(extras=extra_params)
        neighbours = gaia.query_gaia_dr3(self.st, radius=Angle(5.0/3600.0, "deg"))
        print(neighbours)
        print(len(neighbours))
        self.params = {
            'test_event': self.st,
            'test_catalog': neighbours
                       }

    @skip("")
    def test_GAIA_toJHK(self):
        #GAIA_toJHK(G, BpRp)
        pass

    def test_find_companion_stars(self):

        stars_table = interformetry_prediction.find_companion_stars(self.params['test_event'],
                                                                    self.params['test_catalog'])
        print(stars_table)
        print(len(stars_table))

        assert(type(stars_table) == type(Table([])))
        assert(len(stars_table) == len(self.params['test_catalog']))
        assert(len(stars_table.columns) == 4)