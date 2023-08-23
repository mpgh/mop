from django.test import TestCase
from unittest import skip
from tom_targets.tests.factories import SiderealTargetFactory
#from mop.toolbox import interformetry_prediction

class TestInterferometryFunctions(TestCase):
    def setUp(self):
        self.st = SiderealTargetFactory.create()
        self.st.ra = 76.925
        self.st.dec = 24.79888
        extra_params = {'u0': 0.08858,
                        #EU0 = 0.0003,
                        'baseline_magnitude': 17.989,
                        }
        self.st.save(extras=extra_params)
        self.params = {'test_event': self.st}

    @skip("")
    def test_GAIA_toJHK(self):
        #GAIA_toJHK(G, BpRp)
        pass