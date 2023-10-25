from django.test import TestCase
from tom_targets.tests.factories import SiderealTargetFactory
from mop.toolbox import TAP_priority

class TestTAPPriorities(TestCase):
    """
    Class describing unittests for the target priority functions
    """
    def setUp(self):
        self.target = SiderealTargetFactory.create()
        self.target.name = 'Gaia23bvd'
        self.target.ra = 285.8068
        self.target.dec = -30.4068
        self.model_params = {'t': 2462000.0,
                             'te': 30.0,
                             'u0': 0.01,
                             't0': 2462010.0}

    def test_psi_derivatives_squared(self):

        result = TAP_priority.psi_derivatives_squared(self.model_params['t'],
                                                     self.model_params['te'],
                                                     self.model_params['u0'],
                                                     self.model_params['t0'])

        assert(type(result) == type([]))
        assert(len(result) == 3)
        for entry in result:
            assert(type(entry) == type(1.0))

