from django.test import TestCase
from mop.toolbox import utilities
from tom_targets.tests.factories import SiderealTargetFactory
from tom_targets.models import Target
class TestUtilities(TestCase):
    def setUp(self):
        self.target = SiderealTargetFactory.create()
        self.target.name = 'test_target'
        self.target.ra = 262.71041667
        self.target.dec = -28.50847222

    def test_add_gal_coords(self):
        utilities.add_gal_coords(self.target)

        qs = Target.objects.filter(name=self.target.name)

        assert(hasattr(qs[0], 'galactic_lat'))
        assert(hasattr(qs[0], 'galactic_lng'))