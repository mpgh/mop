from django.test import TestCase
from tom_targets.models import Target
from unittest import skip
from tom_targets.tests.factories import SiderealTargetFactory
from mop.brokers import gsc

class TestGSC(TestCase):
    def setUp(self):
        self.st1 = SiderealTargetFactory.create()
        self.st1.name = 'Gaia23aiy'
        self.st1.ra = 265.76672  # 17:43:04.013
        self.st1.dec = -35.25895 # -35:15:32.236

    def test_query_gsc(self):
        result = gsc.query_gsc(self.st1)
