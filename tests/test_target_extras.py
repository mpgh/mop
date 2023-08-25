from django.test import TestCase
from tom_targets.models import Target,TargetExtra
from unittest import skip
from tom_targets.tests.factories import SiderealTargetFactory
import astropy.units as u
from astropy.coordinates import SkyCoord

@skip("")
class TestGaia(TestCase):
    def setUp(self):
        self.st1 = SiderealTargetFactory.create()
        self.st2, created = Target.objects.get_or_create(name='TEST2',
                                                         ra=240.5,
                                                         dec=-35.0,
                                                         type='SIDEREAL',
                                                         epoch=2000)
        self.targets = [self.st1, self.st2]

        # Test data includes two keys that are in the EXTRA_FIELDS list
        # declared in the settings.py file for this TOM
        self.test_extras = {'RPmag': 17.0, 'Extinction_G': 0.3}

    def test_add_target_extra(self):
        for key, value in self.test_extras.items():
            for target in self.targets:

                # Check to see if the extra parameter keys are in the default list
                # of the extra fields for each target
                # This test fails
                assert(key in target.extra_fields.keys())

                # Add the extra field and test to see if the field is added properly.
                # This test also fails
                target.save(extras={key: value})
                assert(target.extra_fields[key] == value)