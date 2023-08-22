from django.test import TestCase
from mop.brokers import ogle
from tom_targets.tests.factories import SiderealTargetFactory

class TestOgleBroker(TestCase):
    def setUp(self):
        self.params = {
                        'years': ['2022', '2023'],
                        'events': {'OGLE-2023-BLG-0001': ('17:53:44.28', '-28:40:04.6'),
                                   'OGLE-2023-BLG-0002': ('17:54:01.37', '-28:26:55.9')}
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