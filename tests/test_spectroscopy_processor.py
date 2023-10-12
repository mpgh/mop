from django.test import TestCase
from tom_targets.tests.factories import SiderealTargetFactory
from tom_dataproducts.models import DataProduct
from django.core.files.uploadedfile import SimpleUploadedFile
from os import path
from mop.processors import spectroscopy_processor

class TestSpectrumProcessor(TestCase):
    def setUp(self):
        st1 = SiderealTargetFactory.create()
        test_file = './tests/data/spectrum_sample.csv'
        self.params = {'target': st1, 'test_file': test_file}

    def test_process_data(self):
        # Load data file contents as a binary stream
        with open(self.params['test_file'], 'rb') as f:
            raw_data = f.read()
            f.close()

        # Create a DataProduct for this data
        dp = DataProduct.objects.create(
            product_id='test_product_id',
            target=self.params['target'],
            data=SimpleUploadedFile(path.basename(self.params['test_file']),
                                    raw_data)
        )

        # Instantiate the custom spectrum processor
        data_processor = spectroscopy_processor.SpectroscopyProcessor()

        # Use the custom processor to parse the data
        data = data_processor.process_data(dp)

        # Verify that the processor returns a list of tuples with three components:
        # (obs_date, datum, source_id)
        assert(type(data) == type([]))
        for dp in data:
            assert(len(dp) == 3)