import mimetypes


from astropy import units
from astropy.io import ascii
from astropy.time import Time, TimezoneInfo


from tom_dataproducts.data_processor import DataProcessor
from tom_dataproducts.exceptions import InvalidFileFormatException
from mop.toolbox import TAP

# This is a custom processor made for AWS S3 bucket compliance

class PhotometryProcessor(DataProcessor):

    def process_data(self, data_product):
        print('IN PROCESSOR with ',data_product)

        mimetype = mimetypes.guess_type(data_product.data.name)[0]
        if mimetype in self.PLAINTEXT_MIMETYPES:
            photometry = self._process_photometry_from_plaintext(data_product)
            return [(datum.pop('timestamp'), datum, datum['filter']) for datum in photometry]
        else:
            raise InvalidFileFormatException('Unsupported file type')

    def _process_photometry_from_plaintext(self, data_product):

        from django.core.files.storage import default_storage
        
        photometry = []

        data_aws = default_storage.open(data_product.data.name, 'r')
        data = ascii.read(data_aws.read(),
                          names=['time', 'filter', 'magnitude', 'error'])

        if len(data) < 1:
            raise InvalidFileFormatException('Empty table or invalid file type')

        for datum in data:
            if float(datum['time'])>2400000.5:
                time = Time(float(datum['time']), format='jd')
            else:
                time = Time(float(datum['time']), format='mjd')
            utc = TimezoneInfo(utc_offset=0*units.hour)
            time.format = 'datetime'
            value = {
                'timestamp': time.to_datetime(timezone=utc),
                'magnitude': datum['magnitude'],
                'filter': datum['filter'],
                'error': datum['error']
            }
            photometry.append(value)

        return photometry
