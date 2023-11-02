from django.test import TestCase
from tom_targets.tests.factories import SiderealTargetFactory

from mop.toolbox import classifier_tools
from mop.management.commands import gaia_classifier

from mop.brokers import gaia as gaia_mop

from os import getcwd, path
import numpy as np
from astropy.time import Time, TimezoneInfo

class TestClassMicrolens(TestCase):
    """
    Class describing unittests for the target planet priority functions
    """
    def setUp(self):
        st1 = SiderealTargetFactory.create()
        st1.name = 'Gaia23cnu'
        st1.ra = 284.1060
        st1.dec = -18.0808
        cwd = getcwd()
        lightcurve_file = path.join(cwd, 'tests/data/Gaia23cnu.csv')
        photometry = generate_test_ReducedDatums(st1, lightcurve_file, 'G')

        self.model_params = {'t': 2460237.0,
                             't0': 2460217.09,
                             'u0': 0.34,
                             'te': 126.4,
                             'Source_magnitude': 16.92,
                             'Blend_magnitude': 16.50,
                             'Baseline_magnitude': 15.94,
                             'chi2': 6.35,
                             'red_chi2': 6135.256,
                             }
        self.params = {
            'target': st1,
            'lightcurve_file': lightcurve_file,
            'photometry': photometry,
            'Latest_data_HJD': 2460207.09
        }

    def test_check_YSO(self):
        coord = SkyCoord(ra=self.ra, dec=self.dec, unit=(u.degree, u.degree), frame='icrs')
        is_YSO = classifier_tools.check_YSO(coord)
        assert (type(is_YSO) == type(True))
        self.assertEqual(is_YSO, False)

def generate_test_ReducedDatums(target, lightcurve_file, tel_label):
    """Taken from test_fittools, by R. Street.
    Method generates a set of ReducedDatums for different telescopes, as is held in the TOM for a
    single target
    """

    data = []
    ts_jds, mags = np.loadtxt(lightcurve_file, delimiter=',', skiprows=2, usecols=(1,2))
    jd = Time(float(ts_jds), format='jd', scale='utc')

    for i in range(len(mags)):
        if(mags[i] != 'untrusted' or mags[i] != 'null'):
            datum = {
                    'magnitude': mags[i],
                    'filter': tel_label,
                    'error': gaia_mop.estimateGaiaError(mags[i])
                    }
            rd, created = ReducedDatum.objects.get_or_create(
                timestamp=jd.to_datetime(timezone=TimezoneInfo()),
                value=datum,
                source_name='Gaia',
                source_location=target.name,
                data_type='photometry',
                target=target)

            if created:
                rd.save()
                data.append(rd)

    return data