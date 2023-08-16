from django.core.management.base import BaseCommand
from tom_targets.models import Target

import datetime
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
class Command(BaseCommand):

    help = 'Add a set of targets to an empty test TOM for testing purposes'

    def handle(self, *args, **options):

        utc_now = datetime.datetime.utcnow()
        jd_now = Time(utc_now, scale='utc')
        covar_array = [[190405.95634386834, -11161.682226967514, 7081585.9694203865, -528121.8073966689, -204845.65284342837, -105956.19650825627, 704169.7264599944], [-11161.682368628773, 656.5436407077727, -402768.7571275531, 31756.590911363186, 12475.519219544338, 6027.565954625822, -40058.79145911312], [7081586.910229315, -402768.80395378673, 335315189.41991, -14982086.244333094, -4888519.729533412, -5009751.425231723, 33291455.641387586], [-528121.85812985, 31756.59349709558, -14982086.646950368, 1766742.2830257271, 745056.7272725725, 224638.47389013524, -1493085.7174158567], [-204845.66759474133, 12475.519960507743, -4888519.774401485, 745056.7233792959, 324014.107430808, 73420.3622829645, -488041.2248857036], [-105956.55440581981, 6027.58621524388, -5009767.67789246, 224639.19715427817, 73420.60004629717, 74848.65558735882, -497394.1897986175], [704172.104996, -40058.92610833448, 33291563.654593002, -1493090.524144607, -488042.8050322431, -497394.1899307398, 3305350.4393726843]]

        target_list = {
                    'Gaia23aiy': {'ra': '16:05:07.332', 'dec': '-56:26:12.26',
                                  'Alive': True, 'Classification': 'Microlensing PSPL',
                                  'extras': {'t0': 2459968.044, 'u0': 0.10824, 'tE': 108.239,
                                           'piEN': -0.12749, 'piEE': 0.29117, 'Source_magnitude': 19.534,
                                           'Baseline_magnitude': 20.106, 'Latest_data_UTC': utc_now,
                                           'Latest_data_HJD': jd_now.jd, 'Fit_covariance': covar_array}},
                    'Gaia22duy': {'ra': '18:41:09.250', 'dec': '-10:23:47.04',
                                  'Alive': True, 'Classification': 'Microlensing PSPL',
                                  'extras': {'t0': 2459894.825, 'u0': 1.80088, 'tE': 45.81,
                                           'piEN': -0.00038, 'piEE': 0.38285, 'Source_magnitude': 14.12,
                                           'Baseline_magnitude': 17.533, 'Latest_data_UTC': utc_now,
                                           'Latest_data_HJD': jd_now.jd, 'Fit_covariance': covar_array}},
                    'Gaia21ccu': {'ra': '09:31:50.899', 'dec': '-55:17:57.23',
                                  'Alive': True, 'Classification': 'Microlensing PSPL',
                                  'extras': {'t0': 2459894.825, 'u0': 2.0, 'tE': 1.4,
                                           'piEN': 2.0, 'piEE': -1.9, 'Source_magnitude': 14.1,
                                           'Baseline_magnitude': 14.292, 'Latest_data_UTC': utc_now,
                                           'Latest_data_HJD': jd_now.jd, 'Fit_covariance': covar_array}},
                    }

        for target_name, params in target_list.items():
            s = SkyCoord(params['ra'], params['dec'], frame='icrs', unit=(u.hourangle, u.deg))
            target, created = Target.objects.get_or_create(name=target_name, ra=s.ra.deg,
                                                           dec=s.ra.deg, type='SIDEREAL', epoch=2000)
            target.save(extras=params['extras'])
