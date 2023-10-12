from django.test import TestCase
from unittest import skip
from tom_targets.tests.factories import SiderealTargetFactory
from mop.toolbox import interferometry_prediction
from astropy.table import Table, Column
from mop.brokers import gaia
from astropy.coordinates import Angle
from astroquery.utils.commons import TableList

class TestInterferomeinterferometry_predictiontryFunctions(TestCase):
    def setUp(self):
        self.st = SiderealTargetFactory.create()
        self.st.ra = 274.2974
        self.st.dec = -22.3452
        extra_params = {'u0': 0.08858,
                        'u0_error': 0.0003,
                        'Baseline_magnitude': 17.989,
                        }
        self.st.save(extras=extra_params)
        neighbours = gaia.query_gaia_dr3(self.st, radius=Angle(5.0/3600.0, "deg"))
        test_star = {'Gmag': 14.0,
                     'BPRP': 0.3,
                     'u0': 0.01,
                     'u0_error': 0.01,
                     'Klens': 9.0,
                     'Kneighbours1': [12.0, 8.5, 10.7, 10.9],
                     'Kneighbours2': [13.5, 11.5, 13.5, 12.6]}
        test_star2 = {'Gmag': 8.0,
                     'BPRP': 0.3,
                     'u0': 0.01,
                     'u0_error': 0.01,
                     'Klens': 9.0,
                     'Kneighbours1': [9.5, 8.0, 10.7, 10.9]}
        self.params = {
            'test_event': self.st,
            'test_catalog': neighbours,
            'test_star': test_star,
            'test_star2': test_star2
                       }
    def test_convert_Gmag_to_JHK(self):

        (J, H, K) = interferometry_prediction.convert_Gmag_to_JHK(self.params['test_catalog'][0]['Gmag'],
                                                                 self.params['test_catalog'][0]['BP-RP'])
        for passband in [J, H, K]:
            assert(len(passband) == len(self.params['test_catalog'][0]['Gmag']))

    def test_find_companion_stars(self):

        stars_table = interferometry_prediction.find_companion_stars(self.params['test_event'],
                                                                    self.params['test_catalog'])

        assert(type(stars_table) == type(Table([])))
        assert(len(stars_table) <= len(self.params['test_catalog'].values()[0]))
        assert(len(stars_table.columns) == 4)

    def test_estimate_target_Gaia_phot_uncertainties(self):
        Gmag_error = interferometry_prediction.estimate_target_Gaia_phot_uncertainties(self.params['test_star']['Gmag'],
                                                                                      self.params['test_star']['u0'],
                                                                                      self.params['test_star']['u0_error'])
        assert(Gmag_error < 2.0)

    def test_interferometry_decision(self):
        (mode1, guide1) = interferometry_prediction.interferometry_decision(self.params['test_star']['Gmag'],
                                                                           self.params['test_star']['BPRP'],
                                                                            self.params['test_star']['Kneighbours1'])
        assert(mode1 == 'Dual Field Wide')
        assert(guide1 == 1)

        (mode2, guide2) = interferometry_prediction.interferometry_decision(self.params['test_star']['Gmag'],
                                                                           self.params['test_star']['BPRP'],
                                                                           self.params['test_star']['Kneighbours2'])
        assert (mode2 == 'No')
        assert (guide2 == 0)

        (mode2, guide2) = interferometry_prediction.interferometry_decision(self.params['test_star2']['Gmag'],
                                                                           self.params['test_star2']['BPRP'],
                                                                           self.params['test_star2']['Kneighbours1'])
        assert (mode2 == 'Single Field')
        assert (guide2 == 0)

    def test_evaluate_target_for_interferometry(self):
        interferometry_prediction.evaluate_target_for_interferometry(self.params['test_event'])
