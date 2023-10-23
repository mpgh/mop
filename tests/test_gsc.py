from django.test import TestCase
from tom_targets.tests.factories import SiderealTargetFactory
from mop.brokers import gsc
from astropy.table import Table, Column
import numpy as np

class TestGSC(TestCase):
    def setUp(self):
        self.st1 = SiderealTargetFactory.create()
        self.st1.name = 'Gaia23aiy'
        self.st1.ra = 265.76672  # 17:43:04.013
        self.st1.dec = -35.25895 # -35:15:32.236
        self.gsc_table = gsc.query_gsc(self.st1)

        # Example output for this test star, verified against ESO's original code
        self.aoft_table = Table([
            Column(name='FTstar', data=np.array(['S8BA011095', 'S8BA393443', 'S8BA011180', 'S8BA086534'])),
            Column(name='SC_separation', data=np.array([0.003631, 0.007347, 0.007644, 0.007724])),
            Column(name='Ksmag', data=np.array([9.750, 8.341, 10.442, 7.911])),
            Column(name='SC_Vloss', data=np.array([79.9297485, 39.9641216, 37.05225546, 36.28624988])),
            Column(name='S8BA011095_SCstrehl', data=np.array([12.703605722016004]*4)),
            Column(name='S8BA011095_FTstrehl', data=np.array([26.62470600000001, 9.067861112603005,
                                                              10.211233460147142, 10.257554843746817])),
            Column(name='S8BA011095_Gmag', data=np.array([13.698877]*4)),
            Column(name='S8BA011095_SC_separation', data=np.array([0.003631]*4)),
            Column(name='S8BA011095_Ksmag', data=np.array([10.191937926243543, 9.952387829175338,
                                                           11.92445446429255, 9.388540359135293])),
            Column(name='S8BA011095_FT_separation', data=np.array([0.0, 0.004991934667016187,
                                                                   0.004395489714960432, 0.004371325972082089]))
        ])

    def test_query_gsc(self):
        result = gsc.query_gsc(self.st1)

        # Query should return a list object with one entry, consisting of a table of non-zero length
        assert(len(result) == 1)
        assert(type(result[0]) == type(Table([])))
        assert(len(result[0]) > 0)

        # Results table should have the following columns:
        expected_column_list = ['_r', 'GSC2', 'RA_ICRS', 'DE_ICRS', 'Gmag', 'e_Gmag', 'Bmag', 'e_Bmag',
                                'Vmag', 'e_Vmag', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Ksmag', 'e_Ksmag',
                                'W1mag', 'e_W1mag', 'pmRA', 'pmDE', 'plx']
        assert(result[0].colnames == expected_column_list)

    def test_verify_Ksmag_data(self):
        initial_missing_entries = np.where(self.gsc_table['Ksmag'].mask)[0]

        amended_table = gsc.verify_Ksmag_data(self.gsc_table)

        final_missing_entries = np.where(amended_table['Ksmag'].mask)[0]

        assert(type(amended_table) == type(Table([])))
        assert(len(amended_table) == len(self.gsc_table))
        assert(len(final_missing_entries) <= len(initial_missing_entries))

    def test_Ksmag_from_phot(self):
        missing_entries = np.where(self.gsc_table['Ksmag'].mask)[0]

        test_functions = [
            (gsc.Ksmag_from_JH, 'Jmag', 'Hmag', None),
            (gsc.Ksmag_from_HW1, 'Hmag', 'W1mag', None),
            (gsc.Ksmag_from_JW1, 'Jmag', 'W1mag', None),
            (gsc.Ksmag_from_JHW1, 'Jmag', 'Hmag', 'W1mag')
                          ]

        for config in test_functions:
            amended_entries = gsc.find_missing_entries_with_phot(self.gsc_table, missing_entries, config[1], config[2],
                                                                 band3=config[3])

            amended_table = config[0](self.gsc_table, missing_entries)

            assert(type(amended_table) == type(Table([])))
            assert(len(amended_table) == len(self.gsc_table))
            assert((amended_table['Ksmag'].mask[amended_entries] == False).all())
            for x in amended_table['Ksmag'].data[amended_entries]:
                assert(type(x) == np.float64)

    def test_select_AO_stars(self):
        amended_table = gsc.select_AO_stars(self.gsc_table)

        assert('AOstar' in amended_table.colnames)
        idx = np.where(amended_table['AOstar'] > 0.0)[0]
        assert(len(idx) > 0)

    def test_select_FT_stars(self):
        amended_table = gsc.select_FT_stars(self.gsc_table)

        assert('FTstar' in amended_table.colnames)
        idx = np.where(amended_table['FTstar'] > 0.0)[0]
        assert(len(idx) > 0)

    def test_create_AOFT_table(self):
        gsc_table = gsc.select_AO_stars(self.gsc_table)
        gsc_table = gsc.select_FT_stars(gsc_table)

        idxFT = np.where(gsc_table['FTstar'] > 0.0)[0]
        idxAO = np.where(gsc_table['AOstar'] > 0.0)[0]
        ncol = 4 + 6*len(idxAO)

        AOFT_table = gsc.create_AOFT_table(gsc_table)
        assert(type(AOFT_table) == type(Table([])))
        assert(len(AOFT_table) == len(idxFT))
        assert(len(AOFT_table.colnames) == ncol)

    def test_calc_mutual_separations(self):
        gsc_table = gsc.select_AO_stars(self.gsc_table)
        gsc_table = gsc.select_FT_stars(gsc_table)

        AOFT_table = gsc.create_AOFT_table(gsc_table)

        AOFT_table = gsc.calc_mutual_separations(gsc_table, AOFT_table)

        for col in AOFT_table.colnames:
            if '_FT_separation' in col:
                assert(AOFT_table[col].sum() > 0.0)

    def test_AOstrehl(self):
        strehl1 = gsc.AOstrehl(9.0, 5.0, 'visible')
        assert(type(strehl1) == np.float64)
        assert(strehl1 < 0.5)

        strehl2 = gsc.AOstrehl(7.0, 5.0, 'ir')
        assert(type(strehl2) == np.float64)
        assert(strehl2 < 0.5)

        strehl3 = gsc.AOstrehl(7.0, 5.0, 'foo')
        assert(np.isnan(strehl3))

    def test_calc_Vloss(self):
        vloss1 = gsc.calc_Vloss(np.arange(1.0,5.0,1.0))
        assert(type(vloss1) == type(np.zeros(1)))
        assert((vloss1 < 1.0).all())

    def test_populate_AOFT_table(self):
        gsc_table = gsc.select_AO_stars(self.gsc_table)
        gsc_table = gsc.select_FT_stars(gsc_table)

        AOFT_table = gsc.create_AOFT_table(gsc_table)

        AOFT_table = gsc.populate_AOFT_table(gsc_table, AOFT_table)

        assert(AOFT_table.colnames == self.aoft_table.colnames)

        for col in self.aoft_table.colnames:
            if col in ['FTstar']:
                assert((self.aoft_table[col].data == AOFT_table[col].data).all())
            else:
                for row in range(0,len(self.aoft_table),1):
                    assert(round(self.aoft_table[col][row], 3) == round(AOFT_table[col][row],3))