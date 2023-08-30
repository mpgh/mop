from django.test import TestCase
from unittest import skip
from tom_targets.tests.factories import SiderealTargetFactory
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target, TargetExtra
import numpy as np
from os import getcwd, path
from mop.toolbox import fittools
from astropy.time import Time, TimeDelta, TimezoneInfo
from astropy import units as u
from astropy.table import Column, Table
from pyLIMA.models import PSPL_model
from pyLIMA import event
from pyLIMA import telescopes
from datetime import datetime

class TestModelingTools(TestCase):
    def setUp(self):
        st1 = SiderealTargetFactory.create()
        st1.name = 'OGLE-2023-BLG-0348'
        st1.ra = 271.1925
        st1.dec = -28.3164
        pevent = self.generate_test_pylima_event(st1)
        cwd = getcwd()
        lightcurve_file = path.join(cwd,'tests/data/OGLE-2023-BLG-0348_phot.dat')
        # Configuration for the different telescopes, dataset_label: [Npts, baseline_mag, median_sigma, tstart, cadence]
        self.tel_configs = {
            'I': [2000, 18.0, 0.01,
                  Time('2023-08-01T00:00:00.0', format='isot'),
                  TimeDelta(0.25*u.day)],
            'G': [100, 17.5, 0.005,
                  Time('2023-08-01T00:00:00.0', format='isot'),
                  TimeDelta(1.0*u.day)]
        }
        photometry = self.generate_test_ReducedDatums(st1)
        model = self.generate_test_pylima_model(pevent)
        model_params = {
            't0': 2560000.5,
            'u0': 0.01,
            'tE': 25.0,
            'piEN': 0.01,
            'piEE': -0.01,
            'Source_magnitude': 20.0,
            'Blend_magnitude': 22.0,
            'Baseline_magnitude': 18.0,
            'Fit_covariance': np.array([[6.28180528e-02, -2.03716118e-04,  3.18785999e-02, -3.98168135e-02],
                             [-2.03716118e-04,  1.11524336e-04,  1.66302104e-03,  3.99922604e-03],
                             [ 3.18785999e-02,  1.66302104e-03,  1.37420159e-01, -2.37194137e-01],
                             [-3.98168135e-02,  3.99922604e-03, -2.37194137e-01,  1.11839732e+01]]),
            'chi2': 3000.0,
            'red_chi2': 1.02
            }
        self.params = {
            'target': st1,
            'lightcurve_file': lightcurve_file,
            'photometry': photometry,
            'pylima_event': pevent,
            'pylima_model': model,
            'model_params': model_params
            }

    def load_test_photometry(self, lightcurve_file):

        data = np.loadtxt(lightcurve_file)
        photometry = []
        for row in data:
            photometry.append([float(row[0]), float(row[1]), float(row[2])])

        datasets = {}
        datasets['I'] = np.array(photometry)

        return datasets

    def generate_test_ReducedDatums(self, target):
        """Method generates a set of ReducedDatums for different telescopes, as is held in the TOM for a
        single target"""

        data = []
        for tel_label, config in self.tel_configs.items():
            mags = np.random.normal(loc=config[1], scale=config[2], size=config[0])
            mag_errs = np.random.normal(loc=config[2], scale=config[2], size=config[0])
            for i in range(0,config[0],1):
                ts = config[3] + i*config[4]
                datum = {'magnitude': mags[i],
                         'filter': tel_label,
                         'error': mag_errs[i]
                         }
                rd, created = ReducedDatum.objects.get_or_create(
                    timestamp=ts.to_datetime(timezone=TimezoneInfo()),
                    value=datum,
                    source_name='OGLE',
                    source_location=target.name,
                    data_type='photometry',
                    target=target)

                if created:
                    rd.save()
                    data.append(rd)

        return data

    def generate_test_pylima_event(self, target):
        pevent = event.Event(ra=target.ra, dec=target.dec)
        pevent.name = target.name
        return pevent
    
    def generate_test_pylima_model(self, pevent):
        pspl = PSPL_model.PSPLmodel(pevent)
        npts = 500
        data = [
            Column(name='time', data=np.linspace(2460000.0, 2460500.0, npts)),
            Column(name='mag', data=np.random.normal(loc=18.0, scale=0.5, size=npts))
        ]
        pspl.lightcurve_magnitude = Table(data)
        return pspl

    def test_fit_pspl_omega2(self):

        datasets = self.load_test_photometry(self.params['lightcurve_file'])

        (model_params, model_lightcurve) = fittools.fit_pspl_omega2(
                self.params['target'].ra, self.params['target'].dec, datasets)

        print('RESULTS: ',model_params)
    def test_repackage_lightcurves(self):

        (datasets, ndata) = fittools.repackage_lightcurves(self.params['photometry'])

        assert(type(datasets) == type({}))
        assert(len(datasets.keys()) == len(self.tel_configs.keys()))

        test_ndata = 0
        for tel, config in self.tel_configs.items():
            test_ndata += config[0]
        assert(ndata == test_ndata)

        for dname, config in self.tel_configs.items():
            datalist = datasets[dname]
            assert(type(datalist) == type(np.array([])))
            assert(len(datalist) == config[0])

    def test_store_model_lightcurve(self):

        fittools.store_model_lightcurve(self.params['target'], self.params['pylima_model'])

        qs = ReducedDatum.objects.filter(
            source_name='MOP',
            data_type='lc_model',
            source_location=self.params['target'].name
        )

        assert(qs.count() > 0)
        
    def test_event_alive(self):
        t0_fit = Time.now() + TimeDelta(2.0*u.day)
        tE_fit = 25.0
        
        status = fittools.check_event_alive(t0_fit.jd, tE_fit)
        assert(status == True)

        t0_fit = Time.now() - TimeDelta(75.0 * u.day)
        tE_fit = 25.0

        status = fittools.check_event_alive(t0_fit.jd, tE_fit)
        assert(status == False)

    def test_store_model_parameters(self):

        fittools.store_model_parameters(
            self.params['target'],
            self.params['model_params'],
            True
        )

        qs = Target.objects.filter(name=self.params['target'].name)
        t = qs[0]

        for key, value in self.params['model_params'].items():
            if key != 'Fit_covariance':
                test_value = t.targetextra_set.get(key=key)
                assert(value == float(test_value.value))

    def test_pylima_telescopes_from_datasets(self):

        (datasets, ndata) = fittools.repackage_lightcurves(self.params['photometry'])

        tel_list = fittools.pylima_telescopes_from_datasets(datasets, emag_limit=None)

        assert(len(tel_list) == len(datasets))
        for tel in tel_list:
            assert(type(tel) == type(telescopes.Telescope()))