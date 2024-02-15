from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target,TargetExtra
from astropy.time import Time
from mop.toolbox import fittools
from datetime import datetime
import json

class MicrolensingEvent(Target):
    """
    Superclass designed to consolidate data associated with a single microlensing event in a 
    form convenient for processing
    """

    def __init__(self, t):
        super().__init__(self, t.name, t.type, t.ra, t.dec, t.galactic_lng, t.galactic_lat)
        self.target = t
        self.red_data = None
        self.extras = None
        self.Last_fit = None
        self.last_observation = None
        self.existing_model = None
        self.need_to_fit = True

    def __str__(self):
        return str(self.name)

    def set_extra_params(self, qs):
        """Extracts the key, value pairs from a QuerySet of ExtraFields and sets them as
        attributes of the Event"""
        self.extras = {}
        for par in qs:
            setattr(self, par.key, par.value)
            self.extras[par.key] = par

    def set_reduced_data(self, qs):
        """Extracts the timeseries data from a QuerySet of ReducedDatums, and
        creates the necessary arrays"""

        # Store the complete set of results
        self.red_data = qs

        # Unpack the lightcurve data:
        (self.datasets, self.ndata) = fittools.repackage_lightcurves(self.red_data)

        # Extract the timestamp of the last observation
        time = [Time(i.timestamp).jd for i in self.red_data if i.data_type == 'photometry']
        self.last_observation = max(time)

        # Identify a pre-existing model lightcurve, if one is available:
        for dset in qs:
            if dset.data_type == 'lc_model':
                self.existing_model = dset
                break

    def check_need_to_fit(self):
        if self.last_observation:
            if self.Last_fit:
                if (float(self.last_observation) < float(self.Last_fit)):
                    self.need_to_fit = False

    def store_model_lightcurve(self, model):
        """Method to store in the TOM the timeseries lihgtcurve corresponding to a fitted model.
        The input is a model fit object from PyLIMA"""

        # Why is this timestamp hardwired?
        model_time = datetime.strptime('2018-06-29 08:15:27.243860', '%Y-%m-%d %H:%M:%S.%f')

        # Extract the model lightcurve timeseries from the PyLIMA fit object
        data = {
            'lc_model_time': model.lightcurve_magnitude['time'].value.tolist(),
            'lc_model_magnitude': model.lightcurve_magnitude['mag'].value.tolist()
        }

        # If there is no existing model for this target, create one
        if not self.existing_model:
            rd = ReducedDatum.objects.create(
                timestamp=model_time,
                value=data,
                source_name='MOP',
                source_location=self.target.name,
                data_type='lc_model',
                target=self.target
            )

            rd.save()
            self.existing_model = rd

        # If no prior model exists, create one
        else:
            self.existing_model.timestamp = self.existing_model.timestamp
            self.existing_model.value = self.existing_model.value
            self.existing_model.source_name = 'MOP'
            self.existing_model.source_location = self.target.name
            self.existing_model.data_type = 'lc_model'
            self.existing_model.target = self.target
            self.existing_model.defaults = {'value': data}
            self.existing_model.save()

    def store_model_parameters(self, model_params):
        """Function to store the fitted model parameters in the TOM"""

        parameters = ['Alive', 'Last_fit',
                      't0', 't0_error', 'u0', 'u0_error', 'tE', 'tE_error',
                      'piEN', 'piEN_error', 'piEE', 'piEE_error',
                      'Source_magnitude', 'Source_mag_error',
                      'Blend_magnitude', 'Blend_mag_error',
                      'Baseline_magnitude', 'Baseline_mag_error',
                      'Fit_covariance', 'chi2', 'red_chi2',
                      'KS_test', 'AD_test', 'SW_test']

        for key in parameters:
            if key == 'Fit_covariance':
                data = json.dumps(model_params['Fit_covariance'].tolist())
            else:
                data = model_params[key]
            self.extras[key].value = data
            self.extras[key].save()