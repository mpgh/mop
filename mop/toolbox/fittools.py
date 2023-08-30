import numpy as np
from pyLIMA import event
from pyLIMA import telescopes
from pyLIMA import toolbox
from pyLIMA.fits import TRF_fit
from pyLIMA.fits import stats
from pyLIMA.models import PSPL_model
from pyLIMA.outputs import pyLIMA_plots
from astropy import units as unit
from astropy.time import Time
import logging
from datetime import datetime
from tom_dataproducts.models import ReducedDatum
import json

logger = logging.getLogger(__name__)

def chi2(params, fit):

    chi2 = np.sum(fit.residuals_LM(params)**2)
    return chi2


def flux_to_mag(flux):

    ZP_pyLIMA = 27.4
    magnitude = ZP_pyLIMA - 2.5 * np.log10(flux)
    return magnitude


def fit_pspl_omega2(ra, dec, datasets, emag_limit=None):
    """
    Fit photometry using pyLIMAv1.9 with a static PSPL TRF fit
    checking if blend is constrained, if so using a soft_l1 loss function
    if not refit unblended model and report parameter

    Parameters
    ----------
    ra : float, Right Ascension in degrees
    dec : float, Declination in degrees
    photometry : array containing all telescope passband light curves
    emag_limit : array, limit on the error

    Returns
    -------
    to_return : list of arrays containing fit parameters, model_telescope and cost function
    """

    # Initialize the new event to be fitted:
    current_event = event.Event(ra=ra, dec=dec)
    current_event.name = 'MOP_to_fit'

    # Using the lightcurves stored in the TOM for this target,
    # create a list of PyLIMA telescopes, and associate them with the event:
    tel_list = pylima_telescopes_from_datasets(datasets, emag_limit=emag_limit)
    for tel in tel_list:
        current_event.telescopes.append(tel)

    # The above function imposes a priority order on the list of lightcurves to model,
    # so the reference dataset will always be the first one
    current_event.find_survey('Tel_0')
    current_event.check_event()

    # MODEL 1: PSPL model without parallax
    pspl = PSPL_model.PSPLmodel(current_event, parallax=['None', 0.])
    pspl.define_model_parameters()
    fit_tap = TRF_fit.TRFfit(pspl, loss_function='soft_l1')
    #delta_t0 = 10.
    #default_t0_lower = fit_tap.fit_parameters["t0"][1][0]
    #default_t0_upper = fit_tap.fit_parameters["t0"][1][1]
    #fit_tap.fit_parameters["t0"][1] = [default_t0_lower, default_t0_upper + delta_t0]
    #fit_tap.fit_parameters["tE"][1] = [1., 3000.]
    #fit_tap.fit_parameters["u0"][1] = [0., 2.0]
    fit_tap.fit()
    cov_fit1 = fit_tap.fit_results["covariance_matrix"]
    best_fit1 = fit_tap.fit_results["best_model"]
    mag_blend_fit = flux_to_mag(best_fit1[4])
    mag_source_fit = flux_to_mag(best_fit1[3])
    mag_baseline_fit = flux_to_mag(best_fit1[3] + best_fit1[4])
    fit_type = 'chi2'
    print('STAGE 1: best fit', best_fit1)

    # MODEL 2: If the initial PSPL fit results indicate a low degree of
    # blend flux, attempt to refit the data without blending
    if (np.abs(best_fit1[4]) < 3. * cov_fit1[4, 4] ** 0.5) or\
            (np.abs(best_fit1[3]) < 3. * cov_fit1[3, 3] ** 0.5) or\
            (np.abs(best_fit1[2]) < 3. * cov_fit1[2, 2] ** 0.5):
        print('Got to clause 1')
        pspl = PSPL_model.PSPLmodel(current_event, parallax=['None', 0.],
                                    blend_flux_parameter='noblend')
        pspl.define_model_parameters()
        fit_tap = TRF_fit.TRFfit(pspl, loss_function='soft_l1')
        fit_tap.fit()
        #fit_tap.fit_parameters["t0"][1] = [default_t0_lower, default_t0_upper + delta_t0]
        #fit_tap.fit_parameters["tE"][1] = [1., 3000.]
        #fit_tap.fit_parameters["u0"][1] = [0., 2.0]
        # default null as in the former implementation
        mag_blend_fit = np.nan
        mag_source_fit = flux_to_mag(fit_tap.fit_results["best_model"][3])
        mag_baseline_fit = flux_to_mag(fit_tap.fit_results["best_model"][3])
        fit_type = 'soft_l1'
    # This appears to be repeating the initial PSPL fit, so commenting out
#    else:
#        print('Got to clause 2')
#        pspl = PSPL_model.PSPLmodel(current_event, parallax=['None', 0.])
#        pspl.define_model_parameters()
#        fit_tap = TRF_fit.TRFfit(pspl)
#        fit_tap.fit_parameters["t0"][1] = [default_t0_lower, default_t0_upper + delta_t0]
#        fit_tap.fit_parameters["tE"][1] = [1., 3000.]
#        fit_tap.fit_parameters["u0"][1] = [0., 2.0]
#        fit_tap.fit()
#        mag_blend_fit = flux_to_mag(fit_tap.fit_results["best_model"][4])
#        mag_source_fit = flux_to_mag(fit_tap.fit_results["best_model"][3])
#        mag_baseline_fit = flux_to_mag(fit_tap.fit_results["best_model"][3] + fit_tap.fit_results["best_model"][4])
#        fit_type = 'chi2'

    # Store the fitted model parameter values
    # the numerical noise threshold implicitly modified the permitted minimum u0 to its value
    model_params = {}
    epsilon_numerical_noise = 1e-5
    if np.abs(fit_tap.fit_parameters["u0"][1][0] - fit_tap.fit_results['best_model'][1]) < epsilon_numerical_noise or\
            np.abs(fit_tap.fit_parameters["u0"][1][1] - fit_tap.fit_results['best_model'][1]) < epsilon_numerical_noise or\
            np.abs(fit_tap.fit_parameters["tE"][1][0] - fit_tap.fit_results['best_model'][2]) < epsilon_numerical_noise or\
            np.abs(fit_tap.fit_parameters["tE"][1][1] - fit_tap.fit_results['best_model'][2]) < epsilon_numerical_noise:
        for key in ['t0', 'u0', 'tE', 'chi2']:
            model_params[key] = np.nan
    else:
        model_params['t0'] = np.around(fit_tap.fit_results["best_model"][0], 3)
        model_params['u0'] = np.around(fit_tap.fit_results["best_model"][1], 5)
        model_params['tE'] = np.around(fit_tap.fit_results["best_model"][2], 3)
        model_params['chi2'] = np.around(fit_tap.fit_results["chi2"], 3)
    model_params['piEN'] = 0.0
    model_params['piEE'] = 0.0
    ndata = 0
    for name, lc in datasets.items():
        ndata += len(lc)
    model_params['red_chi2'] = model_params['chi2'] / float(ndata - 5)
    model_params['Source_magnitude'] = np.around(mag_source_fit, 3)
    model_params['Blend_magnitude'] = np.around(mag_blend_fit, 3)
    model_params['Baseline_magnitude'] = np.around(mag_baseline_fit, 3)
    model_params['Fit_covariance'] = fit_tap.fit_results["covariance_matrix"]

    # Generate the model lightcurve timeseries with the fitted parameters
    pyLIMA_plots.list_of_fake_telescopes = []
    pyLIMA_parameters = pspl.compute_pyLIMA_parameters(fit_tap.fit_results["best_model"])
    model_telescope = pyLIMA_plots.create_telescopes_to_plot_model(pspl, pyLIMA_parameters)[0]
    flux_model = pspl.compute_the_microlensing_model(model_telescope, pyLIMA_parameters)['photometry']
    magnitude = toolbox.brightness_transformation.flux_to_magnitude(flux_model)
    model_telescope.lightcurve_magnitude["mag"] = magnitude * unit.mag
    mask = ~np.isnan(magnitude)
    model_telescope.lightcurve_magnitude = model_telescope.lightcurve_magnitude[mask]

    # Calculate photometric statistics
    try:
        res = fit_tap.model_residuals(fit_tap.fit_results['best_model'])
        model_params['SW_test'] = stats.normal_Shapiro_Wilk((np.ravel(res[0]['photometry']) / np.ravel(res[1]['photometry'])))
        model_params['AD_test'] = stats.normal_Anderson_Darling((np.ravel(res[0]['photometry']) / np.ravel(res[1]['photometry'])))
        model_params['KS_test'] = stats.normal_Kolmogorov_Smirnov((np.ravel(res[0]['photometry']) / np.ravel(res[1]['photometry'])))
        model_params['chi2_dof'] = np.sum((np.ravel(res[0]['photometry']) / np.ravel(res[1]['photometry'])) ** 2) / (
                len(np.ravel(res[0]['photometry'])) - 5)
    except:
        model_params['SW_test'] = np.nan
        model_params['AD_test'] = np.nan
        model_params['KS_test'] = np.nan
        model_params['chi2_dof'] = np.nan

    # chi2_fit is only an actual chi-squared if the loss function was deactivated
    # [np.around(t0_fit, 3), np.around(u0_fit, 5), np.around(tE_fit, 3), np.around(piEN_fit, 5),
    #np.around(piEE_fit, 5),
    #np.around(mag_source_fit, 3), np.around(mag_blend_fit, 3), np.around(mag_baseline_fit, 3),
    #fit_tap.fit_results["covariance_matrix"], model_telescope, np.around(chi2_fit, 3), chi2_dof,
    #sw_test, ad_test, ks_test]

    return model_params, model_telescope


def repackage_lightcurves(qs):
    """Function to sort through a QuerySet of the ReducedDatums for a given event and repackage the data as a
     dictionary of individual lightcurves in PyLIMA-compatible format for different facilities.
     Note that not all of the QuerySet of ReducedDatums may be photometry, so some sorting is required.
     """

    datasets = {}

    for rd in qs:
        if rd.data_type == 'photometry':
            # Identify different lightcurves from the filter label given
            passband = rd.value['filter']
            if passband in datasets.keys():
                lc = datasets[passband]
            else:
                lc = []

            # Append the datapoint to the corresponding dataset
            try:
                lc.append([Time(rd.timestamp).jd, rd.value['magnitude'], rd.value['error']])
            except:
                lc.append([Time(rd.timestamp).jd, rd.value['magnitude'], 1.0])
            datasets[passband] = lc

    # Count the total number of datapoints available, and convert the
    # accumulated lightcurves into numpy arrays:
    ndata = 0
    for passband, lc in datasets.items():
        ndata += len(lc)
        datasets[passband] = np.array(lc)

    return datasets, ndata

def pylima_telescopes_from_datasets(datasets, emag_limit=None):
    """Function to convert the dictionary of datasets retrieved from MOP of the lightcurves for this object,
    and convert them into PyLIMA Telescope objects.
    This function returns a list of Telescope objects containing the lightcurve data, applying an
    order of preference, so that prioritized datasets occur at the start of the list.
    """

    # Sort the available datasets into order, giving preference to main survey datasets
    priority_order = ['I', 'ip', 'i_ZTF', 'r_ZTF', 'R', 'g_ZTF', 'gp', 'G']

    dataset_order = []
    for name in priority_order:
        if name in datasets.keys():
            dataset_order.append(name)

    for name in datasets.keys():
        if name not in priority_order:
            dataset_order.append(name)

    # Loop over all available datasets and create a telescope object for each one
    tel_list = []
    for idx, name in enumerate(dataset_order):
        photometry = datasets[name]

        # Enabling optional filtering for datapoints of low photometric precision
        if emag_limit:

            mask = (np.abs(photometry[:, -2].astype(float)) < emag_limit)

        else:

            mask = (np.abs(photometry[:, -2].astype(float)) < 99.0)

        lightcurve = photometry[mask].astype(float)

        # Treating all sites as ground-based without coordinates
        tel = telescopes.Telescope(name='Tel_'+str(idx), camera_filter=name,
                                         light_curve=photometry[mask],
                                         light_curve_names=['time', 'mag', 'err_mag'],
                                         light_curve_units=['JD', 'mag', 'err_mag'])
        tel_list.append(tel)

    return tel_list

def store_model_lightcurve(target, model):
    """Function to store in the TOM the timeseries lihgtcurve corresponding to a fitted model.
    The input is a model fit object from PyLIMA"""

    # Why is this timestamp hardwired?
    model_time = datetime.strptime('2018-06-29 08:15:27.243860', '%Y-%m-%d %H:%M:%S.%f')

    # Extract the model lightcurve timeseries from the PyLIMA fit object
    data = {
        'lc_model_time': model.lightcurve_magnitude['time'].value.tolist(),
        'lc_model_magnitude': model.lightcurve_magnitude['mag'].value.tolist()
        }

    # Search MOP to see if an existing model is already stored
    existing_model = ReducedDatum.objects.filter(source_name='MOP',data_type='lc_model',
                                                 timestamp=model_time, source_location=target.name)
    logger.info('FIT: Searched for existing models '+repr(existing_model))

    # If there is an existing model for this target, update it
    if existing_model.count() == 0:
        rd, created = ReducedDatum.objects.get_or_create(
                    timestamp=model_time,
                    value=data,
                    source_name='MOP',
                    source_location=target.name,
                    data_type='lc_model',
                    target=target
        )

        rd.save()

    # If no prior model exists, create one
    else:
        rd, created = ReducedDatum.objects.update_or_create(
                    timestamp=existing_model[0].timestamp,
                    value=existing_model[0].value,
                    source_name='MOP',
                    source_location=target.name,
                    data_type='lc_model',
                    target=target,
                    defaults={'value':data}
        )

        rd.save()

def store_model_parameters(target, model_params, event_alive):
    """Function to store the fitted model parameters in the TOM"""

    last_fit = Time(datetime.utcnow()).jd

    extras = {'Alive': event_alive, 'Last_fit': last_fit}

    parameters = ['t0', 'u0', 'tE', 'piEN', 'piEE',
                  'Source_magnitude', 'Blend_magnitude', 'Baseline_magnitude',
                  'Fit_covariance', 'chi2', 'red_chi2']
    for key in parameters:
        if key == 'Fit_covariance':
            data = json.dumps(model_params['Fit_covariance'].tolist())
        else:
            data = model_params[key]
        extras[key] = data

    target.save(extras=extras)

def check_event_alive(t0_fit, tE_fit):
    """Function to evaluate whether or not an event is still actively going on, based on the current time
    relative to the model fit t0 and tE"""

    time_now = Time(datetime.now()).jd

    how_many_tE = (time_now - t0_fit) / tE_fit

    if how_many_tE > 2:

        alive = False

    else:

        alive = True

    return alive