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

def mag_to_flux(mag):

    ZP_pyLIMA = 27.4
    flux = 10**((mag - ZP_pyLIMA) / -2.5)

    return flux

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
    # Fit configuration
    use_boundaries = True
    verbose = True

    # Initialize the new event to be fitted:
    current_event = event.Event(ra=ra, dec=dec)
    current_event.name = 'MOP_to_fit'
    if verbose: logger.info('FITTOOLS: established event')

    # Using the lightcurves stored in the TOM for this target,
    # create a list of PyLIMA telescopes, and associate them with the event:
    tel_list = pylima_telescopes_from_datasets(datasets, emag_limit=emag_limit)
    for tel in tel_list:
        current_event.telescopes.append(tel)
    if verbose: logger.info('FITTOOLS: appended ' + str(len(tel_list)) +' telescopes')

    # The above function imposes a priority order on the list of lightcurves to model,
    # so the reference dataset will always be the first one
    current_event.find_survey('Tel_0')
    current_event.check_event()

    # MODEL 1: PSPL model without parallax
    pspl = PSPL_model.PSPLmodel(current_event, parallax=['None', 0.])
    pspl.define_model_parameters()
    fit_tap = TRF_fit.TRFfit(pspl, loss_function='soft_l1')
    if verbose: logger.info('FITTOOLS: Set model 1, static PSPL')
    if use_boundaries:
        delta_t0 = 10.
        default_t0_lower = fit_tap.fit_parameters["t0"][1][0]
        default_t0_upper = fit_tap.fit_parameters["t0"][1][1]
        fit_tap.fit_parameters["t0"][1] = [default_t0_lower, default_t0_upper + delta_t0]
        fit_tap.fit_parameters["tE"][1] = [1., 3000.]
        fit_tap.fit_parameters["u0"][1] = [0.0, 2.0]
        if verbose: logger.info('FITTOOLS: model 1 fit boundaries: t0: '
                                + repr(fit_tap.fit_parameters["t0"][1])
                                + ' tE: ' + repr(fit_tap.fit_parameters["tE"][1])
                                + ' u0: ' + repr(fit_tap.fit_parameters["u0"][1]))
    fit_tap.fit()
    model1_params = gather_model_parameters(current_event, fit_tap)
    if verbose: logger.info('FITTOOLS: model 1 fitted parameters ' + repr(model1_params))

    # Evaluate the quality of the best-available model.
    # If the fitted values of key parameters are at the boundaries of then they are considered to
    # be unreliable, and the fit parameters are reset to nan
    model1_params = evaluate_model(model1_params)
    if verbose: logger.info('FITTOOLS: model 1 evaluated parameters ' + repr(model1_params))

    # By default, we accept the results of this first model fit as our best model.
    # Then we test whether the initial PSPL fit results indicate a low degree of
    # blend flux.  If so, we attempt to refit the data without blending
    best_model = model1_params
    do_noblend_model = test_quality_of_model_fit(model1_params)
    if verbose: logger.info('FITTOOLS: fit no-blend model? ' + repr(do_noblend_model))

    # MODEL 2: PSPL model without blending or parallax
    if do_noblend_model:
        pspl2 = PSPL_model.PSPLmodel(current_event, parallax=['None', 0.],
                                    blend_flux_parameter='noblend')
        pspl2.define_model_parameters()
        fit_tap2 = TRF_fit.TRFfit(pspl2, loss_function='soft_l1')
        if verbose: logger.info('FITTOOLS: Set model 2, static PSPL without blending')
        if use_boundaries:
            fit_tap2.fit_parameters["t0"][1] = [default_t0_lower, default_t0_upper + delta_t0]
            fit_tap2.fit_parameters["tE"][1] = [1., 3000.]
            fit_tap2.fit_parameters["u0"][1] = [0.0, 2.0]
            if verbose: logger.info('FITTOOLS: model 2 fit boundaries: t0: '
                                + repr(fit_tap2.fit_parameters["t0"][1])
                                + ' tE: ' + repr(fit_tap2.fit_parameters["tE"][1])
                                + ' u0: ' + repr(fit_tap2.fit_parameters["u0"][1]))
        fit_tap2.fit()
        model2_params = gather_model_parameters(current_event, fit_tap2)
        # default null as in the former implementation
        model2_params['Blend_magnitude'] = np.nan
        if verbose: logger.info('FITTOOLS: model 2 fitted parameters ' + repr(model1_params))

        # Evaluate the quality of this model
        model2_params = evaluate_model(model2_params)
        if verbose: logger.info('FITTOOLS: model 2 evaluated parameters ' + repr(model1_params))

        # Decide which fit to accept based on the fitted chi2 in each case:
        if model2_params['chi2'] <= model1_params['chi2']:
            best_model = model2_params
            if verbose: logger.info('FITTOOLS: Using model 2 as best-fit model')

    # Generate the model lightcurve timeseries with the fitted parameters
    if not np.isnan(best_model['tE']):
        model_telescope = generate_model_lightcurve(current_event,best_model)
        if verbose: logger.info('FITTOOLS: generated model lightcurve')
    else:
        model_telescope = None
        if verbose: logger.info('FITTOOLS: cannot generate model lightcurve')

    return best_model, model_telescope


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
        for dataset_id in datasets.keys():
            if name in dataset_id and dataset_id not in dataset_order:
                dataset_order.append(dataset_id)

    for name in datasets.keys():
        if name not in dataset_order:
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
                  'Fit_covariance', 'chi2', 'red_chi2',
                  'KS_test', 'AD_test', 'SW_test']
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

def gather_model_parameters(pevent, model_fit):
    """
    Function to gather the parameters of a PyLIMA fitted model into a dictionary for easier handling.
    """

    # PyLIMA model objects store the fitted values of the model parameters in the fit_results attribute,
    # which is a list of the values pertaining to the model used for the fit.  Since this model can have a
    # variable number of parameters depending on which type of model is used, we use the fit object's built-in
    # list of key indices
    param_keys = list(model_fit.fit_parameters.keys())

    model_params = {}

    for i, key in enumerate(param_keys):
        if key in ['t0' 'tE']:
            ndp = 3
        else:
            ndp = 5
        model_params[key] = np.around(model_fit.fit_results["best_model"][i], ndp)

    model_params['chi2'] = np.around(model_fit.fit_results["best_model"][-1], 3)

    # If the model did not include parallax, zero those parameters
    if 'piEN' not in param_keys:
        model_params['piEN'] = 0.0
        model_params['piEE'] = 0.0

    # Calculate the reduced chi2
    ndata = 0
    for i,tel in enumerate(pevent.telescopes):
        ndata += len(tel.lightcurve_magnitude)
    model_params['red_chi2'] = np.around(model_params['chi2'] / float(ndata - len(param_keys)),3)

    # Retrieve the flux parameters, converting from PyLIMA's key nomenculture to MOPs
    key_map = {
        'fsource_Tel_0': 'Source_magnitude',
        'fblend_Tel_0': 'Blend_magnitude'
    }

    flux_index = []
    for pylima_key,mop_key in key_map.items():
        try:
            idx = param_keys.index(pylima_key)
            model_params[mop_key] = np.around(flux_to_mag(model_fit.fit_results["best_model"][idx]), 3)
            flux_index.append(idx)
        except ValueError:
            model_params[mop_key] = np.nan

    # If the model fitted contains valid entries for both source and blend flux,
    # use these to calculate the baseline magnitude.  Otherwise, use the source magnitude
    if not np.isnan(model_params['Source_magnitude']) \
           and not np.isnan(model_params['Blend_magnitude']):
        unlensed_flux = model_fit.fit_results["best_model"][flux_index[0]] \
                            + model_fit.fit_results["best_model"][flux_index[1]]
        model_params['Baseline_magnitude'] = np.around(flux_to_mag(unlensed_flux), 3)
    else:
        model_params['Baseline_magnitude'] = model_params['Source_magnitude']
    model_params['Fit_covariance'] = model_fit.fit_results["covariance_matrix"]

    model_params['fit_parameters'] = model_fit.fit_parameters

    # Calculate fit statistics
    # The model_fit.model_residuals returns photometric and astrometric residuals as a dictionary
    # while the photometric residuals provides a list of arrays consisting of the
    # photometric residuals, photometric errors, and error_flux
    try:
        res = model_fit.model_residuals(model_fit.fit_results['best_model'])
        sw_test = stats.normal_Shapiro_Wilk(
            (np.ravel(res[0]['photometry'][0]) / np.ravel(res[1]['photometry'][0])))
        model_params['SW_test'] = np.around(sw_test[0],3)
        ad_test = stats.normal_Anderson_Darling(
            (np.ravel(res[0]['photometry'][0]) / np.ravel(res[1]['photometry'][0])))
        model_params['AD_test'] = np.around(ad_test[0],3)
        ks_test = stats.normal_Kolmogorov_Smirnov(
            (np.ravel(res[0]['photometry'][0]) / np.ravel(res[1]['photometry'][0])))
        model_params['KS_test'] = np.around(ks_test[0],3)
        model_params['chi2_dof'] = np.sum((np.ravel(res[0]['photometry'][0]) / np.ravel(res[1]['photometry'][0])) ** 2) / (
                len(np.ravel(res[0]['photometry'][0])) - 5)
    except:
        model_params['SW_test'] = np.nan
        model_params['AD_test'] = np.nan
        model_params['KS_test'] = np.nan
        model_params['chi2_dof'] = np.nan

    return model_params

def test_quality_of_model_fit(model_params):
    """Function to evaluate whether the initial model fit indicates a low degree of
    blend flux.  If so, this criterion is used to determine whether to attempt
    a second model fit without blending"""

    fit_no_blend = False

    cov_fit = model_params['Fit_covariance']

    if (np.abs(model_params['Blend_magnitude']) < 3.0 * cov_fit[4, 4] ** 0.5) or\
            (np.abs(model_params['Source_magnitude']) < 3.0 * cov_fit[3, 3] ** 0.5) or\
            (np.abs(model_params['tE']) < 3. * cov_fit[2, 2] ** 0.5):

        fit_no_blend = True

    return fit_no_blend

def evaluate_model(best_model, verbose=False):
    """Function to evaluate the overall quality of the fitted model.
    The numerical noise threshold implicitly modified the permitted minimum u0 to its value.
    """

    epsilon_numerical_noise = 1e-5
    u0_epsilon = 1e-20

    test1 = np.abs(best_model['fit_parameters']["u0"][1][0] - best_model['u0'])
    test2 = np.abs(best_model['fit_parameters']["u0"][1][1] - best_model['u0'])
    test3 = np.abs(best_model['fit_parameters']["tE"][1][0] - best_model['tE'])
    test4 = np.abs(best_model['fit_parameters']["tE"][1][1] - best_model['tE'])

    if verbose:
        logger.info('FITTOOLS Evaluating model fit:')
        logger.info('Test 1 value='+str(test1)+' criterion >'+str(u0_epsilon))
        logger.info('Test 2 value='+str(test2)+' criterion >'+str(u0_epsilon))
        logger.info('Test 3 value='+str(test3)+' criterion >'+str(epsilon_numerical_noise))
        logger.info('Test 4 value='+str(test3)+' criterion >'+str(epsilon_numerical_noise))

    # The u0 constraints have been removed here after some experimentation because
    # they were found to disallow valid fits for models with low u0 or apparent low u0
    # for models pre-peak. 
    #if test1 < u0_epsilon or \
    #    test2 < u0_epsilon or \
    #    test3 < epsilon_numerical_noise or \
    #    test4 < epsilon_numerical_noise:
    if test3 < epsilon_numerical_noise or \
        test4 < epsilon_numerical_noise:
        for key in ['t0', 'u0', 'tE', 'chi2']:
            best_model[key] = np.nan

        if verbose:
            logger.info('FITTOOLS model failed evaluation')

    return best_model

def generate_model_lightcurve(pevent, model_params):
    """Function to generate a photometric timeseries corresponding to the given model parameters"""

    pyLIMA_plots.list_of_fake_telescopes = []

    # This doesn't include parallax right now, since none of the fitted models do either yet
    pspl = PSPL_model.PSPLmodel(pevent, parallax=['None', 0.])

    params = []
    parameters = ['t0', 'u0', 'tE', 'Source_magnitude', 'Blend_magnitude']
    for key in parameters:
        value = model_params[key]
        if 'magnitude' in key:
            if not np.isnan(value):
                value = mag_to_flux(value)
            else:
                value = 0.0
        params.append(value)
    pyLIMA_parameters = pspl.compute_pyLIMA_parameters(params)

    model_telescope = pyLIMA_plots.create_telescopes_to_plot_model(pspl, pyLIMA_parameters)[0]

    flux_model = pspl.compute_the_microlensing_model(model_telescope, pyLIMA_parameters)['photometry']

    magnitude = toolbox.brightness_transformation.flux_to_magnitude(flux_model)

    model_telescope.lightcurve_magnitude["mag"] = magnitude * unit.mag

    mask = ~np.isnan(magnitude)
    model_telescope.lightcurve_magnitude = model_telescope.lightcurve_magnitude[mask]

    return model_telescope