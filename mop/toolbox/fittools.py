import numpy as np
from pyLIMA import event
from pyLIMA import telescopes
from pyLIMA.fits import TRF_fit
from pyLIMA.models import PSPL_model


def chi2(params, fit):

    chi2 = np.sum(fit.residuals_LM(params)**2)
    return chi2


def flux_to_mag(flux):

    ZP_pyLIMA = 27.4
    magnitude = ZP_pyLIMA -2.5*np.log10(flux)
    return magnitude


def fit_PSPL_omega2(photometry, emag_limit = None, cores = None):

    current_event = event.Event()
    current_event.name = 'MOP_to_fit'
    filters = np.unique(photometry[:, -1])

    for ind, filt in enumerate(filters):

        if emag_limit:

            mask = (photometry[:, -1] == filt) & (np.abs(photometry[:, -2].astype(float)) < emag_limit)

        else:

            mask = (photometry[:, -1] == filt)
        lightcurve = photometry[mask, :-1].astype(float)
        # Treating all sites as ground-based without coordinates
        telescope = telescopes.Telescope(name='Tel_'+str(ind), camera_filter=filt,
                                         light_curve=lightcurve,
                                         light_curve_names=['time', 'mag', 'err_mag'],
                                         light_curve_units=['JD', 'mag', 'err_mag'])
        if len(lightcurve) > 5:
            current_event.telescopes.append(telescope)

    pspl = PSPL_model.PSPLmodel(current_event, parallax=['None', 0.])
    pspl.define_model_parameters()
    fit_tap = TRF_fit.TRFfit(pspl, loss_function='soft_l1')
    fit_tap.fit_parameters["tE"][1] = [1., 3000.]
    fit_tap.fit_parameters["u0"][1] = [0., 1.5]
    fit_tap.fit()
    cov_fit1 = fit_tap.fit_results["covariance_matrix"]
    best_fit1 = fit_tap.fit_results["best_model"]
    mag_blend_fit = flux_to_mag(fit_tap.fit_results["best_model"][4])
    mag_source_fit = flux_to_mag(fit_tap.fit_results["best_model"][3])
    mag_baseline_fit = flux_to_mag(fit_tap.fit_results["best_model"][3] + fit_tap.fit_results["best_model"][4])
    if (np.abs(best_fit1[4]) < 3. * cov_fit1[4, 4] ** 0.5) or\
            (np.abs(best_fit1[3]) < 3. * cov_fit1[3, 3] ** 0.5) or\
            (np.abs(best_fit1[2]) < 3. * cov_fit1[2, 2] ** 0.5):
        pspl = PSPL_model.PSPLmodel(current_event, parallax=['None', 0.],
                                    blend_flux_parameter='noblend')
        pspl.define_model_parameters()
        fit_tap = TRF_fit.TRFfit(pspl, loss_function='soft_l1')
        fit_tap.fit_parameters["tE"][1] = [1., 3000.]
        fit_tap.fit_parameters["u0"][1] = [0., 1.5]
        fit_tap.fit()
        mag_blend_fit = -99.
        mag_source_fit = flux_to_mag(fit_tap.fit_results["best_model"][4])
        mag_baseline_fit = flux_to_mag(fit_tap.fit_results["best_model"][3])

    epsilon_numerical_noise = 1e-6
    if np.abs(fit_tap.fit_parameters["u0"][1][0] - fit_tap.fit_results['best_model'][1]) < epsilon_numerical_noise or\
            np.abs(fit_tap.fit_parameters["u0"][1][1] - fit_tap.fit_results['best_model'][1]) < epsilon_numerical_noise or\
            np.abs(fit_tap.fit_parameters["tE"][1][1] - fit_tap.fit_results['best_model'][2]) < epsilon_numerical_noise or\
            np.abs(fit_tap.fit_parameters["tE"][1][1] - fit_tap.fit_results['best_model'][2]) < epsilon_numerical_noise:
        t0_fit = None
        u0_fit = None
        tE_fit = None
        chi2_fit = None

    t0_fit = fit_tap.fit_results["best_model"][0]
    u0_fit = fit_tap.fit_results["best_model"][1]
    tE_fit = fit_tap.fit_results["best_model"][2]
    chi2_fit = fit_tap.fit_results["chi2"]

    return [t0_fit, u0_fit, tE_fit, mag_source_fit, mag_blend_fit, mag_baseline_fit, chi2_fit]
