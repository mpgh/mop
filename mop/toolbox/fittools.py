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
    magnitude = ZP_pyLIMA - 2.5 * np.log10(flux)
    return magnitude


def fit_PSPL_omega2(ra, dec, photometry, emag_limit=None):
    # Default filter sequence for establishing survey dataset
    filters_order = ['I', 'ip', 'i_ZTF', 'r_ZTF', 'R', 'g_ZTF', 'gp', 'G']

    filters = np.unique(photometry[:, -1])
    order = []
    for fil in filters_order:

        mask = np.where(filters == fil)[0]
        if len(mask) != 0:
            order += mask.tolist()

    for fil in filters:

        if fil not in filters_order:
            mask = np.where(filters == fil)[0]
            if len(mask) != 0:
                order += mask.tolist()

    filters = filters[order]
    current_event = event.Event()
    current_event.name = 'MOP_to_fit'
    filters = np.unique(photometry[:, -1])
    lightcurve = []
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
    delta_t0 = 10.
    default_t0_lower = fit_tap.fit_parameters["t0"][1][0]
    default_t0_upper = fit_tap.fit_parameters["t0"][1][1]
    fit_tap.fit_parameters["t0"][1] = [default_t0_lower, default_t0_upper + delta_t0]
    fit_tap.fit_parameters["tE"][1] = [1., 3000.]
    fit_tap.fit_parameters["u0"][1] = [0., 2.0]
    fit_tap.fit()
    cov_fit1 = fit_tap.fit_results["covariance_matrix"]
    best_fit1 = fit_tap.fit_results["best_model"]
    if (np.abs(best_fit1[4]) < 3. * cov_fit1[4, 4] ** 0.5) or\
            (np.abs(best_fit1[3]) < 3. * cov_fit1[3, 3] ** 0.5) or\
            (np.abs(best_fit1[2]) < 3. * cov_fit1[2, 2] ** 0.5):
        pspl = PSPL_model.PSPLmodel(current_event, parallax=['None', 0.],
                                    blend_flux_parameter='noblend')
        pspl.define_model_parameters()
        fit_tap = TRF_fit.TRFfit(pspl, loss_function='soft_l1')
        fit_tap.fit_parameters["t0"][1] = [default_t0_lower, default_t0_upper + delta_t0]
        fit_tap.fit_parameters["tE"][1] = [1., 3000.]
        fit_tap.fit_parameters["u0"][1] = [0., 2.0]
        fit_tap.fit()
        # default null as in the former implementation
        mag_blend_fit = "null"
        mag_source_fit = flux_to_mag(fit_tap.fit_results["best_model"][3])
        mag_baseline_fit = flux_to_mag(fit_tap.fit_results["best_model"][3])
    else:
        pspl = PSPL_model.PSPLmodel(current_event, parallax=['None', 0.],
                                    blend_flux_parameter='noblend')
        pspl.define_model_parameters()
        fit_tap = TRF_fit.TRFfit(pspl, loss_function='soft_l1')
        fit_tap.fit_parameters["t0"][1] = [default_t0_lower, default_t0_upper + delta_t0]
        fit_tap.fit_parameters["tE"][1] = [1., 3000.]
        fit_tap.fit_parameters["u0"][1] = [0., 2.0]
        fit_tap.fit()
        mag_blend_fit = flux_to_mag(fit_tap.fit_results["best_model"][4])
        mag_source_fit = flux_to_mag(fit_tap.fit_results["best_model"][3])
        mag_baseline_fit = flux_to_mag(fit_tap.fit_results["best_model"][3] + fit_tap.fit_results["best_model"][4])

    epsilon_numerical_noise = 1e-6
    if np.abs(fit_tap.fit_parameters["u0"][1][0] - fit_tap.fit_results['best_model'][1]) < epsilon_numerical_noise or\
            np.abs(fit_tap.fit_parameters["u0"][1][1] - fit_tap.fit_results['best_model'][1]) < epsilon_numerical_noise or\
            np.abs(fit_tap.fit_parameters["tE"][1][1] - fit_tap.fit_results['best_model'][2]) < epsilon_numerical_noise or\
            np.abs(fit_tap.fit_parameters["tE"][1][1] - fit_tap.fit_results['best_model'][2]) < epsilon_numerical_noise:
        t0_fit = np.nan
        u0_fit = np.nan
        tE_fit = np.nan
        chi2_fit = np.nan
    else:
        t0_fit = fit_tap.fit_results["best_model"][0]
        u0_fit = fit_tap.fit_results["best_model"][1]
        tE_fit = fit_tap.fit_results["best_model"][2]
        chi2_fit = fit_tap.fit_results["chi2"]
    piEN_fit = 0.0
    piEE_fit = 0.0
    red_chi2 = chi2_fit / float(len(lightcurve) - 5)
    # model_telescope = toolbox.fake_telescopes.create_a_fake_telescope()
    # pyLIMA_parameters = pspl.compute_pyLIMA_parameters(fit_tap.fit_results["best_model"])
    # model = pspl.compute_the_microlensing_model(model_telescope, pyLIMA_parameters)
    model_telescope = None
    try:
        to_return = [np.around(t0_fit, 3), np.around(u0_fit, 5), np.around(tE_fit, 3), np.around(piEN_fit, 5),
                     np.around(piEE_fit, 5),
                     np.around(mag_source_fit, 3), np.around(mag_blend_fit, 3), np.around(mag_baseline_fit, 3),
                     fit_tap.fit_results["covariance_matrix"], model_telescope, np.around(chi2_fit, 3), red_chi2]
    except:
        to_return = [np.around(t0_fit, 3), np.around(u0_fit, 5), np.around(tE_fit, 3), np.around(piEN_fit, 5),
                     np.around(piEE_fit, 5),
                     np.around(mag_source_fit, 3), mag_blend_fit, np.around(mag_baseline_fit, 3),
                     fit_tap.fit_results["covariance_matrix"], model_telescope, np.around(chi2_fit, 3), red_chi2]
    return to_return
