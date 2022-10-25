import numpy as np
from pyLIMA import event
from pyLIMA import telescopes
from pyLIMA import microlmodels
from pyLIMA import microloutputs
from pyLIMA import microltoolbox

def chi2(params,fit):

     chi2 = np.sum(fit.residuals_LM(params)**2)
     return chi2

def flux_to_mag(flux):

        ZP_pyLIMA = 27.4
        magnitude = ZP_pyLIMA -2.5*np.log10(flux)
        return magnitude

def fit_PSPL(photometry, emag_limit = None, cores = None):

       current_event = event.Event()
       current_event.name = 'MOP_to_fit'
       filters = np.unique(photometry[:,-1])

       for ind,filt in enumerate(filters):

           if emag_limit:

               mask = (photometry[:,-1] == filt) & (np.abs(photometry[:,-2].astype(float))<emag_limit)

           else:

               mask = (photometry[:,-1] == filt)
           lightcurve = photometry[mask,:-1].astype(float)

           telescope = telescopes.Telescope(name='Tel_'+str(ind), camera_filter=filt,
                                            light_curve_magnitude=lightcurve,
                                            light_curve_magnitude_dictionnary={'time': 0, 'mag': 1, 'err_mag': 2},
                                            clean_the_lightcurve='Yes')

           if len(lightcurve)>5:
               current_event.telescopes.append(telescope)

       Model = microlmodels.create_model('PSPL', current_event, parallax=['None', 0])
       Model.parameters_boundaries[0] = [Model.parameters_boundaries[0][0],Model.parameters_boundaries[0][-1]+500]
       Model.parameters_boundaries[1] = [0,2]

       if cores != 0:

           import multiprocessing
           with multiprocessing.Pool(processes=cores) as pool:
                current_event.fit(Model, 'DE',DE_population_size=10,flux_estimation_MCMC = 'polyfit',computational_pool = pool)

       else:

           current_event.fit(Model, 'DE',DE_population_size=10,flux_estimation_MCMC = 'polyfit')




       t0_fit =  current_event.fits[-1].fit_results[0]
       u0_fit =  current_event.fits[-1].fit_results[1]
       tE_fit =  current_event.fits[-1].fit_results[2]
       chi2_fit = current_event.fits[-1].fit_results[-1]

       mag_source_fit = flux_to_mag( current_event.fits[-1].fit_results[3])
       mag_blend_fit = flux_to_mag( current_event.fits[-1].fit_results[3]*current_event.fits[-1].fit_results[4])
       mag_baseline_fit = flux_to_mag( current_event.fits[-1].fit_results[3]*(1+current_event.fits[-1].fit_results[4]))

       if np.isnan(mag_blend_fit):
           mag_blend_fit = "null"




       return [t0_fit,u0_fit,tE_fit,mag_source_fit,mag_blend_fit,mag_baseline_fit,chi2_fit]

def fit_PSPL_parallax(ra,dec,photometry, emag_limit = None, cores = None):

       # Filter orders
       filters_order = ['I','ip','i_ZTF','r_ZTF','R','g_ZTF','gp','G']


       filters = np.unique(photometry[:,-1])
       order = []
       for fil in filters_order:

           mask = np.where(filters==fil)[0]
           if len(mask)!= 0:
                order += mask.tolist()

       for fil in filters:

           if fil not in filters_order:
                mask = np.where(filters==fil)[0]
                if len(mask)!= 0:
                    order += mask.tolist()
       filters = filters[order]



       t0_fit,u0_fit,tE_fit,mag_source_fit,mag_blend_fit,mag_baseline_fit,chi2_fit = fit_PSPL(photometry, emag_limit = None, cores = cores)

       current_event = event.Event()
       current_event.name = 'MOP_to_fit'

       current_event.ra = ra
       current_event.dec = dec

       for ind,filt in enumerate(filters):

           if emag_limit:

               mask = (photometry[:,-1] == filt) & (np.abs(photometry[:,-2].astype(float))<emag_limit)

           else:

               mask = (photometry[:,-1] == filt)
           lightcurve = photometry[mask,:-1].astype(float)

           telescope = telescopes.Telescope(name='Tel_'+str(ind), camera_filter=filt,
                                            light_curve_magnitude=lightcurve,
                                            light_curve_magnitude_dictionnary={'time': 0, 'mag': 1, 'err_mag': 2},
                                            clean_the_lightcurve='Yes')

           if len(lightcurve)>5:
               current_event.telescopes.append(telescope)


       t0_par = t0_fit

       Model_parallax = microlmodels.create_model('PSPL', current_event, parallax=['Full', t0_par])
       Model_parallax.parameters_boundaries[0] = [t0_fit-100,t0_fit+100]

       Model_parallax.parameters_boundaries[1] = [-2,2]
       Model_parallax.parameters_boundaries[2] = [0.1,500]
       #Model_parallax.parameters_boundaries[3] = [-1,1]
       #Model_parallax.parameters_boundaries[4] = [-1,1]
       Model_parallax.parameters_guess = [ t0_fit,u0_fit,tE_fit,0,0]

       if cores !=0:

           import multiprocessing
           with multiprocessing.Pool(processes=cores) as pool:
                current_event.fit(Model_parallax, 'DE',DE_population_size=10,flux_estimation_MCMC = 'polyfit',computational_pool = pool)

       else:

           current_event.fit(Model_parallax, 'DE',DE_population_size=10,flux_estimation_MCMC = 'polyfit')

       #if (chi2_fit-current_event.fits[-1].fit_results[-1])/current_event.fits[-1].fit_results[-1]<0.1:

            #return [t0_fit,u0_fit,tE_fit,None,None,mag_source_fit,mag_blend_fit,mag_baseline_fit]

       t0_fit = current_event.fits[-1].fit_results[0]
       u0_fit =  current_event.fits[-1].fit_results[1]
       tE_fit =  current_event.fits[-1].fit_results[2]
       piEN_fit =  current_event.fits[-1].fit_results[3]
       piEE_fit =  current_event.fits[-1].fit_results[4]
       chi2_fit = current_event.fits[-1].fit_results[-1]

       mag_source_fit = flux_to_mag( current_event.fits[-1].fit_results[5])
       mag_blend_fit = flux_to_mag( current_event.fits[-1].fit_results[5]*current_event.fits[-1].fit_results[6])
       mag_baseline_fit = flux_to_mag( current_event.fits[-1].fit_results[5]*(1+current_event.fits[-1].fit_results[6]))

       if np.isnan(mag_blend_fit):
           mag_blend_fit = "null"


       microloutputs.create_the_fake_telescopes(current_event.fits[0],current_event.fits[0].fit_results[:-1])
       model_telescope = current_event.fits[0].event.fake_telescopes[0]
       pyLIMA_parameters = Model_parallax.compute_pyLIMA_parameters(current_event.fits[0].fit_results[:-1])
       flux_model = current_event.fits[0].model.compute_the_microlensing_model(model_telescope, pyLIMA_parameters)[0]
       magnitude = microltoolbox.flux_to_magnitude(flux_model)
       model_telescope.lightcurve_magnitude[:,1] = magnitude

       mask = ~np.isnan(magnitude)
       model_telescope.lightcurve_magnitude = model_telescope.lightcurve_magnitude[mask]


       try:
            to_return = [np.around(t0_fit,3),np.around(u0_fit,5),np.around(tE_fit,3),np.around(piEN_fit,5),np.around(piEE_fit,5),
               np.around(mag_source_fit,3),np.around(mag_blend_fit,3),np.around(mag_baseline_fit,3),
               current_event.fits[-1].fit_covariance,model_telescope,np.around(chi2_fit,3)]
       except:
            to_return = [np.around(t0_fit,3),np.around(u0_fit,5),np.around(tE_fit,3),np.around(piEN_fit,5),np.around(piEE_fit,5),
               np.around(mag_source_fit,3),mag_blend_fit,np.around(mag_baseline_fit,3),
               current_event.fits[-1].fit_covariance,model_telescope,np.around(chi2_fit,3)]
       return to_return
