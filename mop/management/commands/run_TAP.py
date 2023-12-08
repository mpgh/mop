from django.core.management.base import BaseCommand
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target,TargetExtra,TargetList
from astropy.time import Time, TimeDelta
from mop.toolbox import TAP
from mop.toolbox import TAP_priority
from mop.toolbox import obs_control
from mop.toolbox import omegaII_strategy
from mop.toolbox import interferometry_prediction
import datetime
import json
import numpy as np
import logging

logger = logging.getLogger(__name__)

class Command(BaseCommand):

    help = 'Sort events that need extra observations'

    def add_arguments(self, parser):

        parser.add_argument('target_name', help='name of the event to fit')
        parser.add_argument('observe', help='whether or not to observe (live_obs or none)')

    def handle(self, *args, **options):
        verbose = False

        logger.info("runTAP: Started with options "+repr(options))

        ### Create or load TAP list
        try:

            tap_list = TargetList.objects.filter(name='OMEGAII')[0]

        except:

            tap_list = TargetList(name='OMEGAII')
            tap_list.save()

        if options['target_name'] == 'all':
            list_of_events_alive = Target.objects.filter(targetextra__in=TargetExtra.objects.filter(key='Alive', value=True))
            logger.info('runTAP: Loaded all '+str(len(list_of_events_alive))+' events')

        else:
            qs = Target.objects.filter(name= options['target_name'])
            if qs.count() == 0:
                logger.error('runTAP: Cannot find requested target ' + options['target_name'])
                list_of_events_alive = []
            else:
                target = qs[0]
                list_of_events_alive = [target]
                logger.info('runTAP: Loaded event ' + target.name)
        nalive = str(len(list_of_events_alive[:]))

        KMTNet_fields = TAP.load_KMTNet_fields()

        ## Get list of targets for which there are currently-pending observations already in the LCO Portal.
        response = obs_control.fetch_pending_lco_requestgroups()
        pending_obs = obs_control.parse_lco_requestgroups(response)
        logger.info('runTAP: identified pending observations for ' + str(len(pending_obs)) \
                    + ' targets: ' + repr(pending_obs))

        for k,event in enumerate(list_of_events_alive[:]):
            logger.info('runTAP: analyzing event ' + event.name + ', ' + str(k) + ' out of ' + nalive)

            if 'Microlensing' not in event.extra_fields['Classification']:
               pass

            else:
                try:

                    # Gather the necessary information from the model fit to this event.  Sanity check: if this information
                    # isn't available, skip the event
                    try:
                        time_now = Time(datetime.datetime.now()).jd
                        t0_pspl = event.extra_fields['t0']
                        u0_pspl = event.extra_fields['u0']
                        tE_pspl = event.extra_fields['tE']
                        t0_pspl_error = event.extra_fields['t0_error']
                        tE_pspl_error = event.extra_fields['tE_error']
                        red_chi2 = event.extra_fields['red_chi2']

                        covariance = load_covar_matrix(event.extra_fields['Fit_covariance'])

                        sane = TAP.sanity_check_model_parameters(t0_pspl, t0_pspl_error, u0_pspl,
                                                                 tE_pspl, tE_pspl_error, red_chi2, covariance)

                    except KeyError:
                        logger.warning('runTAP: Insufficent model parameters available for ' + event.name + ', skipping')
                        sane = False

                    if sane:
                        # Categorize the event based on event timescale
                        category = TAP.categorize_event_timescale(event)

                        # Calculate the priority of this event for different science goals
                        planet_priority = TAP_priority.TAP_planet_priority(time_now,t0_pspl,u0_pspl,tE_pspl)
                        planet_priority_error = TAP_priority.TAP_planet_priority_error(time_now,t0_pspl,u0_pspl,tE_pspl,covariance)
                        logger.info('runTAP: Planet priority: ' + str(planet_priority) + ', ' + str(planet_priority_error))

                        # ACTION RAS: Need to calculate this if not already available.
                        # If not available, assume a default 30d period to encourage observations
                        (t_last, t_last_date) = TAP.TAP_time_last_datapoint(event)
                        if not t_last:
                            t_last = Time.now(jd) - TimeDelta(days=30.0)
                        logger.info('runTAP: Last datapoint: ' + str(t_last))

                        mag_now = TAP.TAP_mag_now(event)
                        logger.info('runTAP: Mag now = ' + str(mag_now))

                        long_priority = TAP_priority.TAP_long_event_priority(time_now, t_last, tE_pspl)
                        long_priority_error = TAP_priority.TAP_long_event_priority_error(tE_pspl, covariance)
                        logger.info('runTAP: Long tE priority: ' + str(long_priority) + ' ' + str(long_priority_error))

                        # Storing both types of priority as extra_params and also as ReducedDatums so
                        # that we can track the evolution of the priority as a function of time
                        extras = {'TAP_priority':np.around(planet_priority,5),
                                  'TAP_priority_error': np.around(planet_priority_error, 5),
                                  'TAP_priority_longtE': np.around(long_priority, 5),
                                  'TAP_priority_longtE_error': np.around(long_priority_error, 5)}
                        event.save(extras = extras)

                        data = {'tap_planet': planet_priority,
                                'tap_planet_error': planet_priority_error,
                                'tap_long': long_priority,
                                'tap_long_error': long_priority_error
                                }

                        rd, created = ReducedDatum.objects.get_or_create(
                                  timestamp=datetime.datetime.utcnow(),
                                  value=data,
                                  source_name='MOP',
                                  source_location=event.name,
                                  data_type='TAP_priority',
                                  target=event)
                        if created:
                            rd.save()

                        rd, created = ReducedDatum.objects.get_or_create(
                            timestamp=datetime.datetime.utcnow(),
                            value=data,
                            source_name='MOP',
                            source_location=event.name,
                            data_type='TAP_priority_longtE',
                            target=event)
                        if created:
                            rd.save()

                        # Gather the information required to make a strategy decision
                        # for this target:

                        # Exclude events that are within the High Cadence Zone
                        # event_in_the_Bulge = TAP.event_in_the_Bulge(event.ra, event.dec)
                        event_in_HCZ = TAP.event_in_HCZ(event.ra, event.dec, KMTNet_fields)
                        if event_in_HCZ:
                            sky_location = 'In HCZ'
                        else:
                            sky_location = 'Outside HCZ'
                        logger.info('runTAP: Event in HCZ: ' + str(event_in_HCZ))
                        logger.info('runTAP: Event alive? ' + repr(event.extra_fields['Alive']))

                        # If the event is in the HCZ, set the MOP flag to not observe it
                        if (event_in_HCZ):# (event_in_the_Bulge or)  & (event.extra_fields['Baseline_magnitude']>17):

                            extras = {'Observing_mode':'No', 'Sky_location': sky_location}
                            event.save(extras = extras)
                            logger.info('runTAP: Event in HCZ')

                        # If the event is flagged as not alive, then it is over, and should also not be observed
                        elif not event.extra_fields['Alive']:
                            extras = {'Observing_mode': 'No', 'Sky_location': sky_location}
                            event.save(extras=extras)
                            logger.info('runTAP: Event not Alive')

                        # For Alive events outside the HCZ, the strategy depends on whether it is classified as a
                        # stellar/planetary event, or a long-timescale black hole candidate
                        else:
                            logger.info('runTAP: Event should be observed')
                            # Check target for visibility
                            visible = obs_control.check_visibility(event, Time.now().decimalyear, verbose=False)
                            logger.info('runTAP: Event visible?  ' + repr(visible))

                            if visible:
                                if mag_now:
                                    mag_baseline = event.extra_fields['Baseline_magnitude']
                                    logger.info('runTAP: mag_baseline: ' + str(mag_baseline))
                                    observing_mode = TAP.TAP_observing_mode(planet_priority, planet_priority_error,
                                                                        long_priority, long_priority_error,
                                                                        tE_pspl, tE_pspl_error, mag_now,
                                                                        mag_baseline, red_chi2)

                                else:
                                    observing_mode = None
                                logger.info('runTAP: Observing mode: ' + event.name + ' ' + str(observing_mode))
                                extras = {'Observing_mode': observing_mode, 'Sky_location': sky_location}
                                event.save(extras=extras)

                                if observing_mode in ['priority_stellar_event', 'priority_long_event', 'regular_long_event']:
                                    tap_list.targets.add(event)

                                    # Get the observational configurations for the event, based on the OMEGA-II strategy:
                                    obs_configs = omegaII_strategy.determine_obs_config(event, observing_mode,
                                                                                        mag_now, time_now, t0_pspl, tE_pspl)
                                    logger.info('runTAP: Determined observation configurations: ' + repr(obs_configs))

                                    # Filter this list of hypothetical observations, removing any for which a similar
                                    # request has already been submitted and has status 'PENDING'
                                    obs_configs = obs_control.filter_duplicated_observations(obs_configs, pending_obs)
                                    logger.info('runTAP: Filtered out duplicates: ' + repr(obs_configs))

                                    # Build the corresponding observation requests in LCO format:
                                    obs_requests = obs_control.build_lco_imaging_request(obs_configs)
                                    logger.info('runTAP: Build observation requests: ' + repr(obs_requests))

                                    # Submit the set of observation requests:
                                    # Currently observations are restricted to OGLE events only until the Gaia classifier
                                    # is updated
                                    if 'live_obs' in options['observe'] and ('OGLE' in event.name or 'Gaia' in event.name):
                                        obs_control.submit_lco_obs_request(obs_requests, event)
                                        logger.info('runTAP: SUBMITTING OBSERVATIONS')
                                    else:
                                        logger.warning('runTAP: WARNING: OBSERVATIONS SWITCHED OFF')

                            else:
                                logger.info('runTAP: Target ' + event.name + ' not currently visible')
                                observing_mode = None
                                extras = {'Observing_mode': observing_mode, 'Sky_location': sky_location}
                                event.save(extras=extras)

                        ### Spectroscopy
                        observe_spectro = False
                        if observe_spectro:
                            if (event.extra_fields['Spectras']<1) & (event.extra_fields['Observing_mode'] != 'No'):
                                obs_control.build_and_submit_regular_spectro(event)
                                logger.info('runTAP: Submitted spectroscopic observations for ' + event.name)

                        ### Inteferometry
                        interferometry_prediction.evaluate_target_for_interferometry(event)
                        logger.info('runTAP: Evaluated ' + event.name + ' for interferometry')

                except:
                    logger.warning('runTAP: Cannot perform TAP for target ' + event.name)
        logger.info('runTAP: Completed run')

def load_covar_matrix(raw_covar_data):

    payload = str(raw_covar_data).replace('[[','').replace(']]','').replace('\n','').lstrip()
    array_list = payload.split('] [')

    # Check for older covar matrix format
    if len(array_list) == 1:
        array_list = array_list[0].split('], [')

    data = []
    for entry in array_list:
        try:
            data.append([float(x) for x in entry.split()])
        except ValueError:
            data.append([float(x) for x in entry.split(',')])

    return np.array(data)

