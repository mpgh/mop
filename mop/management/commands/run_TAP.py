from django.core.management.base import BaseCommand
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target,TargetExtra,TargetList
from astropy.time import Time
from mop.toolbox import TAP
from mop.toolbox import obs_control
from mop.toolbox import omegaII_strategy
import datetime
import json
import numpy as np

class Command(BaseCommand):

    help = 'Sort events that need extra observations'

    def add_arguments(self, parser):

        parser.add_argument('target_name', help='name of the event to fit')

    def handle(self, *args, **options):

        ### Create or load TAP list
        try:

            tap_list = TargetList.objects.filter(name='OMEGAII')[0]

        except:

            tap_list = TargetList(name='OMEGAII')
            tap_list.save()

        if options['target_name'] == 'all':


            list_of_events_alive = Target.objects.filter(targetextra__in=TargetExtra.objects.filter(key='Alive', value=True))
        else:

            target, created = Target.objects.get_or_create(name= options['target_name'])
            list_of_events_alive = [target]
        
        KMTNet_fields = TAP.load_KMTNet_fields()

        for event in list_of_events_alive[:]:

            if 'Microlensing' not in event.extra_fields['Classification']:
               pass

            else:
                    try:

                        time_now = Time(datetime.datetime.now()).jd
                        t0_pspl = event.extra_fields['t0']
                        u0_pspl = event.extra_fields['u0']
                        tE_pspl = event.extra_fields['tE']

                        covariance = np.array(json.loads(event.extra_fields['Fit_covariance']))

                        # Categorize the event based on event timescale
                        category = TAP.categorize_event_timescale(event)

                        # Calculate the priority of this event for different science goals
                        planet_priority = TAP.TAP_planet_priority(time_now,t0_pspl,u0_pspl,tE_pspl)
                        planet_priority_error = TAP.TAP_planet_priority_error(time_now,t0_pspl,u0_pspl,tE_pspl,covariance)
                        print('Planet priority: ',planet_priority, planet_priority_error)

                        # ACTION RAS: Need to calculate this if not already available
                        t_last = event.extra_fields['Latest_data_HJD']
                        print('Last datapoint: ',t_last)

                        long_priority = TAP.TAP_long_event_priority(time_now, t_last, tE_pspl)
                        long_priority_error = TAP.TAP_long_event_priority_error(tE_pspl, covariance)
                        print('Long tE priority: ',long_priority, long_priority_error)

                        # Storing both types of priority as extra_params and also as ReducedDatums so
                        # that we can track the evolution of the priority as a function of time
                        extras = {'TAP_priority':np.around(planet_priority,5),
                                  'TAP_priority_longtE': np.around(long_priority, 5)}
                        event.save(extras = extras)
                        print('Starting reduced datum ingest')

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
                        print('Post first rd create')

                        rd, created = ReducedDatum.objects.get_or_create(
                            timestamp=datetime.datetime.utcnow(),
                            value=data,
                            source_name='MOP',
                            source_location=event.name,
                            data_type='TAP_priority_longtE',
                            target=event)
                        if created:
                            rd.save()
                        print('Got here')
                        # Gather the information required to make a strategy decision
                        # for this target:

                        # Exclude events that are within the High Cadence Zone
                        # event_in_the_Bulge = TAP.event_in_the_Bulge(event.ra, event.dec)
                        event_not_in_OMEGA_II = TAP.event_not_in_OMEGA_II(event.ra, event.dec, KMTNet_fields)
                        print('Event NOT in OMEGA: ', event_not_in_OMEGA_II)

                        # If the event is in the HCZ, set the MOP flag to not observe it
                        if (event_not_in_OMEGA_II):# (event_in_the_Bulge or)  & (event.extra_fields['Baseline_magnitude']>17):

                            extras = {'Observing_mode':'No'}
                            event.save(extras = extras)
                            print('Event not in OMEGA-II')

                        # If the event is flagged as not alive, then it is over, and should also not be observed
                        elif not event.extra_fields['Alive']:
                            extras = {'Observing_mode': 'No'}
                            event.save(extras=extras)
                            print('Event not Alive')

                        # For Alive events outside the HCZ, the strategy depends on whether it is classified as a
                        # stellar/planetary event, or a long-timescale black hole candidate
                        else:
                            # CHECK HERE FOR VISIBILITY
                            print('Would be scheduling obs')
                            mag_now = TAP.TAP_mag_now(event)
                            print('Mag now = ',mag_now)
                            mag_now = 17.5          ### REMOVE THIS TESTING ONLY
                            mag_baseline = 18.5     ### REMOVE THIS TESTING ONLY
                            if mag_now:
                                mag_baseline = event.extra_fields['Baseline_magnitude']
                                print('mag_baseline: ', mag_baseline)
                                observing_mode = TAP.TAP_observing_mode(planet_priority, planet_priority_error,
                                                                    long_priority, long_priority_error,
                                                                    mag_now, mag_baseline)
                            else:
                                observing_mode = None
                            print(event.name, observing_mode)
                            extras = {'Observing_mode': observing_mode}
                            event.save(extras=extras)

                            observe = False
                            if observing_mode in ['priority_stellar_event', 'priority_long_event', 'regular_long_event']:
                                tap_list.targets.add(event)

                                # Get the observational configurations for the event, based on the OMEGA-II strategy:
                                obs_configs = omegaII_strategy.determine_obs_config(event, observing_mode,
                                                                                    mag_now, time_now, t0_pspl, tE_pspl)

                                # Filter this list of hypothetical observations, removing any for which a similar
                                # request has already been submitted and has status 'PENDING'
                                obs_configs = obs_control.filter_duplicated_observations(obs_configs)

                                # Build the corresponding observation requests in LCO format:
                                obs_requests = obs_control.build_lco_imaging_request(obs_configs)

                                # Submit the set of observation requests:
                                if observe:
                                    submit_lco_obs_requests(obs_requests)

                        ### Spectroscopy
                        if (event.extra_fields['Spectras']<1) & (event.extra_fields['Observing_mode'] != 'No'):
                            obs_control.build_and_submit_regular_spectro(event)

                    except:

                        print('Can not perform TAP for this target')
