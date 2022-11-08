from django.core.management.base import BaseCommand
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target,TargetExtra,TargetList
from astropy.time import Time
from mop.toolbox import TAP
from mop.toolbox import obs_control
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

            tap_list = TargetList.objects.filter(name='TAP')[0]

        except:

            tap_list = TargetList(name='TAP')
            tap_list.save()

        if options['target_name'] == 'all':


            list_of_events_alive = Target.objects.filter(targetextra__in=TargetExtra.objects.filter(key='Alive', value=True))
        else:

            target, created = Target.objects.get_or_create(name= options['target_name'])
            list_of_events_alive = [target]

        for event in list_of_events_alive[:]:

            if 'Microlensing' not in event.extra_fields['Classification']:
               pass

            ######
            elif 'Gaia' in event.name:
               pass

            else:
                    try:

                        time_now = Time(datetime.datetime.now()).jd
                        t0_pspl = event.extra_fields['t0']
                        u0_pspl = event.extra_fields['u0']
                        tE_pspl = event.extra_fields['tE']


                        covariance = np.array(json.loads(event.extra_fields['Fit_covariance']))

                        planet_priority = TAP.TAP_planet_priority(time_now,t0_pspl,u0_pspl,tE_pspl)
                        planet_priority_error = TAP.TAP_planet_priority_error(time_now,t0_pspl,u0_pspl,tE_pspl,covariance)

                        #psi_deriv = TAP.psi_derivatives_squared(time_now,t0_pspl,u0_pspl,tE_pspl)
                        #error = (psi_deriv[2] * covariance[2,2] + psi_deriv[1] * covariance[1,1] + psi_deriv[0] * covariance[0,0])**0.5
                        ### need to create a reducedatum for planet priority


                        data = {'tap': planet_priority,
                                'tap_error': planet_priority_error
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
                        extras = {'TAP_priority':np.around(planet_priority,5)}
                        event.save(extras = extras)




                        ### Regular obs

                        event_in_the_Bulge = TAP.event_in_the_Bulge(event.ra, event.dec)

                        if (event_in_the_Bulge):# & (event.extra_fields['Baseline_magnitude']>17):

                               extras = {'Observing_mode':'No'}
                               event.save(extras = extras)


                        else:

                           if event.extra_fields['Alive'] == True:
                               extras = {'Observing_mode':'Regular'}
                               event.save(extras = extras)
                               obs_control.build_and_submit_regular_phot(event)
                           else:
                               extras = {'Observing_mode':'No'}
                               event.save(extras = extras)

                        # Priority mode

                        mag_now = TAP.TAP_mag_now(event)
                        mag_baseline = event.extra_fields['Baseline_magnitude']
                        new_observing_mode = TAP.TAP_observing_mode(planet_priority,planet_priority_error,mag_now,mag_baseline)

                        if new_observing_mode:
                           tap_list.targets.add(event)


                           extras = {'Observing_mode':new_observing_mode}
                           event.save(extras = extras)
                           print(planet_priority,planet_priority_error)
                           obs_control.build_and_submit_priority_phot(event)

                        ### Spectroscopy
                        if (event.extra_fields['Spectras']<1) & (event.extra_fields['Observing_mode'] != 'No'):
                            obs_control.build_and_submit_regular_spectro(event)

                    except:

                        print('Can not perform TAP for this target')
