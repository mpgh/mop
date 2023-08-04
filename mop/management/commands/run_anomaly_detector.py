from django.core.management.base import BaseCommand
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target,TargetExtra,TargetList
from astropy.time import Time
from mop.toolbox import anomaly_detector

import datetime
import json
import numpy as np

class Command(BaseCommand):

    help = 'Assess anomalous status of events'

    def add_arguments(self, parser):

        parser.add_argument('target_name', help='name of the event to fit')

    def handle(self, *args, **options):

        ### Create or load Anomalous list
        try:

            #tap_list = TargetList.objects.filter(name='OMEGAII')[0]
            
            anomalous_list = TargetList.objects.filter(name='Anomalous')[0]

        except:

            anomalous_list = TargetList(name='Anomalous')
            anomalous_list.save()

        if options['target_name'] == 'all':

            list_of_events_alive = Target.objects.filter(targetextra__in=TargetExtra.objects.filter(key='Alive', value=True))
            
        else:

            target, created = Target.objects.get_or_create(name= options['target_name'])
            list_of_events_alive = [target]
        
        KMTNet_fields = TAP.load_KMTNet_fields()

        for event in list_of_events_alive[:]:

            if 'Microlensing' not in event.extra_fields['Classification']:
               pass

            ######
            #elif 'Gaia' in event.name:
            #   pass

            else:
                    try:


                        KS_test = event.extra_fields['KS_test']
                        SW_test = event.extra_fields['SW_test']
                        chi2dof = event.extra_fields['chi2dof']

                        anomaly_status = anomaly_detector.assess_anomaly(KS_test,SW_test,chi2dof)
     
                        #Do we want to stre the flags time serie? Similar to below
                        
##                        data = {'tap': planet_priority,
##                                'tap_error': planet_priority_error
##                                }

##                        rd, created = ReducedDatum.objects.get_or_create(
##                                  timestamp=datetime.datetime.utcnow(),
##                                  value=data,
##                                  source_name='MOP',
##                                  source_location=event.name,
##                                  data_type='TAP_priority',
##                                  target=event)

#                        if created:
#                            rd.save()
#                        extras = {'TAP_priority':np.around(planet_priority,5)}
#                        event.save(extras = extras)
                         
                         if anomaly_status == 'Anomalous':
                         
                            anomalous_list.targets.add(event)
                         
                    except:

                        print('Can not perform Anomaly assesment for this target')
