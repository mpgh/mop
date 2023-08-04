from django.core.management.base import BaseCommand
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target,TargetExtra
from django.db import transaction
from django.db import models
from astropy.time import Time
from mop.toolbox import fittools
from mop.brokers import gaia as gaia_mop
from django.db.models import Q
import numpy as np
from mop.toolbox import logs
import traceback
import datetime
import random
import json
import sys
import os


def run_fit(target, cores):
    print(f'Working on {target.name}')

    # Start logging process:
    log = logs.start_log()
    log.info('Fitting needed event: '+target.name)

    try:
       if 'Gaia' in target.name:

           gaia_mop.update_gaia_errors(target)

       # Add photometry model

       if 'Microlensing' not in target.extra_fields['Classification']:
           alive = False

           extras = {'Alive':alive}
           target.save(extras = extras)

       else:


           datasets = ReducedDatum.objects.filter(target=target)
           time = [Time(i.timestamp).jd for i in datasets if i.data_type == 'photometry']

           phot = []
           for data in datasets:
               if data.data_type == 'photometry':
                  try:
                       phot.append([data.value['magnitude'],data.value['error'],data.value['filter']])

                  except:
                       # Weights == 1
                       phot.append([data.value['magnitude'],1,data.value['filter']])


           photometry = np.c_[time,phot]

           t0_fit, u0_fit, tE_fit, piEN_fit, piEE_fit, mag_source_fit, mag_blend_fit, mag_baseline_fit, cov, model, chi2_fit, red_chi2, sw_test, ad_test, ks_test = fittools.fit_pspl_omega2(target.ra, target.dec, photometry)

           # Add photometry model

           model_time = datetime.datetime.strptime('2018-06-29 08:15:27.243860', '%Y-%m-%d %H:%M:%S.%f')
           data = {'lc_model_time': model.lightcurve_magnitude[:,0].tolist(),
           'lc_model_magnitude': model.lightcurve_magnitude[:,1].tolist()
                   }
           existing_model =   ReducedDatum.objects.filter(source_name='MOP',data_type='lc_model',
                                                         timestamp=model_time,source_location=target.name)


           if existing_model.count() == 0:
               rd, created = ReducedDatum.objects.get_or_create(
                                                                   timestamp=model_time,
                                                                   value=data,
                                                                   source_name='MOP',
                                                                   source_location=target.name,
                                                                   data_type='lc_model',
                                                                   target=target)

               rd.save()

           else:
               rd, created = ReducedDatum.objects.update_or_create(
                                                                   timestamp=existing_model[0].timestamp,
                                                                   value=existing_model[0].value,
                                                                   source_name='MOP',
                                                                   source_location=target.name,
                                                                   data_type='lc_model',
                                                                   target=target,
                                                                   defaults={'value':data})

               rd.save()


           time_now = Time(datetime.datetime.now()).jd
           how_many_tE = (time_now-t0_fit)/tE_fit


           if how_many_tE>2:

              alive = False

           else:

              alive = True

           last_fit = Time(datetime.datetime.utcnow()).jd


           extras = {'Alive':alive, 't0':t0_fit,'u0':u0_fit,'tE':tE_fit,
                         'piEN':piEN_fit,'piEE':piEE_fit,
                         'Source_magnitude':mag_source_fit,
                         'Blend_magnitude':mag_blend_fit,
                         'Baseline_magnitude':mag_baseline_fit,
                         'Fit_covariance':json.dumps(cov.tolist()),
                         'chi2':chi2_fit,
                         'red_chi2': red_chi2,
                         'Last_Fit':last_fit}
           log.info('Fitted parameters for '+target.name+': '+repr(extras))

           target.save(extras = extras)

       logs.stop_log(log)

    except:
        print(f'Job failed: {target.name}')
        logs.info(f'ERROR: Job failed: {target.name}')
        traceback.print_exc()
        logs.stop_log(log)
        return None

class Command(BaseCommand):
    help = 'Fit events with PSPL and parallax, then ingest fit parameters in the db'

    def add_arguments(self, parser):
        parser.add_argument('--cores', help='Number of workers (CPU cores) to use', default=os.cpu_count(), type=int)
        parser.add_argument('--run-every', help='Run each Fit every N hours', default=4, type=int)

    def handle(self, *args, **options):
        # The TOM Toolkit project does not automatically create key/value pairs
        # in the "Extra Fields" during a database migration. We use this silly
        # method to automatically add this field to any Target objects which were
        # created before this field was added to the database.
        # Adding Last_fit if dos not exist
        list_of_targets = Target.objects.filter()
        for target in list_of_targets:
            try:
                last_fit = target.extra_fields['Last_fit']
            except:
                last_fit = 2446756.50000
                extras = {'Last_fit':last_fit}
                target.save(extras = extras)

        # Run until all objects which need processing have been processed
        while True:
            # One instance of our database model to process (if found)
            element = None

            # Select the first available un-claimed object for processing. We indicate
            # ownership of the job by advancing the timestamp to the current time. This
            # ensures that we don't have two workers running the same job. A beneficial
            # side effect of this implementation is that a job which crashes isn't retried
            # for another N hours, which limits the potential impact.
            #
            # The only time this system breaks down is if a single processing fit takes
            # more than N hours. We'll instruct Kubernetes that no data processing Pod
            # should run for that long. That'll protect us against that overrun scenario.
            #
            # The whole thing is wrapped in a database transaction to protect against
            # collisions by two workers. Very unlikely, but we're good software engineers
            # and will protect against that.
            with transaction.atomic():

                # Cutoff date: N hours ago (from the "--run-every=N" hours command line argument)
                cutoff = Time(datetime.datetime.utcnow() - datetime.timedelta(hours=options['run_every'])).jd

                # Find any objects which need to be run
                # https://docs.djangoproject.com/en/3.0/ref/models/querysets/#select-for-update
                queryset = Target.objects.select_for_update(skip_locked=True)
                queryset = queryset.filter(targetextra__in=TargetExtra.objects.filter(key='Last_fit', value__lte=cutoff))
                queryset = queryset.filter(targetextra__in=TargetExtra.objects.filter(key='Alive', value=True))

                # Inform the user how much work is left in the queue
                print('Target(s) remaining in processing queue:', queryset.count())

                # Retrieve the first element which meets the condition
                element = queryset.first()

                need_to_fit = True

                try:
                    last_fit = element.extra_fields['Last_fit']
                    datasets = ReducedDatum.objects.filter(target=element)
                    time = [Time(i.timestamp).jd for i in datasets if i.data_type == 'photometry']
                    last_observation = max(time)
                    existing_model = ReducedDatum.objects.filter(source_name='MOP',data_type='lc_model',source_location=element.name)

                    if (last_observation<last_fit) & (existing_model.count() != 0) :
                        need_to_fit = False
                except:

                    pass

                # Element was found. Claim the element for this worker (mark the fit as in
                # the "RUNNING" state) by setting the Last_fit timestamp. This method has
                # the beneficial side effect such that if a fit crashes, it won't be re-run
                # (retried) for another N hours. This limits the impact of broken code on the cluster.
                if element is not None:
                    last_fit = Time(datetime.datetime.utcnow()).jd
                    extras = {'Last_fit':last_fit}
                    element.save(extras = extras)


            # If there are no more objects left to process, then the job is finished.
            # Inform Kubernetes of this fact by exiting successfully.
            if element is None:
                print('Job is finished, no more objects left to process! Goodbye!')
                sys.exit(0)

            # Now we know for sure we have an element to process, and we haven't locked
            # a row (object) in the database. We're free to process this for up to four hours.

            # Check if the fit is really needed, i.e. is there new data since the last fit?


            if need_to_fit:
                result = run_fit(element, cores=options['cores'])


if __name__ == '__main__':
    main()
