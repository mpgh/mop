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
import traceback
import datetime
import random
import json
import sys
import os
import logging

logger = logging.getLogger(__name__)

def run_fit(target, cores):
    logger.info('Fitting event: '+target.name)

    #try:
    if 'Gaia' in target.name:

        gaia_mop.update_gaia_errors(target)

    # Add photometry model

    if 'Microlensing' not in target.extra_fields['Classification']:
        alive = False

        extras = {'Alive':alive}
        target.save(extras = extras)

    else:

        # Retrieve all available ReducedDatum entries for this target.  Note that this may include data
        # other than lightcurve photometry, so the data are then filtered and repackaged for later
        # convenience
        red_data = ReducedDatum.objects.filter(target=target).order_by("timestamp")

        (datasets, ndata) = fittools.repackage_lightcurves(red_data)

        logger.info('FIT: Found '+str(len(datasets))+' datasets and a total of '
                    +str(ndata)+' datapoints to model for event '+target.name)

        if ndata > 10:
            (model_params, model_telescope) = fittools.fit_pspl_omega2(target.ra, target.dec, datasets)
            logger.info('FIT: completed modeling process for '+target.name)

            # Store model lightcurve
            if model_telescope:
                fittools.store_model_lightcurve(target, model_telescope)
                logger.info('FIT: Stored model lightcurve for event '+target.name)
            else:
                logger.warning('FIT: No valid model fit produced so not model lightcurve for event '+target.name)

            # Determine whether or not an event is still active based on the
            # current time relative to its t0 and tE
            alive = fittools.check_event_alive(model_params['t0'], model_params['tE'])

            # Store model parameters
            last_fit = Time(datetime.datetime.utcnow()).jd

            extras = {'Alive':alive, 'Last_fit': last_fit}
            store_keys = ['t0', 'u0', 'tE', 'piEN', 'piEE',
                          'Source_magnitude', 'Blend_magnitude', 'Baseline_magnitude',
                          'Fit_covariance', 'chi2', 'red_chi2',
                          'KS_test', 'AD_test', 'SW_test']
            for key in store_keys:
                extras[key] = model_params[key]
            logger.info('Fitted parameters for '+target.name+': '+repr(extras))

            target.save(extras = extras)
            logger.info('FIT: Stored model parameters for event ' + target.name)

        else:
            logger.info('Insufficient lightcurve data available to model event '+target.name)

    #except:
    #    logger.error('Job failed: '+target.name)
    #    return None

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
                logger.info('Target(s) remaining in processing queue: '+str(queryset.count()))

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
                logger.info('Job is finished, no more objects left to process! Goodbye!')
                sys.exit(0)

            # Now we know for sure we have an element to process, and we haven't locked
            # a row (object) in the database. We're free to process this for up to four hours.

            # Check if the fit is really needed, i.e. is there new data since the last fit?


            if need_to_fit:
                result = run_fit(element, cores=options['cores'])


if __name__ == '__main__':
    main()
