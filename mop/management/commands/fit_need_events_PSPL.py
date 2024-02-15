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

from django.db import connection
from django.db.models import prefetch_related_objects

def run_fit(target, red_data, cores=0):
    """
    Function to perform a microlensing model fit to timeseries photometry.

    Parameters:
        target   Target object
        red_data QuerySet of ReducedDatums for the target
        cores integer, optional number of processing cores to use
    """

    logger.info('Fitting event: '+target.name)

    t5 = datetime.datetime.utcnow()
    print('N DB connections run_fit 1: ' + str(len(connection.queries)))

#try:
    t6 = datetime.datetime.utcnow()
    print('N DB connections run_fit 2: ' + str(len(connection.queries)))
    print('Time taken: ' + str(t6 - t5))

    # Retrieve all available ReducedDatum entries for this target.  Note that this may include data
    # other than lightcurve photometry, so the data are then filtered and repackaged for later
    # convenience
    (datasets, ndata) = fittools.repackage_lightcurves(red_data)

    logger.info('FIT: Found '+str(len(datasets))+' datasets and a total of '
                +str(ndata)+' datapoints to model for event '+target.name)

    t7 = datetime.datetime.utcnow()
    print('N DB connections run_fit 3: ' + str(len(connection.queries)))
    print('Time taken: ' + str(t7 - t6))

    if ndata > 10:
        (model_params, model_telescope) = fittools.fit_pspl_omega2(target.ra, target.dec, datasets)
        logger.info('FIT: completed modeling process for '+target.name)

        t8 = datetime.datetime.utcnow()
        print('N DB connections run_fit 3: ' + str(len(connection.queries)))
        print('Time taken: ' + str(t8 - t7))

        # Store model lightcurve
        if model_telescope:
            fittools.store_model_lightcurve(target, model_telescope)
            logger.info('FIT: Stored model lightcurve for event '+target.name)
        else:
            logger.warning('FIT: No valid model fit produced so not model lightcurve for event '+target.name)

        t9 = datetime.datetime.utcnow()
        print('N DB connections run_fit 4: ' + str(len(connection.queries)))
        print('Time taken: ' + str(t9 - t8))

        # Determine whether or not an event is still active based on the
        # current time relative to its t0 and tE
        alive = fittools.check_event_alive(model_params['t0'], model_params['tE'])

        t10 = datetime.datetime.utcnow()
        print('N DB connections run_fit 5: ' + str(len(connection.queries)))
        print('Time taken: ' + str(t10 - t9))

        # Store model parameters
        last_fit = Time(datetime.datetime.utcnow()).jd

        extras = {'Alive':alive, 'Last_fit': last_fit}
        store_keys = ['t0', 't0_error', 'u0', 'u0_error', 'tE', 'tE_error',
                      'piEN', 'piEN_error', 'piEE', 'piEE_error',
                      'Source_magnitude', 'Source_mag_error',
                      'Blend_magnitude', 'Blend_mag_error',
                      'Baseline_magnitude', 'Baseline_mag_error',
                      'Fit_covariance', 'chi2', 'red_chi2',
                      'KS_test', 'AD_test', 'SW_test']
        for key in store_keys:
            extras[key] = model_params[key]
        logger.info('Fitted parameters for '+target.name+': '+repr(extras))

        target.save(extras = extras)
        logger.info('FIT: Stored model parameters for event ' + target.name)

        t11 = datetime.datetime.utcnow()
        print('N DB connections run_fit 6: ' + str(len(connection.queries)))
        print('Time taken: ' + str(t11 - t10))

    else:
        logger.info('Insufficient lightcurve data available to model event '+target.name)

        # Return True because no further processing is required
        return True

    #except:
    #    logger.error('Job failed: '+target.name)
    #    return False

class Command(BaseCommand):
    help = 'Fit events with PSPL and parallax, then ingest fit parameters in the db'

    def add_arguments(self, parser):
        parser.add_argument('--cores', help='Number of workers (CPU cores) to use', default=os.cpu_count(), type=int)
        parser.add_argument('--run-every', help='Run each Fit every N hours', default=4, type=int)

    def handle(self, *args, **options):

        with transaction.atomic():

            t1 = datetime.datetime.utcnow()
            print('N DB connections 1: ' + str(len(connection.queries)))

            # Cutoff date: N hours ago (from the "--run-every=N" hours command line argument)
            cutoff = Time(datetime.datetime.utcnow() - datetime.timedelta(hours=options['run_every'])).jd
            print(cutoff)
            # Find alive microlensing targets for which the fits need to be updated
            # https://docs.djangoproject.com/en/3.0/ref/models/querysets/#select-for-update

            # Prefetch allows us to load data from other tables related to the queryset results, by reversing
            # a foreign key.  TargetExtras have Targets as a foreign key, but Targets don't have TargetExtras, so this
            # query has to search first on the TargetExtras table.

            t = Target.objects.get(name='Gaia23aiy')
            print('Extra fields for ' + t.name + ' include classification='
                  + t.extra_fields['Classification'] + ' and Alive='
                  + repr(t.extra_fields['Alive']))

            ts = TargetExtra.objects.filter(
                key='Classification', value__icontains='Microlensing binary'
                ).filter(
                key='Alive', value=True
                )
            print(ts)

            # ADD BACK CUTOFF CRITERION
            #targetset = Target.objects.select_for_update(skip_locked=True).select_related('extra_field','reduceddatum')
            #targetset.filter(
            #    targetextra__in=TargetExtra.objects.filter(key='Last_fit', value__lte=cutoff)
            #    )
            #targetset.filter(
            #    targetextra__in=TargetExtra.objects.filter(key='Alive', value=True)
            #)
            #targetset.filter(
            #    targetextra__in=TargetExtra.objects.filter(key='Classification', value__icontains='microlensing')
            #)

            #print(targetset)
            print('N DB connections 2: ' + str(len(connection.queries)))
            #print(ts[0].extra_fields)
            print('N DB connections 2: ' + str(len(connection.queries)))
            exit()

            t2 = datetime.datetime.utcnow()
            print('N DB connections 2: ' + str(len(connection.queries)))
            print('Time for targets query: ' + str(t2 - t1))

            # Inform the user how much work is left in the queue
            logger.info('Target(s) remaining in processing queue: '+str(targetset.count()))

            # Loop through all targets in the set
            for element in targetset:

                need_to_fit = True

                try:
                    last_fit = element.extra_fields['Last_fit']
                    red_data = ReducedDatum.objects.filter(target=element).order_by("timestamp")
                    time = [Time(i.timestamp).jd for i in red_data if i.data_type == 'photometry']
                    last_observation = max(time)

                    existing_model = None
                    for dset in red_data:
                        if dset.data_type == 'lc_model':
                            existing_model = dset
                            break

                    if (last_observation<last_fit) & (existing_model.count() != 0) :
                        need_to_fit = False
                except:

                    pass


                t3 = datetime.datetime.utcnow()
                print('N DB connections 3: ' + str(len(connection.queries)))
                print('Time for ReducedDatums query for 1 target: ' + str(t3 - t2))

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

            t4 = datetime.datetime.utcnow()
            print('N DB connections 4: ' + str(len(connection.queries)))
            print('Time for storing Last_fit timestamp: ' + str(t4 - t3))

            if need_to_fit:
                result = run_fit(element, red_data, cores=options['cores'])


if __name__ == '__main__':
    main()
