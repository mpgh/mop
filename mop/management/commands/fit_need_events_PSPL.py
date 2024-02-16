from django.core.management.base import BaseCommand
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target,TargetExtra
from django.db import transaction
from astropy.time import Time
from mop.toolbox import fittools
from mop.toolbox.mop_classes import MicrolensingEvent
import datetime
import os
import logging

logger = logging.getLogger(__name__)

from django.db import connection

def run_fit(mulens, cores=0, verbose=False):
    """
    Function to perform a microlensing model fit to timeseries photometry.

    Parameters:
        target   Target object
        red_data QuerySet of ReducedDatums for the target
        cores integer, optional number of processing cores to use
    """

    logger.info('Fitting event: '+mulens.name)

    t5 = datetime.datetime.utcnow()
    if verbose: logger.info('N DB connections run_fit 1: ' + str(len(connection.queries)))

    try:
        t6 = datetime.datetime.utcnow()
        if verbose: logger.info('N DB connections run_fit 2: ' + str(len(connection.queries)))
        if verbose: logger.info('Time taken: ' + str(t6 - t5))

        # Retrieve all available ReducedDatum entries for this target.  Note that this may include data
        # other than lightcurve photometry, so the data are then filtered and repackaged for later
        # convenience
        logger.info('FIT: Found '+str(len(mulens.datasets))+' datasets and a total of '
                    +str(mulens.ndata)+' datapoints to model for event '+mulens.name)

        t7 = datetime.datetime.utcnow()
        if verbose: logger.info('N DB connections run_fit 3: ' + str(len(connection.queries)))
        if verbose: logger.info('Time taken: ' + str(t7 - t6))

        if mulens.ndata > 10:
            (model_params, model_telescope) = fittools.fit_pspl_omega2(mulens.ra, mulens.dec, mulens.datasets)
            logger.info('FIT: completed modeling process for '+mulens.name)

            t8 = datetime.datetime.utcnow()
            if verbose: logger.info('N DB connections run_fit 3: ' + str(len(connection.queries)))
            if verbose: logger.info('Time taken: ' + str(t8 - t7))

            # Store model lightcurve
            if model_telescope:
                mulens.store_model_lightcurve(model_telescope)
                logger.info('FIT: Stored model lightcurve for event '+mulens.name)
            else:
                logger.warning('FIT: No valid model fit produced so not model lightcurve for event '+mulens.name)

            t9 = datetime.datetime.utcnow()
            if verbose: logger.info('N DB connections run_fit 4: ' + str(len(connection.queries)))
            if verbose: logger.info('Time taken: ' + str(t9 - t8))

            # Determine whether or not an event is still active based on the
            # current time relative to its t0 and tE
            alive = fittools.check_event_alive(model_params['t0'], model_params['tE'])

            t10 = datetime.datetime.utcnow()
            if verbose: logger.info('N DB connections run_fit 5: ' + str(len(connection.queries)))
            if verbose: logger.info('Time taken: ' + str(t10 - t9))

            # Store model parameters
            model_params['Last_fit'] = Time(datetime.datetime.utcnow()).jd
            model_params['Alive'] = alive
            mulens.store_model_parameters(model_params)
            logger.info('FIT: Stored model parameters for event ' + mulens.name)

            t11 = datetime.datetime.utcnow()
            if verbose: logger.info('N DB connections run_fit 6: ' + str(len(connection.queries)))
            if verbose: logger.info('Time taken: ' + str(t11 - t10))

        else:
            logger.info('Insufficient lightcurve data available to model event '+mulens.name)

            # Return True because no further processing is required
            return True

    except:
        logger.error('Job failed: '+target.name)
        return False

class Command(BaseCommand):
    help = 'Fit events with PSPL and parallax, then ingest fit parameters in the db'

    def add_arguments(self, parser):
        parser.add_argument('--cores', help='Number of workers (CPU cores) to use', default=os.cpu_count(), type=int)
        parser.add_argument('--run-every', help='Run each Fit every N hours', default=4, type=int)

    def handle(self, *args, **options):

        with transaction.atomic():

            t1 = datetime.datetime.utcnow()
            logger.info('FIT_NEED_EVENTS: Starting with N DB connections: ' + str(len(connection.queries)))

            # Cutoff date: N hours ago (from the "--run-every=N" hours command line argument)
            cutoff = Time(datetime.datetime.utcnow() - datetime.timedelta(hours=options['run_every'])).jd

            # Find alive microlensing targets for which the fits need to be updated
            # We really want the intersection of three querysets: those targets that have separate
            # TargetExtras meeting all criteria.  Using prefetch_related targets here because
            # it is massively more efficient in time and DB connections, together with select_for_update
            # to skip any Target being accessed by another process.
            ts1 = TargetExtra.objects.prefetch_related('target').select_for_update(skip_locked=True).filter(
                key='Classification', value__icontains='Microlensing'
                )
            ts2 = TargetExtra.objects.prefetch_related('target').select_for_update(skip_locked=True).filter(
                key='Alive', value__icontains=True
                )
            ts3 = TargetExtra.objects.prefetch_related('target').select_for_update(skip_locked=True).filter(
                key='Last_fit', value__lte=cutoff
            )
            ts4 = TargetExtra.objects.prefetch_related('target').select_for_update(skip_locked=True).filter(
                key='Latest_data_HJD', value__gt=cutoff
            )
            logger.info('FIT_NEED_EVENTS: Initial queries selected '
                        + str(ts1.count()) + ' events classified as microlensing, '
                        + str(ts2.count()) + ' events currently Alive, '
                        + str(ts3.count()) + ' events last modeled before ' + repr(cutoff)
                        + ', and with data added since then')

            # This doesn't directly produce a queryset of targets, instead it returns a queryset of target IDs.
            # So we have to extract the corresponding targets:
            targets1 = [x.target for x in ts1]
            targets2 = [x.target for x in ts2]
            targets3 = [x.target for x in ts3]
            targets4 = [x.target for x in ts4]
            target_list = list(set(targets1).intersection(
                set(targets2)
            ).intersection(
                set(targets3)
            ).intersection(
                set(targets4))
            )

            logger.info('FIT_NEED_EVENTS: Initial target list has ' + str(len(target_list)) + ' entries')

            target_extras = TargetExtra.objects.filter(
                target__in=target_list
            )
            datums = ReducedDatum.objects.filter(target__in=target_list).order_by("timestamp")

            t2 = datetime.datetime.utcnow()
            logger.info('FIT_NEED_EVENTS: Retrieved associated data for ' + str(len(target_list)) + ' Targets')
            logger.info('FIT_NEED_EVENTS: N DB connections 2: ' + str(len(connection.queries)))
            logger.info('FIT_NEED_EVENTS: Time taken: ' + str(t2 - t1))

            logger.info('FIT_NEED_EVENTS: Reviewing target list to identify those that need remodeling')
            target_data = {}
            for t in target_list:
                mulens = MicrolensingEvent(t)
                mulens.set_extra_params(target_extras.filter(target=t))
                mulens.set_reduced_data(datums.filter(target=t))
                mulens.check_need_to_fit()

                if mulens.need_to_fit:
                    target_data[t] = mulens

            t3 = datetime.datetime.utcnow()
            logger.info('FIT_NEED_EVENTS: Collated data for ' + str(len(target_data)) + ' targets in ' + str(t3 - t2))
            logger.info('FIT_NEED_EVENTS: N DB connections 3: ' + str(len(connection.queries)))

            # Loop through all targets in the set
            for i, (t, mulens) in enumerate(target_data.items()):
                logger.info('FIT_NEED_EVENTS: modeling target ' + mulens.name + ', '
                            + str(i) + ' out of ' + str(len(target_data)))

                t4 = datetime.datetime.utcnow()

                result = run_fit(mulens, cores=options['cores'])

                t5 = datetime.datetime.utcnow()
                logger.info('FIT_NEED_EVENTS: Completed modeling of ' + mulens.name + ' in ' + str(t5 - t4))
                logger.info('FIT_NEED_EVENTS: N DB connections 4: ' + str(len(connection.queries)))

            t6 = datetime.datetime.utcnow()
            logger.info('FIT_NEED_EVENTS: Finished modeling set of targets in ' + str(t6 - t1))
            logger.info('FIT_NEED_EVENTS: N DB connections 5: ' + str(len(connection.queries)))

if __name__ == '__main__':
    main()
