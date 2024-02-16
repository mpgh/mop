from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target,TargetExtra,TargetList
from mop.toolbox.mop_classes import MicrolensingEvent
import logging
import datetime
from django.db import connection

logger = logging.getLogger(__name__)

def get_alive_events(option):
    """
    Function to retrieve Targets that are classified as microlensing events that are currently ongoing.

    Parameters:
        option  str   Can be 'all', in which case the whole database is searched for targets or
                        the name of a specific event.
    """

    if option == 'all':
        t1 = datetime.datetime.utcnow()

       # Search TargetExtras to identify alive microlensing events
        ts1 = TargetExtra.objects.prefetch_related('target').select_for_update(skip_locked=True).filter(
            key='Classification', value__icontains='Microlensing'
        )
        ts2 = TargetExtra.objects.prefetch_related('target').select_for_update(skip_locked=True).filter(
            key='Alive', value__icontains=True
        )
        logger.info('queryTools: Initial queries selected '
                    + str(ts1.count()) + ' events classified as microlensing, '
                    + str(ts2.count()) + ' events currently Alive')

        # This doesn't directly produce a queryset of targets, instead it returns a queryset of target IDs.
        # So we have to extract the corresponding targets:
        targets1 = [x.target for x in ts1]
        targets2 = [x.target for x in ts2]
        target_list = list(set(targets1).intersection(
            set(targets2)
        ))

        logger.info('queryTools: Initial target list has ' + str(len(target_list)) + ' entries')

        # Now gather any TargetExtra and ReducedDatums associated with these targets.
        # This is managed as a dictionary of MicrolensingEvent objects
        target_data = fetch_data_for_targetset(target_list, check_need_to_fit=False)

        t3 = datetime.datetime.utcnow()
        logger.info('queryTools: Collated data for ' + str(len(target_data)) + ' targets in ' + str(t3 - t2))
        logger.info('queryTools: N DB connections: ' + str(len(connection.queries)))

    # Search for a specific target
    else:
        qs = Target.objects.filter(name=option)

        if qs.count() == 1:
            target_list = [qs[0]]
            target_data = fetch_data_for_targetset(target_list, check_need_to_fit=False)

        elif qs.count() == 0:
            logger.error('runTAP: Cannot find requested target ' + option)
            target_data = {}

        else:
            target_data = {}
            logger.error('runTAP: Found multiple events by the name ' + target.name)

    return target_data

def fetch_data_for_targetset(target_list, check_need_to_fit=True):
    """
    Function to retrieve all TargetExtra and ReducedDatums associated with a set of targets
    """
    t1 = datetime.datetime.utcnow()

    # Perform the search for associated data
    target_extras = TargetExtra.objects.filter(
        target__in=target_list
    )
    datums = ReducedDatum.objects.filter(target__in=target_list).order_by("timestamp")

    t2 = datetime.datetime.utcnow()
    logger.info('queryTools: Retrieved associated data for ' + str(len(target_list)) + ' Targets')
    logger.info('queryTools: N DB connections: ' + str(len(connection.queries)))
    logger.info('queryTools: Time taken: ' + str(t2 - t1))

    # Create microlensing event instances for the selected targets, associating all of the
    # data products for later use
    logger.info('queryTools: collating data on microlensing event set')
    target_data = {}
    for i, t in enumerate(target_list):
        logger.info('queryTools: evaluating target ' + t.name + ', '
                    + str(i) + ' out of ' + str(len(target_list)))

        mulens = MicrolensingEvent(t)
        mulens.set_extra_params(target_extras.filter(target=t))
        mulens.set_reduced_data(datums.filter(target=t))
        (status, reason) = mulens.check_need_to_fit()
        logger.info('queryTools: Need to fit: ' + repr(status) + ', reason: ' + reason)

        if check_need_to_fit and mulens.need_to_fit:
            target_data[t] = mulens

        if not check_need_to_fit:
            target_data[t] = mulens

    t3 = datetime.datetime.utcnow()
    logger.info('queryTools: Collated data for ' + str(len(target_data)) + ' targets in ' + str(t3 - t2))
    logger.info('queryTools: N DB connections: ' + str(len(connection.queries)))

    return target_data