from django.core.management.base import BaseCommand
from mop.brokers import ogle
import logging

logger = logging.getLogger(__name__)

class Command(BaseCommand):

    help = 'Downloads OGLE data for all events for a given years list'
    def add_arguments(self, parser):
        parser.add_argument('years', help='years you want to harvest, separated by ,')
        parser.add_argument('events', help='name of a specific event, all or an integer number')

    def handle(self, *args, **options):
        logger.info('Starting run of OGLE event harvester')

        Ogle = ogle.OGLEBroker()

        if not str(options['events']).isnumeric() or str(options['events']).lower() == 'all':
            list_of_targets = Ogle.fetch_alerts(years=[options['years']], events=options['events'])
        else:
            list_of_targets = Ogle.fetch_alerts(years=[options['years']], events='all')
        logger.info('Identified and ingested '+str(len(list_of_targets))+' target(s) from OGLE survey')

        # The following random selection is made to avoid the harvesting process taking so long
        # that the Kubernetes pod times out.  By randomizing the target selection, we ensure
        # that all targets should be updated quite often through more frequent runs of this harvester
        if str(options['events']).isnumeric():
            selected_targets = Ogle.select_random_targets(list_of_targets, ntargets=int(options['events']))
        else:
            selected_targets = Ogle.sort_target_list(list_of_targets)
        logger.info('Ingesting data from '+str(len(selected_targets))+' randomly-selected targets')

        Ogle.find_and_ingest_photometry(selected_targets)

        logger.info('Completed run of OGLE event harvester')


   
