from django.core.management.base import BaseCommand
from mop.brokers import ogle
import logging

logger = logging.getLogger(__name__)

class Command(BaseCommand):

    help = 'Downloads OGLE data for all events for a given years list'
    def add_arguments(self, parser):
        parser.add_argument('years', help='years you want to harvest, separated by ,')
        parser.add_argument('events', help='name of a specific event or all')

    def handle(self, *args, **options):
        logger.info('Starting run of OGLE event harvester')

        Ogle = ogle.OGLEBroker()

        list_of_targets = Ogle.fetch_alerts(years=[options['years']], events=options['events'])
        logger.info('Identified and ingested '+str(len(list_of_targets))+' target(s) from OGLE survey')

        Ogle.find_and_ingest_photometry(list_of_targets)

        logger.info('Completed run of OGLE event harvester')


   
