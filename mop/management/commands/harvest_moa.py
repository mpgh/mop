from django.core.management.base import BaseCommand
from mop.brokers import moa
from mop.brokers import gaia as gaia_mop
import logging

logger = logging.getLogger(__name__)

class Command(BaseCommand):

    help = 'Downloads MOA data for all events for a given years list'

    def add_arguments(self, parser):
        parser.add_argument('years', help='years you want to harvest, spearted by ,')

    def handle(self, *args, **options):
        
        Moa = moa.MOABroker()
        (list_of_targets, new_targets) = Moa.fetch_alerts('./data/',[options['years']])
        logger.info('MOA HARVESTER: Found '+str(len(list_of_targets))+' targets')

        Moa.find_and_ingest_photometry(list_of_targets)
        logger.info('MOA HARVESTER: Ingested photometry for MOA targets')

        for target in new_targets:
            gaia_mop.fetch_gaia_dr3_entry(target)
        logger.info('MOA HARVESTER: Retrieved Gaia photometry for MOA targets')
