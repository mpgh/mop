from django.core.management.base import BaseCommand
from mop.brokers import ogle



class Command(BaseCommand):

    help = 'Downloads OGLE data for all events for a given years list'
    def add_arguments(self, parser):
        parser.add_argument('years', help='years you want to harvest, spearted by ,')

    def handle(self, *args, **options):
        
        Ogle = ogle.OGLEBroker()
        list_of_targets = Ogle.fetch_alerts([options['years']])
        Ogle.find_and_ingest_photometry(list_of_targets)


   
