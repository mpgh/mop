from django.core.management.base import BaseCommand

from mop.asassn import ASASSNBroker

BROKER_URL = 'http://www.astronomy.ohio-state.edu/asassn/transients.html'
photometry = 'https://asas-sn.osu.edu/photometry'


class Command(BaseCommand):

    help = 'Downloads ASAS-SN data for all new microlensing events from the transient table'

    def add_arguments(self, parser):
        pass

    def handle(self):
        Asassn = ASASSNBroker()
        table = Asassn.retrieve_transient_table()
        list_of_events = Asassn.retrieve_microlensing_coordinates(table)
        list_of_targets = Asassn.fetch_alerts(list_of_events)
        Asassn.find_and_ingest_photometry(list_of_events, list_of_targets)
