from django.core.management.base import BaseCommand
from tom_alerts.models import BrokerQuery

class Command(BaseCommand):
    help = "Command to clear the table of broker queries"

    def handle(self, *args, **options):
        qs = BrokerQuery.objects.all()

        print('Found ' + str(qs.count()) + ' stored broker queries to remove')

        for query in qs:
            print('Removing query ' + query.name + ' broker ' + query.broker)
            query.delete()
            