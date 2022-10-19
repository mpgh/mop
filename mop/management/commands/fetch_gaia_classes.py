from django.core.management.base import BaseCommand
from tom_targets.models import Target, TargetExtra
from os import path

class Command(BaseCommand):
    help = 'Retrieve a list of all alive Gaia alerts and their current classifications'

    def add_arguments(self, parser):
        parser.add_argument('output_file', help='Path to the output data file')

    def handle(self, *args, **options):

        # Retrieve a list of Gaia Targets that are flagged as Alive:
        targets = Target.objects.filter(name__contains='Gaia',
                                        targetextra__in=TargetExtra.objects.filter(key='Alive', value=True))

        # Log output to file:
        f = open(options['output_file'], 'w')
        f.write('# Event   Classification\n')

        # Record the status of all matching Targets:
        for event in targets:
            f.write(event.name+'  '+event.extra_fields['Classification']+'\n')

        f.close()
