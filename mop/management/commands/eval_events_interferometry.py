from django.core.management.base import BaseCommand
from tom_targets.models import Target
from mop.toolbox import interferometry_prediction
import logging

logger = logging.getLogger(__name__)

class Command(BaseCommand):

    help = 'Evaluate one or more events as candidates for interferometry'

    def add_arguments(self, parser):
        parser.add_argument('target_selection', help='name of the event to fit or all')

    def handle(self, *args, **options):
        logger.info('INTERFERO: Starting evaluation of events as potential interferometry targets')

        # Select events to analyze:
        if str(options['target_selection']).lower() == 'all':
            target_list = Target.objects.all()
        else:
            qs = Target.objects.filter(name=options['target_selection'])

            if len(qs) == 0:
                logger.info('INTERFERO: No events found to match selection criterion: '+str(options['target_selection']))
            else:
                target_list = [qs[0]]

        logger.info('INTERFERO: Found '+str(len(target_list))
                    +' event(s) found to match selection criterion: ' + str(options['target_selection']))

        for target in target_list:
            interferometry_prediction.evaluate_target_for_interferometry(target)

        logger.info('INTERFERO: Completed evaluation of event set')