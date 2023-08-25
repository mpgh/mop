from django.core.management.base import BaseCommand
from tom_targets.models import Target
from mop.brokers import gaia
from astropy.coordinates import Angle
import astropy.units as u
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
            qs = Target.objects.filter(name=options['target_selction'])

            if len(qs) == 0:
                logger.info('INTERFERO: No events found to match selection criterion: '+str(options['target_selection']))
            else:
                target_list = [qs[0]]

        logger.info('INTERFERO: Found '+str(len(target_list))
                    +' event(s) found to match selection criterion: ' + str(options['target_selection']))


        # Search radius is set by the interferometry requirement to have bright neighbouring stars
        # that are accessible to the GRAVITY instrument
        neighbour_radius = Angle(20.0/3600.0, "deg")

        for target in target_list:
            # Search the Gaia catalog for all stars neighbouring the target
            star_catalog = gaia.query_gaia_dr3(target, radius=neighbour_radius)

            # Estimate the JHK photometry for all neighbouring stars
