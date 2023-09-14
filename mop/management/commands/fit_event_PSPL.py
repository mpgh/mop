from django.core.management.base import BaseCommand
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target
from astropy.time import Time
from mop.toolbox import fittools
from mop.brokers import gaia as gaia_mop

import json
import numpy as np
import datetime
import os
import logging
from mop.management.commands.fit_need_events_PSPL import run_fit

logger = logging.getLogger(__name__)

class Command(BaseCommand):

    help = 'Fit an event with PSPL and parallax, then ingest fit parameters in the db'

    def add_arguments(self, parser):
        parser.add_argument('target_name', help='name of the event to fit')
        parser.add_argument('--cores', help='Number of workers to use', default=os.cpu_count(), type=int)


    def handle(self, *args, **options):

        target, created = Target.objects.get_or_create(name= options['target_name'])
        logger.info('Fitting single event: '+target.name)

        #try:

        if 'Gaia' in target.name:
            gaia_mop.update_gaia_errors(target)

        result = run_fit(target, cores=options['cores'])

        #except:
        #    logger.warning('Fitting event '+target.name+' hit an exception')
