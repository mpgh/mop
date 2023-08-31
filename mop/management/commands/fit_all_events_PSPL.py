from django.core.management.base import BaseCommand
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target,TargetExtra
from astropy.time import Time
from mop.toolbox import fittools
from mop.brokers import gaia as gaia_mop
from mop.toolbox import logs
import numpy as np
import datetime
import random
import json
import datetime
from mop.management.commands.fit_need_events_PSPL import run_fit
import os

import logging

logger = logging.getLogger(__name__)

class Command(BaseCommand):

    help = 'Fit an event with PSPL and parallax, then ingest fit parameters in the db'

    def add_arguments(self, parser):

        parser.add_argument('events_to_fit', help='all, alive, need or [years]')
        parser.add_argument('--cores', help='Number of workers to use', default=os.cpu_count(), type=int)


    def handle(self, *args, **options):

        logger.info('Fitting all events')

        all_events = options['events_to_fit']

        if all_events == 'all':
            list_of_targets = Target.objects.filter()
        if all_events == 'alive':
            list_of_targets = Target.objects.filter(targetextra__in=TargetExtra.objects.filter(key='Alive', value=True))
        if all_events == 'need':

            four_hours_ago = Time(datetime.datetime.utcnow() - datetime.timedelta(hours=4)).jd
            list_of_targets = Target.objects.exclude(targetextra__in=TargetExtra.objects.filter(key='Last_Fit', value_lte=four_hours_ago))

        if all_events[0] == '[':

            years = all_events[1:-1].split(',')
            events = Target.objects.filter()
            list_of_targets = []
            for year in years:

                list_of_targets =  [i for i in events if year in i.name]

                list_of_targets = list(list_of_targets)
                random.shuffle(list_of_targets)

                logger.info('Found '+str(len(list_of_targets))+' targets to fit')

        for target in list_of_targets:
            # if the previous job has not been started by another worker yet, claim it

            print('Working on '+target.name)
            logger.info('Fitting data for '+target.name)
            try:
                if 'Gaia' in target.name:
                    gaia_mop.update_gaia_errors(target)


                if 'Microlensing' not in target.extra_fields['Classification']:
                    alive = False

                    extras = {'Alive':alive}
                    target.save(extras = extras)
                    logger.info(target.name+' not classified as microlensing')

                else:
                    result = run_fit(target, cores=options['cores'])

            except:
                logger.warning('Fitting event '+target.name+' hit an exception')
