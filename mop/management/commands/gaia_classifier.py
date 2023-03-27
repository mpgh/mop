from django.core.management.base import BaseCommand
from django.db.models import Q
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target, TargetExtra
from astropy.time import Time
from mop.toolbox import fittools
from mop.brokers import gaia as gaia_mop
from mop.toolbox import logs

import json
import numpy as np
import datetime
import os
from functools import reduce

class Command(BaseCommand):

    help = 'Identify microlensing events from Gaia alerts'

    def handle(self, *args, **options):
        classifier = 1

        # Start logging process:
        log = logs.start_log()
        log.info('Gaia classifier started run - version check')

        # Retrieve a list of Gaia Targets that are flagged as Alive:
        targets = Target.objects.filter(name__contains='Gaia',
                                        targetextra__in=TargetExtra.objects.filter(key='Alive', value=True))

        log.info('Found '+str(len(targets))+' alive Gaia targets')

        if classifier == 1:
            # Evaluate each selected Target:
            for event in targets:

                # The expectation is that the lightcurve data for them will have a model
                # fit by a separate process, which will have stored the resulting model
                # parameters in the EXTRA_PARAMs for each Target.  Targets with no
                # fit parameters are ignored until they are model fitted.
                # Fitted targets will have their class set to microlensing by default
                if event.extra_fields['u0'] != 0.0 \
                    and event.extra_fields['t0'] != 0.0 \
                    and event.extra_fields['tE'] != 0.0 \
                    and 'Microlensing' in event.extra_fields['Classification']:

                    # Retrieve the Gaia photometry for this Target:
                    photometry = retrieve_target_photometry(event)

                    # Test for an invalid blend magnitude:
                    valid_blend_mag = True
                    if event.extra_fields['Blend_magnitude'] == None \
                        or event.extra_fields['Blend_magnitude'] == 0.0:
                        valid_blend_mag = False

                    # Test for a suspiciously large u0:
                    valid_u0 = True
                    if abs(event.extra_fields['u0']) > 0.5:
                        valid_u0 = False

                    # Test for low-amplitude change in photometry:
                    if len(photometry) > 0:
                        peak_mag = photometry[:,1].min()
                        delta_mag = event.extra_fields['Baseline_magnitude'] - peak_mag
                        valid_dmag = True
                        if delta_mag < 0.5:
                            valid_dmag = False
                    else:
                        valid_dmag = False

                    # Test for suspicious reduced chi squared value
                    if 'red_chi2' in event.extra_fields.keys():
                        valid_chisq = True
                        if event.extra_fields['red_chi2'] > 50.0 \
                            or event.extra_fields['red_chi2'] < 0.0:
                            valid_chisq = False

                    # If a target fails all three criteria, set its classification
                    # to 'Unclassified variable'.  Note that TAP will consider scheduling
                    # observations for any object with 'microlensing' in the
                    # classification
                    if not valid_blend_mag or not valid_u0 or not valid_dmag:
                        event.save(extras={'Classification': 'Unclassified variable'})
                        log.info(event.name+': Reset as unclassified variable')
                    if 'red_chi2' in event.extra_fields.keys():
                        if not valid_chisq:
                            event.save(extras={'Classification': 'Unclassified poor fit'})
                            log.info(event.name+': Reset as unclassified poor fit')


        elif classifier == 2:
            for event in targets:
                event.save(extras={'Classification': 'Unclassified Gaia target'})
                log.info(event.name+': Reset as unclassified Gaia target')

        logs.stop_log(log)

def retrieve_target_photometry(target):
    """Function to retrieve all available photometry for a target, combining
    all datasets.  Based on code by E. Bachelet."""

    datasets = ReducedDatum.objects.filter(target=target)
    time = [Time(i.timestamp).jd for i in datasets if i.data_type == 'photometry']

    phot = []
    for data in datasets:
        if data.data_type == 'photometry':
           try:
                phot.append([float(data.value['magnitude']),
                             float(data.value['error'])])

           except:
                # Weights == 1
                phot.append([float(data.value['magnitude']),
                             1])


    photometry = np.c_[time,phot]

    return photometry
