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

from astropy.coordinates import SkyCoord
from astropy import units as u

from mop.toolbox.classifier_tools import check_YSO, check_QSO, check_SN

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

                    # Check target in catalogs
                    coord = SkyCoord(ra=event.ra, dec=event.dec, unit=(u.degree, u.degree), frame='icrs')

                    is_YSO, is_QSO, is_SN = False, False, False

                    if 'is_YSO' in event.extra_fields.keys():
                        if event.extra_fields['is_YSO']:
                            is_YSO = True
                    else:
                        is_YSO = check_YSO(coord)
                        event.save(extras={'is_YSO': is_YSO})
                        log.info(event.name + ': Checked if YSO for the first time.')

                    if 'is_QSO' in event.extra_fields.keys():
                        if event.extra_fields['is_QSO']:
                            is_QSO = True
                    else:
                        is_QSO = check_QSO(coord)
                        event.save(extras={'is_QSO': is_QSO})
                        log.info(event.name + ': Checked if QSO for the first time.')

                    if 'is_SN' in event.extra_fields.keys():
                        if event.extra_fields['is_SN']:
                            is_SN = True
                    else:
                        is_SN = check_SN(coord)
                        event.save(extras={'is_SN': is_SN})
                        log.info(event.name + ': Checked if SN for the first time.')

                    # Save classification based on catalogs
                    if is_YSO:
                        event.save(extras={'Classification': 'Possible YSO'})
                        log.info(event.name + ': set as a possible YSO')

                    if is_QSO:
                        event.save(extras={'Classification': 'QSO'})
                        log.info(event.name + ': set as a QSO')
                    elif is_SN:
                        # Some QSOs from Miliquas are in GLADE+ catalog,
                        # so we don't want them "misclassified"
                        event.save(extras={'Classification': 'Possible SN'})
                        log.info(event.name + ': set as a possible SN')

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
