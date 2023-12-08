from django.core.management.base import BaseCommand
from django.db.models import Q
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target, TargetExtra
from astropy.time import Time
from mop.toolbox import fittools
from mop.brokers import gaia as gaia_mop
from mop.brokers import tns

import json
import numpy as np
import datetime
import os
from functools import reduce

from astropy.coordinates import SkyCoord
from astropy import units as u

from mop.toolbox import classifier_tools

import logging

logger = logging.getLogger(__name__)

class Command(BaseCommand):

    help = 'Identify microlensing events from Gaia alerts'

    def handle(self, *args, **options):
        classifier = 1

        logger.info('Gaia classifier started run - version check')

        # Retrieve a list of Gaia Targets that are flagged as Alive:
        targets = Target.objects.filter(name__contains='Gaia',
                                        targetextra__in=TargetExtra.objects.filter(key='Alive', value=True))

        logger.info('Found '+str(len(targets))+' alive Gaia targets')

        if classifier == 1:
            # Evaluate each selected Target:
            for event in targets:
                logger.info('Classifier evaluating ' + event.name)

                # The expectation is that the lightcurve data for them will have a model
                # fit by a separate process, which will have stored the resulting model
                # parameters in the EXTRA_PARAMs for each Target.  Targets with no
                # fit parameters are ignored until they are model fitted.
                # Fitted targets will have their class set to microlensing by default
                print(event.extra_fields['u0'], event.extra_fields['t0'], event.extra_fields['tE'],
                      event.extra_fields['Classification'], event.ra, event.dec)
                if event.extra_fields['u0'] != 0.0 \
                    and event.extra_fields['t0'] != 0.0 \
                    and event.extra_fields['tE'] != 0.0 \
                    and 'Microlensing' in event.extra_fields['Classification'] \
                    and event.ra != None and event.dec != None:

                    # Retrieve the Gaia photometry for this Target:
                    photometry = retrieve_target_photometry(event)

                    # Test for an invalid blend magnitude:
                    valid_blend_mag = classifier_tools.check_valid_blend(event.extra_fields['Blend_magnitude'])

                    # Test for a suspiciously large u0:
                    valid_u0 = classifier_tools.check_valid_u0(event.extra_fields['u0'])

                    # Test for low-amplitude change in photometry:
                    valid_dmag = classifier_tools.check_valid_dmag(event.extra_fields['Baseline_magnitude'], photometry)

                    # Test for suspicious reduced chi squared value
                    valid_chisq = classifier_tools.check_valid_chi2sq(event.extra_fields)

                    # Check target in catalogs
                    coord = SkyCoord(ra=event.ra, dec=event.dec, unit=(u.degree, u.degree), frame='icrs')

                    is_YSO, is_QSO, is_galaxy = False, False, False

                    if 'is_YSO' in event.extra_fields.keys():
                        if event.extra_fields['is_YSO'] == "True":
                            is_YSO = True
                        elif event.extra_fields['is_YSO'] == "False":
                            is_YSO = False
                    else:
                        is_YSO = classifier_tools.check_YSO(coord)
                        event.save(extras={'is_YSO': is_YSO})
                        logger.info(event.name + ': Checked if YSO for the first time.')

                    if 'is_QSO' in event.extra_fields.keys():
                        if event.extra_fields['is_QSO'] == "True":
                            is_QSO = True
                        if event.extra_fields['is_QSO'] == "False":
                            is_QSO = False
                    else:
                        is_QSO = classifier_tools.check_QSO(coord)
                        event.save(extras={'is_QSO': is_QSO})
                        logger.info(event.name + ': Checked if QSO for the first time.')

                    if 'is_galaxy' in event.extra_fields.keys():
                        if event.extra_fields['is_galaxy'] == "True":
                            is_galaxy = True
                        if event.extra_fields['is_galaxy'] == "False":
                            is_galaxy = False
                    else:
                        is_galaxy = classifier_tools.check_galaxy(coord)
                        event.save(extras={'is_galaxy': is_galaxy})
                        logger.info(event.name + ': Checked if SN for the first time.')

                    # Check target's TNS classification
                    if 'TNS_name' in event.extra_fields.keys() \
                        and event.extra_fields['TNS_name'] != 'None':
                        if event.extra_fields['TNS_class'] == 'None':
                            parameters = {
                                'objname': name
                            }
                            tns_object = tns.Custom_TNS
                            tns_class = tns.Custom_TNS.fetch_tns_class(tns_object, parameters)

                            for entry in tns_class:
                                if entry != None:
                                    event.extra_fields['TNS_class'] = str(entry)
                            logger.info(event.name + ': Known TNS, checked TNS class.')
                    else:

                        parameters = {
                            'ra': event.ra,
                            'dec': event.dec,
                            'radius': 1.0,
                            'units': 'arcsec',
                        }
                        tns_object = tns.Custom_TNS
                        tns_name = tns.Custom_TNS.fetch_tns_name(tns_object, parameters)

                        logger.info(event.name + ': Unknown TNS, checked TNS name.')

                        for name in tns_name:
                            parameters = {
                                'objname' : name
                            }
                            tns_class = tns.Custom_TNS.fetch_tns_class(tns_object, parameters)
                            logger.info(event.name + ': Unknown TNS, checked TNS class.')
                            if len(tns_name)>1 :
                                logger.info(event.name + ': More than one TNS entry within 1 arcsec radius.')
                                if(tns_class != None):
                                    event.save(extras={'TNS_name' : str(name),
                                                       'TNS_class' : str(tns_class)})
                                    logger.info(event.name + ': Saved a classified TNS entry.')
                            else:
                                event.save(extras={'TNS_name': str(name),
                                                   'TNS_class': str(tns_class)})
                                logger.info(event.name + ': Saved TNS entry.')

                    # Save classification based on catalogs and tests
                    if is_YSO:
                        event.save(extras={'Classification': 'Possible YSO'})
                        logger.info(event.name + ': set as a possible YSO')

                    if (event.extra_fields['TNS_class'] != 'None' \
                            and 'Other' not in event.extra_fields['TNS_class']):
                        # TNS has many classes, many of them related to SN, but not all.
                        # Classified microlensing events land in class "Other", however
                        # this class contains events that can have H alpha in emission.
                        event.save(extras={'Classification': 'Known TNS transient.'})
                        logger.info(event.name + ': set as a known SN')
                    elif is_QSO:
                        event.save(extras={'Classification': 'QSO'})
                        logger.info(event.name + ': set as a QSO')
                    elif is_galaxy:
                        # Some QSOs from Miliquas are in GLADE+ catalog,
                        # so we don't want them "misclassified"
                        event.save(extras={'Classification': 'Possible SN'})
                        logger.info(event.name + ': set as a possible SN')

                    # If a target fails all three criteria, set its classification
                    # to 'Unclassified variable'.  Note that TAP will consider scheduling
                    # observations for any object with 'microlensing' in the
                    # classification
                    if not valid_blend_mag or not valid_u0 or not valid_dmag:
                        event.save(extras={'Classification': 'Unclassified variable'})
                        logger.info(event.name+': Reset as unclassified variable')
                    if 'red_chi2' in event.extra_fields.keys():
                        if not valid_chisq:
                            event.save(extras={'Classification': 'Unclassified poor fit'})
                            logger.info(event.name+': Reset as unclassified poor fit')


        elif classifier == 2:
            for event in targets:
                event.save(extras={'Classification': 'Unclassified Gaia target'})
                logger.info(event.name+': Reset as unclassified Gaia target')


def retrieve_target_photometry(target):
    """Function to retrieve all available photometry for a target, combining
    all datasets.  Based on code by E. Bachelet."""

    datasets = ReducedDatum.objects.filter(target=target)

    phot = []
    time = []
    for data in datasets:
        if data.data_type == 'photometry':
            if 'magnitude' in data.value.keys():
                try:
                    phot.append([float(data.value['magnitude']),
                    float(data.value['error'])])
                    time.append(Time(data.timestamp).jd)
                except:
                    # Weights == 1
                    phot.append([float(data.value['magnitude']),
                    1])
                    time.append(Time(data.timestamp).jd)

    photometry = np.c_[time,phot]

    return photometry
