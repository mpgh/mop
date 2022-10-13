from django.core.management.base import BaseCommand
from django.db.models import Q
from tom_dataproducts.models import ReducedDatum
from tom_targets.models import Target
from astropy.time import Time
from mop.toolbox import fittools
from mop.brokers import gaia as gaia_mop


import json
import numpy as np
import datetime
import os
from functools import reduce

class Command(BaseCommand):

    help = 'Identify microlensing events from Gaia alerts'

    def handle(self, *args, **options):

        # Retrieve a list of Gaia Targets that are flagged as Alive:
        targets = Target.objects.filter(name__includes='Gaia',
                                        targetextra__in=TargetExtra.objects.filter(key='Alive', value=True))


        # Evaluate each selected Target:
        for event in targets:

            # The expectation is that the lightcurve data for them will have a model
            # fit by a separate process, which will have stored the resulting model
            # parameters in the EXTRA_PARAMs for each Target.  Targets with no
            # fit parameters are ignored until they are model fitted.
            # Fitted targets will have their class set to microlensing by default
            if 'Blend_magnitude' in target.extra_fields.keys() \
                and 'u0' in target.extra_fields.keys() \
                and 't0' in target.extra_fields.keys() \
                and 'tE' in target.extra_fields.keys():

                # Retrieve the Gaia photometry for this Target:
                photometry = retrieve_target_photometry(target)

                # Test for an invalid blend magnitude:
                valid_blend_mag = True
                if target.extra_fields['Blend_magnitude'] == None:
                    valid_blend_mag = False

                # Test for a suspiciously large u0:
                valid_u0 = True
                if abs(target.extra_fields['u0']) > 0.5:
                    valid_u0 = False

                # Test for low-amplitude change in photometry:
                peak_mag = photometry[:,1].min()
                delta_mag = target.extra_fields['Baseline_magnitude'] - peak_mag
                valid_dmag = True
                if delta_mag < 0.5:
                    valid_dmag = False

                # If a target fails all three criteria, set its classification
                # to 'Unclassified variable'.  Note that TAP will consider scheduling
                # observations for any object with 'microlensing' in the
                # classification
                if not valid_blend_mag and not valid_u0 and not valid_dmag:
                    target.save(extras={'Classification': 'Unclassified variable'})

def retrieve_target_photometry(target):
    """Function to retrieve all available photometry for a target, combining
    all datasets.  Based on code by E. Bachelet."""

    datasets = ReducedDatum.objects.filter(target=target)
    time = [Time(i.timestamp).jd for i in datasets if i.data_type == 'photometry']

    phot = []
    for data in datasets:
        if data.data_type == 'photometry':
           try:
                phot.append([data.value['magnitude'],data.value['error'],data.value['filter']])

           except:
                # Weights == 1
                phot.append([data.value['magnitude'],1,data.value['filter']])


    photometry = np.c_[time,phot]

    return photometry
