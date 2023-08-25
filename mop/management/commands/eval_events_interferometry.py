from django.core.management.base import BaseCommand
from tom_targets.models import Target
from tom_dataproducts.models import ReducedDatum
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
            logger.info('INTERFERO: Evaluating event '+target.name)

            # If sufficient information is available for this target, extract it's Gaia photometry and
            # estimate JHK-band photometry with uncertainties.  This requires a valid model
            if target.extra_fields['Gmag'] > 0.0 and target.extra_fields['BP-RP'] > 0.0 \
                and target.extra_fields['u0'] > 0.0:
                (Jtarget, Htarget, Ktarget) = interformetry_prediction.convert_Gmag_to_JHK(target.extra_fields['Gmag'],
                                                                     target.extra_fields['BP-RP'])
                logger.info('-> Calculated J='+str(Jtarget)+'mag, H='+str(Htarget)+'mag K='+str(Ktarget))

                Gmag_error = interformetry_prediction.estimate_target_Gaia_phot_uncertainties(
                    target.extra_fields['Gmag'], target.extra_fields['u0'], 0.01)

                target.save(extras = {
                    'Gmag_error': Gmag_error,
                    'Mag_peak_J': Jtarget,
                    'Mag_peak_H': Htarget,
                    'Mag_peak_K': Ktarget
                })
                logger.info('-> Calculated uncertainties Gmag='+str(target.extra_fields['Gmag'])
                            +'+/-'+str(Gmag_error)+'mag')

                # Search the Gaia catalog for all stars neighbouring the target
                star_catalog = gaia.query_gaia_dr3(target, radius=neighbour_radius)
                neighbours = interformetry_prediction.find_companion_stars(target, star_catalog)
                logger.info('-> Identified '+str(len(neighbours)-1)+' neighbouring stars')

                # Estimate the JHK photometry for all neighbouring stars
                (J, H, K) = interformetry_prediction.convert_Gmag_to_JHK(neighbours['Gmag'],
                                                                         neighbours['BP-RP'])
                logger.info('-> Computed JHK photometry for neighbouring stars')

                # Repackage data into a convenient form for storage
                datum = {
                    'Gaia_Source_ID': [x for x in neighbours['Source']],
                    'Gmag': [x for x in neighbours['Gmag']],
                    'BP-RP': [x for x in neighbours['BP-RP']],
                    'Jmag': [x for x in J],
                    'Hmag': [x for x in H],
                    'Kmag': [x for x in K],
                    'Separation':  [x for x in neighbours['Separation']],
                    }
                tnow = Time.now()
                rd, created = ReducedDatum.objects.get_or_create(
                    timestamp=tnow.to_datetime(timezone=TimezoneInfo()),
                    value=datum,
                    source_name='Interferometry_predictor',
                    source_location=target.name,
                    data_type='photometry',
                    target=target)
                logger.info('-> Stored neighbouring star data in MOP')

                # Determine whether or not this is a candidate target for interferometry
                (mode, guide) = interformetry_prediction.interferometry_decision(Ktarget, np.array(K))
                target.save(extras={
                    'Interferometry_mode': mode,
                    'Inteferometry_guide_star': guide
                })
                logger.info('-> Evaluation for interferometry for '+target.name+': '+str(mode)+' guide='+str(guide))

            else:
                logger.info('-> target '+target.name
                            +' missing one of Gmag, BP-RP or u0, cannot evaluate for interferometry')

        logger.info('INTERFERO: Completed evaluation of event set')