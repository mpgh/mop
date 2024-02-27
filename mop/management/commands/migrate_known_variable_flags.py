from django.core.management.base import BaseCommand
from tom_targets.models import Target, TargetExtra
from mop.toolbox import classifier_tools
from astropy import units as u
from astropy.coordinates import SkyCoord
import logging

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    help = 'Evaluate whether each Target is a known YSO, QSO or galaxy, and set the extra_field parmaeters accordingly'

    def handle(self, *args, **options):

        # Parameters to check, their mapping to the new parameter set, and the function
        # that should be used to verify the setting
        param_mapping = {
            'is_YSO': {'new_key': 'YSO', 'check_function': classifier_tools.check_YSO},
            'is_QSO': {'new_key': 'QSO', 'check_function': classifier_tools.check_QSO},
            'is_galaxy': {'new_key': 'galaxy', 'check_function': classifier_tools.check_galaxy}
        }

        # Walk over all Targets
        qs = Target.objects.all()
        nt = qs.count()

        for i,t in enumerate(qs):
            logger.info('Evaluating target ' + t.name + ', ' + str(i) + ' out of ' + str(nt))

            for old_key, config in param_mapping.items():
                logger.info('Checking old key ' + old_key)

                extra_pars = t.extra_fields
                try:
                    s = SkyCoord(t.ra, t.dec, frame='icrs', unit=(u.deg, u.deg))

                    # If the new key is already set, we assume it is correct
                    if config['new_key'] not in extra_pars.keys():
                        if old_key in extra_pars.keys() and config['new_key'] not in extra_pars.keys():
                            status, new_value = self.interpret_entry(extra_pars[old_key])
                            logger.info(' -> Found existing entry ' + old_key + '=' + str(new_value) + ', status=' + status)

                            # If the existing value is unset, or has an invalid entry, double check
                            # what it should be
                            if status == 'unknown':
                                new_value = config['check_function'](s)
                                logger.info(' -> Unknown status, so checked catalog with result ' + str(new_value))

                        elif old_key not in extra_pars.keys() and config['new_key'] not in extra_pars.keys():
                            new_value = config['check_function'](s)
                            logger.info(' -> No existing entry, so checked catalog with result ' + str(new_value))

                        ep = TargetExtra.objects.create(
                            target=t,
                            key=config['new_key'],
                            value=new_value
                        )
                        ep.save()
                        logger.info(' -> Created new key ' + config['new_key'] + '=' + str(new_value))

                    else:
                        logger.info(' -> Found existing new key entry for ' + config['new_key']
                                    + '=' + str(extra_pars[config['new_key']]))

                    # Remove the old key entries if available
                    if old_key in extra_pars.keys():
                        ep = TargetExtra.objects.get(target=t, key=old_key)
                        ep.delete()
                        logger.info(' -> Removed old key ' + old_key + ' entry ' + repr(ep))

                except TypeError:
                    logger.info('Target ' + t.name + ' has no valid coordinates, skipping')
                    
                #opt = input('Continue?')
            #if opt != 'y':
            #    exit()

    def interpret_entry(self, value):
        """Method to interpret an existing entry for one of the old-style TargetExtra parameters
        for is_YSO, is_QSO and is_galaxy.  These are defined to be strings, but may be empty,
        set to 'true' or 'false' or have boolean values or None."""

        if isinstance(value, str):
            # If the entry is an empty string, we don't know what the setting should be
            if len(value) == 0:
                return 'unknown', None

            # If the entry string contains either 'true' or 'false', assume it has been set
            # and return the corresponding boolean.  We can't interpret any other entry
            # so return Unknown to force the code to double-check
            else:
                if 'false' in str(value).lower():
                    return 'set', False
                elif 'true' in str(value).lower():
                    return 'set', True
                else:
                    return 'unknown', None

        elif isinstance(value, bool):
            # If the value is a boolean, then it must have already been set, so return the value
            return 'set', value

        elif isinstance(value, None):
            # If the entry is None, the value hasn't been set
            return 'unknown', None

        else:
            # In any other case, go and double-check the entry
            return 'unknown', None