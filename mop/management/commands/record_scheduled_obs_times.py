from django.core.management.base import BaseCommand
from tom_observations.models import ObservationRecord

class Command(BaseCommand):

    help = "Command to back-fill the scheduled start and end times of observation groups"

    def handle(self, *args, **options):

        obs_list = ObservationRecord.objects.all()

        for obs_record in obs_list:
            if not obs_record.scheduled_start or not obs_record.scheduled_end:
                if 'windows' in obs_record.parameters['requests'][0].keys():
                    tstart = obs_record.parameters['requests'][0]['windows'][0]['start']
                    tend = obs_record.parameters['requests'][0]['windows'][0]['start']
                    for req in obs_record.parameters['requests']:
                        for window in req['windows']:
                            if window['start'] < tstart:
                                tstart = window['start']
                            if window['end'] > tend:
                                tend = window['end']
                    obs_record.scheduled_start = tstart
                    obs_record.scheduled_end = tend
                    obs_record.save()
            else:
                print(obs_record.target, obs_record.scheduled_start, obs_record.scheduled_end)