from django.core.management.base import BaseCommand
from tom_observations.models import ObservationRecord
from datetime import datetime
import pytz

class Command(BaseCommand):

    help = "Command to back-fill the scheduled start and end times of observation groups"

    def handle(self, *args, **options):

        tz = pytz.timezone('utc')

        obs_list = ObservationRecord.objects.all()

        for obs_record in obs_list:
            if not obs_record.scheduled_start or not obs_record.scheduled_end:
                tstart = None
                tend = None
                if 'requests' in obs_record.parameters:
                    # Handle current (2023) format of JSON parameter field
                    if 'windows' in obs_record.parameters['requests'][0].keys():
                        tstart = obs_record.parameters['requests'][0]['windows'][0]['start']
                        tend = obs_record.parameters['requests'][0]['windows'][0]['start']
                        for req in obs_record.parameters['requests']:
                            for window in req['windows']:
                                if window['start'] < tstart:
                                    tstart = window['start']
                                if window['end'] > tend:
                                    tend = window['end']

                # Handle older format of parameter dictionary
                else:
                    if 'start' in obs_record.parameters and 'end' in obs_record.parameters:
                        tstart = obs_record.parameters['start']
                        tend = obs_record.parameters['end']

                if tstart and tend:
                    tstart = self.convert_to_datetime(tstart, tz)
                    tend = self.convert_to_datetime(tend, tz)

                # This is a separate if clause to handle the case where the conversion to datetime fails
                if tstart and tend:
                    obs_record.scheduled_start = tstart
                    obs_record.scheduled_end = tend
                    obs_record.save()

            else:
                print(obs_record.target, obs_record.scheduled_start, obs_record.scheduled_end)

    def convert_to_datetime(self, date_string, tz):
        """Convert dates and times from strings to datetime objects, when they may be in different string formats"""

        # Check whether the input is already a datetime object, and return if it is
        if type(date_string) == type(datetime.utcnow()):
            t = date_string.replace(tzinfo=tz)
            return t

        # Parse date/time strings in different formats
        t = None
        if 'T' in date_string:
            formats = ['%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S']
        else:
            formats = ['%Y-%m-%dT %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S']

        for format in formats:
            try:
                t = datetime.strptime(date_string, format)
            except ValueError:
                pass

        # Ensure the timezone information is set correctly
        if t:
            t = t.replace(tzinfo=tz)

        return t