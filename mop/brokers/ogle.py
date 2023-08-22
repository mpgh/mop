from django.core.management.base import BaseCommand
from django.core.exceptions import ObjectDoesNotExist
from tom_alerts.alerts import GenericBroker, GenericQueryForm
from django import forms
from tom_targets.models import Target
from tom_observations import facility
from tom_dataproducts.models import ReducedDatum

from astropy.coordinates import SkyCoord
import astropy.units as unit
import ftplib
import os
import numpy as np
import requests

from astropy.time import Time, TimezoneInfo


BROKER_URL = 'ftp.astrouw.edu.pl'


class OGLEQueryForm(GenericQueryForm):
    target_name = forms.CharField(required=False)
    cone = forms.CharField(
        required=False,
        label='Cone Search',
        help_text='RA,Dec,radius in degrees'
    )

    def clean(self):
        if len(self.cleaned_data['target_name']) == 0 and \
                        len(self.cleaned_data['cone']) == 0:
            raise forms.ValidationError(
                "Please enter either a target name or cone search parameters"
                )

class OGLEBroker(GenericBroker):
    name = 'OGLE'
    form = OGLEQueryForm

    def fetch_alerts(self, years = []):
        """Fetch data on microlensing events discovered by OGLE"""
        
        # Read the lists of events for the given years
        ogle_events = self.fetch_lens_model_parameters(years)

        #ingest the TOM db
        list_of_targets = self.ingest_events(ogle_events)

        return list_of_targets

    def fetch_lens_model_parameters(self, years):
        """Method to retrieve the text file of the model parameters for fits by the OGLE survey"""

        events = {}
        for year in years:
            par_file_url = os.path.join('https://www.astrouw.edu.pl/ogle/ogle4/ews',year,'lenses.par')
            response = requests.request('GET', par_file_url)
            if response.status_code == 200:
                for line in response.iter_lines():
                    line = str(line)
                    if 'StarNo' not in line and len(line) > 5:      # Skip the file header
                        entries = line.split()
                        name = 'OGLE-'+entries[0].replace("b'","")
                        ra = entries[3]
                        dec = entries[4]
                        events[name] = (ra,dec)

        return events

    def ingest_events(self, ogle_events):
        """Function to ingest the targets into the TOM database"""

        list_of_targets = []

        for event_name, event_params in ogle_events.items():
            s = SkyCoord(event_params[0], event_params[1], unit=(unit.hourangle, unit.deg), frame='icrs')
            target, created = Target.objects.get_or_create(name=event_name, ra=s.ra.deg, dec=s.dec.deg,
                                                           type='SIDEREAL', epoch=2000)
            if created:
                target.save()

            list_of_targets.append(target)

        return list_of_targets

    def find_and_ingest_photometry(self, targets):
        
        ftp_tunnel = ftplib.FTP( BROKER_URL ) 
        ftp_tunnel.login()
        ftp_file_path = os.path.join( 'ogle', 'ogle4', 'ews' )
        ftp_tunnel.cwd(ftp_file_path)

        previous_year = targets[0].name.split('-')[1]
        ftp_tunnel.cwd(previous_year)
        for target in targets:
            
            year = target.name.split('-')[1]
            event = target.name.split('-')[2]+'-'+target.name.split('-')[3]

            if year == previous_year:
        
               pass
        
            else:

               ftp_tunnel.cwd('../../')
               ftp_tunnel.cwd(year)
            
            ftp_tunnel.cwd(event.lower())
            ftp_tunnel.retrbinary('RETR phot.dat',open('./data/ogle_phot.dat', 'wb').write)
            photometry = np.loadtxt('./data/ogle_phot.dat')
            photometry = photometry[photometry[:,0].argsort()[::-1],] #going backward to save time on ingestion
            ftp_tunnel.cwd('../')
            for index,point in enumerate(photometry):
        
                jd = Time(point[0], format='jd', scale='utc')
                jd.to_datetime(timezone=TimezoneInfo())
                data = {   'magnitude': point[1],
                           'filter': 'OGLE_I',
                           'error': point[2]
                       }
                rd, created = ReducedDatum.objects.get_or_create(
                timestamp=jd.to_datetime(timezone=TimezoneInfo()),
                value=data,
                source_name='OGLE',
                source_location=target.name,
                data_type='photometry',
                target=target)
               
                if created:

                    rd.save()

                else:
                    # photometry already there (I hope!)
                    break
            os.remove('./data/ogle_phot.dat')

    def to_generic_alert(self, alert):
        pass
  
