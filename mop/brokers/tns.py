from tom_alerts.brokers.tns import TNSBroker

class Custom_TNS(TNSBroker):
    def fetch_tns_class(cls, parameters):
        '''
        Modified version of fetch_alert from original TOM Toolkit.
        '''
        broker_feedback = ''

        data = {
            'api_key': settings.BROKERS['TNS']['api_key'],
            'data': json.dumps({
                'ra': parameters['ra'],
                'dec': parameters['dec'],
                'radius': parameters['radius'],
                'units': parameters['units'],
            }
            )
        }
        response = requests.post(TNS_SEARCH_URL, data, headers=cls.tns_headers())
        response.raise_for_status()
        transients = response.json()
        classes = []
        for transient in transients['data']['reply']:
            data = {
                'api_key': settings.BROKERS['TNS']['api_key'],
                'data': json.dumps({
                    'objname': transient['objname'],
                    'photometry': 1,
                    'spectroscopy': 0,
                }
                )
            }
            response = requests.post(TNS_OBJECT_URL, data, headers=cls.tns_headers())
            response.raise_for_status()
            alert = response.json()['data']['reply']
            tns_class = alert['object_type']['name']
            classes.append(tns_class)

        return classes
