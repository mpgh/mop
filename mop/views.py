from django.shortcuts import redirect
from django.urls import reverse
from django.core.management import call_command
from io import StringIO
from tom_targets.views import TargetDetailView
from tom_targets.models import Target, TargetExtra
from tom_observations.models import ObservationRecord
from tom_observations.views import ObservationFilter
from mop.toolbox.TAP import set_target_sky_location
from django_filters.views import FilterView
from guardian.mixins import PermissionListMixin
from guardian.shortcuts import get_objects_for_user
from datetime import datetime, timedelta
from django.views.generic.list import ListView

class MOPTargetDetailView(TargetDetailView):

    def get(self, request, *args, **kwargs):
        # Ensure that the target's location flag is set
        target = self.get_object()
        set_target_sky_location(target)

        fit_event = request.GET.get('fit_event', False)
        if fit_event:
            target_id = self.get_object().id
            target_name = self.get_object().name
            out = StringIO()
            call_command('fit_event_PSPL', target_name, cores=0, stdout=out)
            return redirect(reverse('tom_targets:detail', args=(target_id,)))

        TAP_event = request.GET.get('tap_event', False)
        if TAP_event:
            target_id = self.get_object().id
            target_name = self.get_object().name
            out = StringIO()
            call_command('run_TAP', target_name, stdout=out)
            return redirect(reverse('tom_targets:detail', args=(target_id,)))

        return super().get(request, *args, **kwargs)

class ActiveObsView(ListView):
    template_name = 'active_obs_list.html'
    paginate_by = 25
    strict = False
    model = ObservationRecord
    filterset_class = ObservationFilter
    permission_required = 'tom_targets.view_target'
    ordering = ['-created']

    def get_context_data(self, *args, **kwargs):
        """Method to retrieve data on targets which TAP has recommended for observation within the last 24hrs,
        and their associated observation records

        Returns:
            context     dict    Context data for webpage template
        """
        # Parameters to extract from the observation records
        target_keys = ['tE', 't0', 'Mag_now', 'Category', 'TAP_priority', 'TAP_priority_longtE']

        context = super().get_context_data(*args, **kwargs)

        # Provide a count of the number of targets returned to enable pagination of the results
        context['target_count'] = context['paginator'].count

        # Retrieve a list of Targets for which observations have recently been requested
        if self.request.user.is_authenticated:
            utc24 = datetime.utcnow() - timedelta(days=1.0)
            obs_qs = ObservationRecord.objects.filter(created__gt=utc24)

            # Build lists of the ObservationRecords and Targets, while respecting
            # user access permissions
            targets = {}
            for obs in obs_qs:
                if self.request.user.has_perm('tom_observations.view_observation'):
                    if 'requests' in obs.parameters.keys():
                        for request in obs.parameters['requests']:
                            obs_data = {}
                            obs_data['start'] = request['windows'][0]['start']
                            obs_data['end'] = request['windows'][0]['end']
                            obs_data['facility'] = 'LCO'
                            obs_data['observation_type'] = obs.parameters['observation_type']
                            obs_data['ipp_value'] = obs.parameters['ipp_value']
                            for c, config in enumerate(request['configurations']):
                                obs_data['c_' + str(c+1) + '_instrument_type'] = config['instrument_type']
                                for ic,inst_config in enumerate(config['instrument_configs']):
                                    obs_data['c_' + str(c+1) + '_ic_' + str(ic+1) + '_filter'] = inst_config['optical_elements']['filter']
                                    obs_data['c_' + str(c+1) + '_ic_' + str(ic+1) + '_exposure_time'] = inst_config['exposure_time']
                                    obs_data['c_' + str(c+1) + '_ic_' + str(ic+1) + '_exposure_count'] = inst_config['exposure_count']

                        if obs.target.name in targets.keys():
                            target_data = targets[obs.target.name]
                        else:
                            target_data = {'name': obs.target.name, 'obs_list': []}
                            for key in target_keys:
                                try:
                                    target_data[key] = obs.target.extra_fields[key]
                                except KeyError:
                                    target_data[key] = 'None'

                        target_data['obs_list'].append(obs_data)
                        targets[obs.target.name] = target_data

            context['targets'] = targets.values()
            
        else:
            context['targets'] = []
            context['obs_list'] = []

        context['query_string'] = self.request.META['QUERY_STRING']

        return context

