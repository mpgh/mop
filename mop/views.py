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
        obs_keys = ['facility', 'observation_type', 'observation_mode', 'start', 'end', 'period', 'ipp_value',
                    'c_1_instrument_type', 'c_1_ic_1_filter', 'c_1_ic_1_exposure_time',
                    'c_1_ic_2_filter', 'c_1_ic_2_exposure_time']
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
                    obs_data = {}
                    for key in obs_keys:
                        try:
                            obs_data[key] = obs.parameters[key]
                        except KeyError:
                            obs_data[key] = 'None'

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

