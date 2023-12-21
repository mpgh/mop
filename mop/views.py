from django.shortcuts import redirect
from django.urls import reverse
from django.core.management import call_command
from io import StringIO
from tom_targets.views import TargetDetailView
from tom_targets.models import Target, TargetExtra
from tom_observations.models import ObservationRecord
from tom_observations.views import ObservationFilter
from mop.toolbox.TAP import set_target_sky_location
from mop.toolbox.obs_control import fetch_all_lco_requestgroups, parse_lco_requestgroups
from mop.forms import TargetClassificationForm
from django_filters.views import FilterView
from guardian.mixins import PermissionListMixin
from guardian.shortcuts import get_objects_for_user
from datetime import datetime, timedelta
from django.views.generic.list import ListView
import numpy as np

class MOPTargetDetailView(TargetDetailView):

    def get_context_data(self, *args, **kwargs):
        """
        Adds the ``TargetClassificationForm`` to the context and prepopulates the hidden fields.

        :returns: context object
        :rtype: dict
        """
        context = super().get_context_data(*args, **kwargs)
        context['class_form'] = TargetClassificationForm()
        return context

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
            utc_threshold = datetime.utcnow() - timedelta(days=7.0)
            utcnow = datetime.utcnow()
            obs_qs = ObservationRecord.objects.filter(scheduled_start__gt=utc_threshold,
                                                      scheduled_end__gt=utcnow)

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

        # Retrieve a list of pending observations from the LCO Portal
        response = fetch_all_lco_requestgroups()
        pending_obs = parse_lco_requestgroups(response, short_form=False, pending_only=False)

        pending_obs_list = []
        for target, obs_info in pending_obs.items():
            for config in obs_info:
                config['target'] = target
                config['filters'] = ','.join(config['filters'])
                config['exposure_times'] = ','.join([str(x) for x in config['exposure_times']])
                config['exposure_counts'] = ','.join([str(x) for x in config['exposure_counts']])
                pending_obs_list.append(config)

        context['pending_obs'] = pending_obs_list

        context['query_string'] = self.request.META['QUERY_STRING']

        return context

class PriorityTargetsView(ListView):
    template_name = 'priority_targets_list.html'
    paginate_by = 25
    strict = False
    model = TargetExtra
    permission_required = 'tom_targets.view_target'

    def get_context_data(self, *args, **kwargs):
        """Method to retrieve data on targets which TAP has prioritized

        Returns:
            context     dict    Context data for webpage template
        """

        # Configure context
        context = super().get_context_data(*args, **kwargs)
        context['target_count'] = context['paginator'].count

        # Retrieve two lists of Targets, which have been prioritized as stellar and Black Hole candidates
        if self.request.user.is_authenticated:

            # Query for matching TargetExtra entries returns a list of Target PKs
            qs_stars = TargetExtra.objects.filter(
                key='TAP_priority', float_value__gt=10.0
            ).exclude(
                value=np.nan
            ).exclude(
                value__exact=''
            ).exclude(
                value__exact='None'
            ).values_list('target').distinct()
            qs_bh = TargetExtra.objects.filter(
                key='TAP_priority_longtE', float_value__gt=10.0
            ).exclude(
                value=np.nan
            ).exclude(
                value__exact=''
            ).exclude(
                value__exact='None'
            ).filter(
                key='Alive', bool_value=True
            ).values_list('target').distinct()

            # Repackage the two lists to extract the parameters to display in the table.
            # This also checks to see if a Target is alive and not flagged as a known
            # variable before including it
            context['stellar_targets'] = self.extract_target_parameters(qs_stars, 'stellar')
            context['bh_targets'] = self.extract_target_parameters(qs_bh, 'bh')

        # If user is not logged in, return empty lists:
        else:
            context['stellar_targets'] = []
            context['bh_targets'] = []

        return context

    def extract_target_parameters(self, qs, target_category):

        key_list = ['t0', 't0_error', 'u0', 'u0_error', 'tE', 'tE_error', 'Mag_now', 'Baseline_magnitude']

        target_data = []
        priority = []
        for target_id in qs:
            target = Target.objects.filter(pk=target_id[0])[0]
            target_info = {'name': target.name, 'id': target_id[0]}
            if self.check_classification(target) and self.check_valid_target(target):
                if target_category == 'stellar':
                    target_info['priority'] = round(target.extra_fields['TAP_priority'],3)
                    # Not all entries have an uncertainty set, due to older versions of the code not storing it
                    try:
                        target_info['priority_error'] = round(target.extra_fields['TAP_priority_error'],3)
                    except KeyError:
                        target_info['priority_error'] = np.nan
                else:
                    target_info['priority'] = round(target.extra_fields['TAP_priority_longtE'],3)
                    try:
                        target_info['priority_error'] = round(target.extra_fields['TAP_priority_longtE_error'],3)
                    except KeyError:
                        target_info['priority_error'] = np.nan

                for key in key_list:
                    try:
                        target_info[key] = target.extra_fields[key]
                    except KeyError:
                        target_info[key] = np.nan
                    if key == 't0':
                        target_info[key] = round((target_info[key] - 2460000.0), 3)
                    else:
                        # Round floating point values where possible to save space in the table, catching
                        # NaN entries.  Skip in the event that the value is None or a string.
                        try:
                            if not np.isnan(target_info[key]):
                                target_info[key] = round(target_info[key], 3)
                        except:
                            pass

                if 'Outside HCZ' in target.extra_fields['Sky_location']:
                    target_data.append(target_info)
                    priority.append(target_info['priority'])

        # Sort the returned list
        priority = np.array(priority)
        idx = np.argsort(priority)[::-1]
        sorted_targets = [target_data[x] for x in idx]

        return sorted_targets


    def check_classification(self, target):
        """Method to check that the listed events are actually microlensing"""

        if 'microlensing' in str(target.extra_fields['Classification']).lower():
            criteria = [True]
        else:
            return False

        # Records True for each condition if it is NOT the classification
        bool_keys = ['is_YSO', 'is_QSO', 'is_galaxy']
        for key in bool_keys:
            if key in target.extra_fields.keys():
                if 'false' in str(target.extra_fields[key]).lower():
                    value = True
                else:
                    value = False
            else:
                value = True
            criteria.append(value)

        # Check for a TNS classification or name
        keys = ['TNS_name', 'TNS_class']
        for key in keys:
            if key in target.extra_fields.keys():
                if len(str(target.extra_fields[key]).lower().replace(' ','')) > 0:
                    value = True
                else:
                    value = False
                criteria.append(value)

        return all(criteria)

    def check_valid_target(self, target):
        """
        Method to verify that a Target is Alive and not flagged as a known variable before it is
        included in the Priority Targets table
        """

        if 'Alive' not in target.extra_fields.keys():
            return False

        if not target.extra_fields['Alive']:
            return False

        if 'is_YSO' in target.extra_fields.keys() and target.extra_fields['is_YSO']:
            return False

        if 'is_QSO' in target.extra_fields.keys() and target.extra_fields['is_QSO']:
            return False

        if 'is_galaxy' in target.extra_fields.keys() and target.extra_fields['is_galaxy']:
            return False

        return True