from django.shortcuts import redirect
from django.urls import reverse
from django.core.management import call_command
from io import StringIO
from tom_targets.views import TargetDetailView


class MOPTargetDetailView(TargetDetailView):

    def get(self, request, *args, **kwargs):
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
