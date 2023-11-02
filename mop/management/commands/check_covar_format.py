from django.core.management.base import BaseCommand
from tom_targets.models import Target, TargetExtra
from mop.management.commands import run_TAP
class Command(BaseCommand):
    help = "Verify the format used to store the covarience matrix"

    def handle(self, *args, **options):
        qs = Target.objects.all()

        for target in qs:
            print(target.name)
            if 'Fit_covariance' in target.extra_fields.keys():
                covar = run_TAP.load_covar_matrix(target.extra_fields['Fit_covariance'])
                print(target.name + ' covar=' + repr(covar))