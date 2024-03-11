from django.core.management.base import BaseCommand
from tom_targets.models import Target
from django.contrib.auth.models import Group, User
from guardian.shortcuts import assign_perm
from tom_targets.models import Target



class Command(BaseCommand):

    help = 'Add targets to the OMEGA users list'
    # def add_arguments(self, parser):
    #     parser.add_argument('target_name', help='name of the event to compute errors')


    def handle(self, *args, **options):

        omega_group = Group.objects.filter(name='OMEGA').first()
        if not omega_group:
            omega_group = Group.objects.filter(name='OMEGA-II').first()

        for target in Target.objects.all():
            assign_perm('tom_targets.view_target', omega_group, target)
            assign_perm('tom_targets.change_target', omega_group, target)
            #assign_perm('tom_targets.delete_target', omega_group, t)
