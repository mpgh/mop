from django import forms
from tom_targets.models import TargetList
from tom_observations import facility

class TargetClassificationForm(forms.Form):

    DEFAULT_CLASSES = [
        ('',''),
        ('Microlensing PSPL', 'Microlensing PSPL'),
        ('Microlensing binary', 'Microlensing binary'),
        ('Unclassified poor fit', 'Unclassified poor fit'),
        ('Variable star', 'Variable star'),
        ('Extra-galactic variable', 'Extra-galactic variable'),
    ]
    DEFAULT_CATEGORIES = [
        ('',''),
        ('Microlensing stellar/planet', 'Microlensing stellar/planet'),
        ('Microlensing long-tE', 'Microlensing long-tE'),
        ('Unclassified', 'Unclassified'),
        ('Eclipsing binary', 'Eclipsing binary'),
        ('Nova/supernova', 'Nova/supernova'),
        ('Stellar activity', 'Stellar activity'),
    ]

    classification = forms.ChoiceField(choices=DEFAULT_CLASSES, required=False, widget=forms.Select)
    category = forms.ChoiceField(choices=DEFAULT_CATEGORIES, required=False, widget=forms.Select)
    text_class = forms.CharField(required=False, widget=forms.TextInput)
    text_category = forms.CharField(required=False, widget=forms.TextInput)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class TargetSelectionForm(forms.Form):
    """
    Form for selecting the targets from a pre-existing TargetList
    """
    target_list = forms.ModelChoiceField(
        TargetList.objects.all(),
        required=False)
    facilities = [(x,x) for x in facility.get_service_classes()]
    observatory = forms.ChoiceField(required=True, choices=facilities)
    date = forms.DateTimeField(required=True, help_text='YYYY-MM-DD')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

