import factory

from tom_dataproducts.models import ReducedDatum

class ReducedDatumFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = ReducedDatum

    timestamp = factory.Faker('datetime')
    value = factory.Faker('pydict')
    source_name = 'MOP'
    source_location = factory.Faker('pystr')
    data_type = 'lightcurve'
