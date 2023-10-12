from tom_targets.models import Target
def fetch_extra_param(target, key):
    string_keys = ['Classification', 'Category', 'Observing_mode', 'Sky_location',
                   'Gaia_Source_ID', 'Interferometry_mode', 'Mag_now_passband']
    if key in target.extra_fields.keys():
        if key in string_keys:
            value = target.extra_fields[key]
        else:
            try:
                value = float(target.extra_fields[key])
            except ValueError:
                value = None
    elif key in target.tags.keys():
        if key in string_keys:
            value = target.extra_fields[key]
        else:
            try:
                value = float(target.tags[key])
            except ValueError:
                value = None
    else:
        value = None

    return value