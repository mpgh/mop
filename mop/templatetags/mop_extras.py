from django import template
from django.conf import settings
from django.core.paginator import Paginator
from datetime import datetime

from plotly import offline
import plotly.graph_objs as go

from tom_dataproducts.models import DataProduct, ReducedDatum
from tom_dataproducts.processors.data_serializers import SpectrumSerializer
from mop.toolbox import utilities
from mop.forms import TargetClassificationForm
from astropy.time import Time
from datetime import datetime
import numpy as np
import logging

logger = logging.getLogger(__name__)

register = template.Library()

@register.inclusion_tag('tom_dataproducts/partials/photometry_for_target.html')
def mop_photometry(mulens):
    """
    Renders a photometric plot for a target.
    This templatetag requires all ``ReducedDatum`` objects with a data_type of ``photometry`` to be structured with the
    following keys in the JSON representation: magnitude, error, filter
    """

    logger.info('MOP PHOTOMETRY: Got ' + str(mulens.ndata) + ' datasets for target ' + mulens.name)

    plot_data = [
        go.Scatter(
            x=dataset_values[:,0] - 2460000.0,
            y=dataset_values[:,1], mode='markers',
            name=dataset_name,
            error_y=dict(
                type='data',
                array=dataset_values[:,2],
                visible=True
            )
        ) for dataset_name, dataset_values in mulens.datasets.items()]

    layout = go.Layout(
        yaxis=dict(autorange='reversed'),
        height=600,
        width=700,
             
    )

    fig = go.Figure(data=plot_data, layout=layout)
    current_time =  Time.now().jd-2460000
    fig.add_shape(
        # Line Vertical
        dict(
             type="line",
             x0=current_time,
             y0=0,
             x1=current_time,
             y1=1,
             yref='paper',
             layer='below',
             line=dict(
                 color="Black",
                 width=1,
                 dash='dash',
                 )

    ))

    ### Try to plot model if exist
    if mulens.existing_model:

        fig.add_trace(go.Scatter(x = np.array(mulens.existing_model.value['lc_model_time']) - 2460000,
                                 y = np.array(mulens.existing_model.value['lc_model_magnitude']),
                                 mode = 'lines',
                                 name = 'Model',
                                 opacity = 0.5,
                                 line = dict(color = 'rgb(128,128,128)',
                                             width = 5,),
                                 )
                      )

    fig.update_layout(

    annotations=[
        dict(
             x=current_time,
             xanchor="left",
             y=0.05,
             yref="paper",
             text="JD now : "+str(np.round(current_time,3))+" ("+str(Time.now().value).split(' ')[0]+")",
             showarrow=False,
             textangle=-90,)
    ],		
             xaxis_title="HJD-2460000",
             yaxis_title="Mag",
    )

    return {
        'target': mulens.target,
        'plot': offline.plot(fig, output_type='div', show_link=False)
    }


@register.inclusion_tag('tom_dataproducts/partials/interferometry_data.html')
def interferometry_data(mulens):
    context = {}
    t1 = datetime.utcnow()
    logger.info('MOP INTERFEROMETRY: ' + str(t1))
    utilities.checkpoint()

    # Gather extra_param key values for target's Gaia catalog data
    #u0 = utilities.fetch_extra_param(target, 'u0')
    #u0_error = utilities.fetch_extra_param(target, 'u0_error')
    u0 = float(mulens.u0)
    u0_error = float(mulens.u0_error)
    if u0 == None or u0_error == None or u0 == 0.0 or u0_error == 0.0:
        context['model_valid'] = False
    else:
        context['model_valid'] = True

    def clean_key_string(key):
        return key.replace('[','').replace(']','').replace('/','_').replace('(','_').replace(')','_').replace('-','_')

    key_list = ['Gaia_Source_ID',
                'Gmag', 'Gmag_error', 'RPmag', 'RPmag_error', 'BPmag', 'BPmag_error', 'BP-RP', 'BP-RP_error',
                'Reddening(BP-RP)', 'Extinction_G', 'Distance', 'Teff', 'logg', '[Fe/H]', 'RUWE',
                'Interferometry_mode', 'Interferometry_guide_star', 'Interferometry_interval',
                'Mag_base_J', 'Mag_base_H', 'Mag_base_K', 'Mag_peak_J', 'Mag_peak_J_error',
                'Mag_peak_H', 'Mag_peak_H_error', 'Mag_peak_K', 'Mag_peak_K_error',
                't0', 't0_error', 'Interferometry_candidate'
                ]
    bool_keys = ['Interferometry_candidate']

    for key in key_list:
        clean_key = clean_key_string(key)
        value = mulens.extras[key].value
        if key in bool_keys and value != None:
            if value:
                context[clean_key] = 'True'
            else:
                context[clean_key] = 'False'
        else:
            context[clean_key] = value

    if mulens.existing_model:
        context['t0_date'] = str(convert_JD_to_UTC(mulens.t0))
    else:
        context['t0_date'] = None

    t2 = datetime.utcnow()
    logger.info('MOP INTERFEROMETRY chk 2: ' + str(t2 - t1))
    utilities.checkpoint()

    # Gather the ReducedData for the neighbouring stars
    # Unpack the QuerySet returned into a more convenient format for display
    if mulens.neighbours:
        rd = mulens.neighbours
        nstars = len(rd.value['Gaia_Source_ID'])
        context['neighbours'] = [{'Gaia_Source_ID': rd.value['Gaia_Source_ID'][i],
                         'Gmag': np.around(rd.value['Gmag'][i],3),
                         'Gmag_error': np.around(rd.value['Gmag_error'][i],3),
                         'BPmag': np.around(rd.value['BPmag'][i],3),
                         'BPmag_error': np.around(rd.value['BPmag_error'][i],3),
                         'RPmag': np.around(rd.value['RPmag'][i],3),
                         'RPmag_error': np.around(rd.value['RPmag_error'][i],3),
                         'BP_RP': np.around(rd.value['BP-RP'][i],3),
                         'BP_RP_error': np.around(rd.value['BP-RP_error'][i],3),
                         'Jmag': np.around(rd.value['Jmag'][i],3),
                         'Hmag': np.around(rd.value['Hmag'][i],3),
                         'Kmag': np.around(rd.value['Kmag'][i],3),
                         'Reddening_BP_RP': np.around(rd.value['Reddening(BP-RP)'][i],3),
                         'Extinction_G': np.around(rd.value['Extinction_G'][i],3),
                         'Distance': rd.value['Distance'][i],
                         'Teff': rd.value['Teff'][i],
                         'logg': rd.value['logg'][i],
                         'Fe_H': rd.value['[Fe/H]'][i],
                         'RUWE': np.around(rd.value['RUWE'][i],3),
                         'Separation': np.around((rd.value['Separation'][i]*3600.0),3)}
                             for i in range(0,nstars,1)]

    else:
        context['neighbours'] = []

    t3 = datetime.utcnow()
    logger.info('MOP INTERFEROMETRY chk 3: ' + str(t3 - t2))
    utilities.checkpoint()

    # Retrieve and unpacked the ReducedData for the nearby GSC stars
    if mulens.gsc_results:
        rd = mulens.gsc_results
        nstars = len(rd.value['GSC2'])
        context['gsc_table'] = [{'GSC2': rd.value['GSC2'][i],
                 'Separation': np.around(rd.value['Separation'][i], 3),
                 'Jmag': np.around(rd.value['Jmag'][i], 3),
                 'Hmag': np.around(rd.value['Hmag'][i], 3),
                 'Ksmag': np.around(rd.value['Ksmag'][i], 3),
                 'W1mag': np.around(rd.value['W1mag'][i], 3),
                 'RA': np.around(rd.value['RA'][i], 6),
                 'Dec': np.around(rd.value['Dec'][i], 6),
                 'plx': np.around(rd.value['plx'][i], 4),
                 'pmRA': np.around(rd.value['pmRA'][i], 3),
                 'pmDE': np.around(rd.value['pmDE'][i], 3),
                 'AOstar': rd.value['AOstar'][i],
                 'FTstar': rd.value['FTstar'][i]}
                                for i in range(0,nstars,1)]
    else:
        context['gsc_table'] = []
    t4 = datetime.utcnow()
    logger.info('MOP INTERFEROMETRY chk 4: ' + str(t4 - t3))
    utilities.checkpoint()

    # Retrieve and unpack the ReducedData for the AOFT table
    maxStrehl = 40
    AOFT_table = []
    FTstar_columns = ['FTstar', 'SC_separation', 'Ksmag', 'SC_Vloss']
    AOstar_columns = ['_FTstrehl', '_Ksmag', '_FT_separation']
    AOstars = []
    if mulens.aoft_table:
        rd = mulens.aoft_table

        # Review the dictionary keys(=column names) to identify the AOstar names
        AOstars = []
        for col in rd.value.keys():
            if '_SCstrehl' in col:
                if len(rd.value[col]) > 0:
                    AOstars.append({'name': col.split('_')[0]})
        context['naostars'] = len(AOstars)

        # Review the datum entries and extract the AO star parameters into columnar data
        suffices = ['_SCstrehl', '_Gmag', '_SC_separation']
        for star in AOstars:
            for suffix in suffices:
                try:
                    star[suffix[1:]] = np.around(rd.value[star['name']+suffix][0], 3)

                    if 'SCstrehl' in suffix:
                        star['scstrehl_colour'] = colour_percent(star[suffix[1:]], norm=maxStrehl)
                except KeyError:
                    star[suffix[1:]] = 0.0

        AOFT_table = []
        nrows = len(rd.value['FTstar'])
        context['nftstars'] = nrows
        col_names = rd.value.keys()
        for i in range(0, nrows, 1):
            row = {'aostars': []}
            for col in FTstar_columns:
                # Extract the parameters relating directly to each FT star
                if col in ['FTstar']:
                    row[col] = rd.value[col][i]
                elif col in ['SC_Vloss', 'SC_separation']:
                    row[col] = np.around(rd.value[col][i], 3)
                    row[str(col).replace('_', '').lower() + '_colour'] = colour_percent(rd.value[col][i])
                elif col in ['Ksmag']:
                    row[col] = np.around(rd.value[col][i], 3)
                    row[str(col).replace('_', '').lower() + '_colour'] = brightness_shader(rd.value[col][i])

            # Extract the parameters stored for each AO star, per FT star
            for aostar in AOstars:
                params = {'name': aostar['name']}
                for suffix in AOstar_columns:
                    col_name = aostar['name']+suffix
                    if col_name in rd.value.keys():
                        if '_FTstrehl' in suffix:
                            params['FTstrehl'] = (np.around(rd.value[col_name][i], 3),
                                                    colour_percent(rd.value[col_name][i], norm=maxStrehl))
                        elif '_Ksmag' in suffix:
                            params['Ksmag'] = (np.around(rd.value[col_name][i], 3),
                                                    brightness_shader(rd.value[col_name][i]))
                        elif '_FT_separation' in suffix:
                            params['FT_separation'] = (np.around(rd.value[col_name][i], 3),
                                                        distance_shader(rd.value[col_name][i]))
                    else:
                        params['FTstrehl'] = (0.0, '#fcfafa')
                        parmas['Ksmag'] = (0.0, '#fcfafa')
                        params['FT_separation'] = (0.0, '#fcfafa')
                row['aostars'].append(params)
            AOFT_table.append(row)

        context['AOFT_table'] = AOFT_table
        context['naostars'] = len(rd.value.keys()) - 4
        context['nftstars'] = nrows
        context['AOstars'] = AOstars
    else:
        context['AOFT_table'] = []
        context['naostars'] = 0
        context['nftstars'] = 0
        context['AOstars'] = []

    t5 = datetime.utcnow()
    logger.info('MOP INTERFEROMETRY chk 5: ' + str(t5 - t4))
    utilities.checkpoint()
    logger.info('MOP INTERFEROMETRY took ' + str(t5 - t1))

    return context

def colour_percent(value, thresholds=None, norm=100):
    """
    Function to perform a look-up of table cell background colour, based on the value of that cell,
    for values in the form of percentages.  Adapted from code by Antoine Merand.
    """
    # Colours: Magenta, red, yellow, green, cyan, blue
    ascii_cols = ['\033[45m', '\033[41m', '\033[43m', '\033[42m', '\033[46m', '\033[44m']
    hex_cols = ['#de47fc', '#f70a0e', '#e8db23', '#23e851', '#23e8db', '#3e96fa']

    if thresholds is None:
        thresholds = norm * np.linspace(0, 1, len(hex_cols) + 1)[1:]

    for i, v in enumerate(thresholds):
        if value <= v:
            return hex_cols[i]

    return hex_cols[-1]

def brightness_shader(Ksmag, kmax=10.5):
    # Default is white
    col = '#fcfafa'

    if Ksmag < kmax - 2:
        #col = '\033[44m'
        col = '#3e96fa'
    elif Ksmag < kmax - 1:
        #col = '\033[46m'
        col = '#23e8db'
    elif Ksmag < kmax:
        #col = '\033[42m'
        col = '#23e851'
    elif Ksmag < kmax + 1:
        #col = '\033[43m'
        col = '#e8db23'
    else:
        #col = '\033[41m'
        col = '#f70a0e'

    return col

def distance_shader(dist, distmax=30):
    """
    Function returns the table cell background colour based on the angular separation in arcsec.
    """
    # Default is purple
    col = '#de47fc'

    if dist == 0:
        #cold.append('\033[44m')
        col = '#3e96fa'
    elif dist < distmax / 4:
        #cold.append('\033[46m')
        col = '#23e8db'
    elif dist < distmax / 2:
        #cold.append('\033[42m')
        col = '#23e851'
    elif dist < 3 * distmax / 4:
        #cold.append('\033[43m')
        col = '#e8db23'
    elif dist <= distmax:
        #cold.append('\033[41m')
        col = '#f70a0e'

    return col

def convert_JD_to_UTC(jd):
    try:
        t = Time(jd, format='jd')
        t = t.utc
        t.format = 'iso'
        t.out_subfmt = 'date'
        ts = t.value
    except TypeError:
        ts = None
    except ValueError:
        ts = None
    return ts

@register.inclusion_tag('tom_dataproducts/partials/gaia_neighbours_data.html')
def gaia_neighbours_data(mulens):
    context = {}

    # Gather the ReducedData for the neighbouring stars
    # Unpack the QuerySet returned into a more convenient format for display
    neighbours = []
    if mulens.neighbours:
        rd = mulens.neighbours
        nstars = len(rd.value['Gaia_Source_ID'])
        neighbours = [
                    {'Gaia_Source_ID': rd.value['Gaia_Source_ID'][i],
                     'Gmag': np.around(rd.value['Gmag'][i],3),
                     'Gmag_error': np.around(rd.value['Gmag_error'][i],3),
                     'BPmag': np.around(rd.value['BPmag'][i],3),
                     'BPmag_error': np.around(rd.value['BPmag_error'][i],3),
                     'RPmag': np.around(rd.value['RPmag'][i],3),
                     'RPmag_error': np.around(rd.value['RPmag_error'][i],3),
                     'BP_RP': np.around(rd.value['BP-RP'][i],3),
                     'BP_RP_error': np.around(rd.value['BP-RP_error'][i],3),
                     'Jmag': np.around(rd.value['Jmag'][i],3),
                     'Hmag': np.around(rd.value['Hmag'][i],3),
                     'Kmag': np.around(rd.value['Kmag'][i],3),
                     'Reddening_BP_RP': np.around(rd.value['Reddening(BP-RP)'][i],3),
                     'Extinction_G': np.around(rd.value['Extinction_G'][i],3),
                     'Distance': rd.value['Distance'][i],
                     'Teff': rd.value['Teff'][i],
                     'logg': rd.value['logg'][i],
                     'Fe_H': rd.value['[Fe/H]'][i],
                     'RUWE': np.around(rd.value['RUWE'][i],3),
                     'Separation': np.around((rd.value['Separation'][i]*3600.0),3)}
                    for i in range(0, nstars, 1) ]

    context['neighbours'] = neighbours

    print('END ', datetime.utcnow())

    return context
@register.inclusion_tag('tom_targets/partials/current_timestamp.html')
def current_timestamp():

    context = {}
    utc_now = datetime.utcnow()
    context['utc_now'] = utc_now.strftime("%Y-%m-%d %H:%M:%S")
    context['jd_now'] = str(np.around(Time(utc_now).jd,3))

    return context


@register.inclusion_tag('tom_targets/partials/mulens_target_data.html')
def mulens_target_data(mulens, request):
    """
    Displays the data of a target.
    """

    t1 = datetime.utcnow()
    utilities.checkpoint()
    target_data = {
        'target': mulens.target,
        'extras': mulens.extras,
        'request': request
    }
    utilities.checkpoint()
    t2 = datetime.utcnow()
    logger.info('MULENS TARGET DATA time taken: ' + str(t2 - t1))

    return target_data

@register.inclusion_tag('tom_targets/partials/target_class_form.html', takes_context=True)
def classification_form(context):
    """Embedded form to enable user to change a Target's classification extra_field parameters"""

    target = context['target']
    try:
        request = context['request']
    except KeyError:
        request = None

    class_form = TargetClassificationForm()
    default_classes = [x[0] for x in class_form.fields['classification'].choices]
    default_categories = [x[0] for x in class_form.fields['category'].choices]

    # If there is no request associated with the call to this function, try and extract the classification
    # parameters from the target extras
    class_data = {'classification': None, 'text_class': None, 'category': None, 'text_category': None}
    if not request:
        for key in ['Classification', 'Category']:
            if key in context['extras'].keys():
                class_data[key.lower()] = context['extras'][key]

    # If the user has set values for any of the form fields, then populate the class_form with those values
    else:
        if any(request.GET.get(x) for x in class_data.keys()):
            class_data['classification'] = get_request_param('classification', 'Classification', request, context['extras'])
            class_data['text_class'] = request.GET.get('text_class', '')
            class_data['category'] = get_request_param('category', 'Category', request, context['extras'])
            class_data['text_category'] = request.GET.get('text_category', '')
        else:
            for key in ['Classification', 'Category']:
                if key in context['extras'].keys():
                    class_data[key.lower()] = context['extras'][key]

    class_form = TargetClassificationForm(class_data)

    if class_form.is_valid() and request:
        extras = {
            'Classification': class_data['classification'],
            'Category': class_data['category']
        }

        # If the user has entered their own text into the text fields, this takes priority, otherwise
        # take the option set from the select menu
        if len(class_form.cleaned_data['text_class']) > 0:
            extras['Classification'] = class_form.cleaned_data['text_class']
        elif len(class_form.cleaned_data['classification']) > 0:
            extras['Classification'] = class_form.cleaned_data['classification']

        if len(class_form.cleaned_data['text_category']) > 0:
            extras['Category'] = class_form.cleaned_data['text_category']
        elif len(class_form.cleaned_data['category']) > 0:
            extras['Category'] = class_form.cleaned_data['category']

        # Save the updated extra_field parameters
        target.save(extras=extras)

        # Return a refreshed, empty form:
        class_form = TargetClassificationForm()

    result = {
        'form': class_form,
        'target': target,
        'categories': default_categories,
        'classifications': default_classes,
    }
    if 'Category' in target.extra_fields.keys():
        result['current_category'] = target.extra_fields['Category']
    if 'Classification' in target.extra_fields.keys():
        result['current_classification'] = target.extra_fields['Classification']

    return result

def get_request_param(request_key, extras_key, request, extras):
    """Fetch a parameter value from the request object, but fall back to the already-set extra-parameter
    value if the result is an empty string"""

    request_value = request.GET.get(request_key, '')
    if len(request_value) == 0:
        request_value = extras[extras_key]

    return request_value
