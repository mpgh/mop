from django import template
from django.conf import settings
from django.core.paginator import Paginator
from datetime import datetime

from plotly import offline
import plotly.graph_objs as go

from tom_dataproducts.models import DataProduct, ReducedDatum
from tom_dataproducts.processors.data_serializers import SpectrumSerializer
from mop.toolbox import utilities
from astropy.time import Time
import numpy as np

register = template.Library()



@register.inclusion_tag('tom_dataproducts/partials/photometry_for_target.html')
def mop_photometry(target):
    """
    Renders a photometric plot for a target.
    This templatetag requires all ``ReducedDatum`` objects with a data_type of ``photometry`` to be structured with the
    following keys in the JSON representation: magnitude, error, filter
    """
    photometry_data = {}
    for datum in ReducedDatum.objects.filter(target=target, data_type=settings.DATA_PRODUCT_TYPES['photometry'][0]):
        values = datum.value
        try:
           
                photometry_data.setdefault(values['filter'], {})
                photometry_data[values['filter']].setdefault('time', []).append(Time(datum.timestamp).jd-2450000)
                photometry_data[values['filter']].setdefault('magnitude', []).append(values.get('magnitude'))
                photometry_data[values['filter']].setdefault('error', []).append(values.get('error'))

        except:
                
                photometry_data.setdefault(values['filter'], {})
                photometry_data[values['filter']].setdefault('time', []).append(Time(datum.timestamp).jd-2450000)
                photometry_data[values['filter']].setdefault('magnitude', []).append(values.get('magnitude'))
                photometry_data[values['filter']].setdefault('error', []).append(values.get('error'))
    plot_data = [
        go.Scatter(
            x=filter_values['time'],
            y=filter_values['magnitude'], mode='markers',
            name=filter_name,
            error_y=dict(
                type='data',
                array=filter_values['error'],
                visible=True
            )
        ) for filter_name, filter_values in photometry_data.items()]
   


    layout = go.Layout(
        yaxis=dict(autorange='reversed'),
        height=600,
        width=700,
             
    )
   
    fig = go.Figure(data=plot_data, layout=layout)
    current_time =  Time.now().jd-2450000
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
    try:
       for datum in ReducedDatum.objects.filter(target=target, data_type=settings.DATA_PRODUCT_TYPES['lc_model'][0]).order_by('-id')[:1]:

            values =datum.value
            time = np.array(values['lc_model_time'])-2450000
            mag = np.array(values['lc_model_magnitude'])
            fig.add_trace(go.Scatter(x=time, y=mag,
                    mode='lines',
                    name='Model',
                    opacity=0.5,
                    line = dict(color='rgb(128,128,128)',
                                width=5,
                                ),
                     ))
    except:

       pass
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
             xaxis_title="HJD-2450000",
             yaxis_title="Mag",
    )

    return {
        'target': target,
        'plot': offline.plot(fig, output_type='div', show_link=False)
    }


@register.inclusion_tag('tom_dataproducts/partials/interferometry_data.html')
def interferometry_data(target):
    context = {}

    # Gather extra_param key values for target's Gaia catalog data
    key_list = ['Gaia_Source_ID',
                'Gmag', 'Gmag_error', 'RPmag', 'RPmag_error', 'BPmag', 'BPmag_error', 'BP-RP', 'BP-RP_error',
                'Reddening(BP-RP)', 'Extinction_G', 'Distance', 'Teff', 'logg', '[Fe/H]', 'RUWE',
                'Interferometry_mode', 'Interferometry_guide_star', 'Mag_peak_J', 'Mag_peak_H', 'Mag_peak_K']
    for key in key_list:
        clean_key = key.replace('[','').replace(']','').replace('/','_').replace('(','_').replace(')','_').replace('-','_')
        value = utilities.fetch_extra_param(target, key)
        if key in ['Gmag', 'Gmag_error', 'RPmag', 'RPmag_error', 'BPmag', 'BPmag_error', 'BP-RP', 'BP-RP_error',
                'Reddening(BP-RP)', 'Extinction_G', 'Mag_peak_J', 'Mag_peak_H', 'Mag_peak_K'] and value != None:
            context[clean_key] = np.around(value,3)
        else:
            context[clean_key] = value

    # Gather the ReducedData for the neighbouring stars
    # Unpack the QuerySet returned into a more convenient format for display
    qs = ReducedDatum.objects.filter(target=target, data_type='tabular', source_name='Interferometry_predictor')
    neighbours = []
    if len(qs) > 0:
        rd = qs[0]
        nstars = len(rd.value['Gaia_Source_ID'])
        for i in range(0,nstars,1):
            neighbours.append(
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
            )
    context['neighbours'] = neighbours

    return context