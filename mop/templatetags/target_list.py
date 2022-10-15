from django import template
from django.conf import settings
import dash
import dash_html_components as html
import dash_table
from django_plotly_dash import DjangoDash
import pandas as pd

register = template.Library()

@register.inclusion_tag('oss/partials/target_table.html')
def target_table(targets):
    """Produces a Dash interactive table of targets"""

    table_columns = [dict(name='Name', id='Name', type='text', presentation='markdown'),
                     dict(name='Type', id='Type', type='text', presentation='markdown'),
                     dict(name='Observations', id='Observations', type='text', presentation='markdown'),
                     dict(name='Saved Data', id='Saved_Data')]

    table_data = []
    for target in targets:
        table_data.append( dict(Name=target.names.join(', '),
                                Type=target.get_type_display,
                                Observations=target.observationrecord_set.count,
                                Saved_Data=target.dataproduct_set.count) )

    app = build_target_table_app(table_columns, table_data)

    return {'request': app}

def build_target_table_app(table_columns, table_data):

    app = DjangoDash('TargetsTable')

    app.layout = html.Div( dash_table.DataTable(
                    id='TargetsTable',
                    columns=table_columns,
                    data=table_data,
                    sort_action="native",
                    filter_action="native",
                    style_table={'height': '600px', 'overflowY': 'auto'},
                    style_cell={'fontSize':18, 'font-family':'sans-serif'},
                    style_cell_conditional=[],
                    ) )

    return app
