from argparse import Namespace
from collections import OrderedDict
from copy import deepcopy
from hashlib import sha256
import functools as ft
import json
import logging
import lxml
import re

from shapely.algorithms.polylabel import polylabel
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import numpy as np
import pandas as pd
import semantic_version as sv

from .core import ureg
from .svg import shape_points

__all__ = ['draw', 'load', 'to_unit', 'as_obj']

logger = logging.getLogger(__name__)


def load(chip_file):
    '''
    Parameters
    ----------
    chip_file : str
        Path to chip design file.

    Returns
    -------
    dict


    .. versionchanged:: 0.3.0
        Read design ID and test routes from ``<dmf:ChipDesign>`` tag.
        See https://github.com/sci-bots/dmf-chip/issues/1 for more information.
    '''
    info = {'__metadata__': {}}

    with open(chip_file, 'rb') as input_:
        info['__metadata__']['sha256'] = sha256(input_.read()).hexdigest()

    root = lxml.etree.parse(chip_file).getroot()
    NSMAP = {k: v for k, v in root.nsmap.items() if k}
    # Short-hand to xpath using namespaces referenced in file.
    _xpath = ft.wraps(root.xpath)(ft.partial(root.xpath, namespaces=NSMAP))

    try:
        inkscape_version = _xpath('/svg:svg/@inkscape:version')[0]
        semantic_version = sv.Version(re.split('\s+', inkscape_version)[0],
                                      partial=True)
        if semantic_version >= sv.Version('0.92.0'):
            logger.info('Detected Inkscape 0.92+; using 96 pixels per inch')
            info['__metadata__']['ppi'] = 96
        else:
            logger.info('Detected Inkscape <0.92; using 90 pixels per inch')
            info['__metadata__']['ppi'] = 90
    except Exception:
        logger.exception('')

    if 'ppi' not in info['__metadata__']:
        logger.info('Using 96 pixels per inch')
        info['__metadata__']['ppi'] = 96

    if all(k in root.attrib for k in ('width', 'height')):
        ppi = info['__metadata__']['ppi'] * ureg.pixel / ureg.inch
        shape = {k: ureg.parse_expression(root.attrib[k])
                 for k in ('width', 'height')}
        for k, v in shape.items():
            if isinstance(v, ureg.Quantity):
                shape[k] = (v * ppi).to('pixel').magnitude
        info['__metadata__'].update(shape)

    # Read design-id from `<dmf:ChipDesign><dmf:DesignID>` tag.
    ids = _xpath('//dmf:ChipDesign/dmf:DesignId')
    if ids:
        info['__metadata__']['design-id'] = ids[0].text

    # Read test routes.
    test_route_elements = \
        _xpath('//dmf:ChipDesign/dmf:TestRoutes/dmf:TestRoute[@id!=""]')
    test_routes = []
    for route in test_route_elements:
        xpath_ = ft.wraps(route.xpath)(ft.partial(route.xpath,
                                                  namespaces=NSMAP))
        route_ = dict(route.attrib.items())
        route_['waypoints'] = [w.text for w in xpath_('dmf:Waypoint')]
        test_routes.append(route_)
    info['__metadata__']['test-routes'] = test_routes

    # Extract electrode information.
    device_layer = _xpath('//svg:g[@inkscape:label="Device"]')[0]
    _device_xpath = ft.partial(device_layer.xpath, namespaces=NSMAP)
    electrode_elements = _device_xpath('.//svg:path | .//svg:polygon')

    electrodes = []

    for p in electrode_elements:
        electrode_info =  {'points': shape_points(p),
                           '__tag__': p.tag,
                           '__sourceline__': p.sourceline}
        if 'data-channels' in p.attrib and p.attrib['data-channels']:
            electrode_info['channels'] = \
                map(int, re.split(r'\s*,\s*', p.attrib['data-channels']))
        electrode_info.update(p.attrib)
        electrodes.append(electrode_info)

    info['electrodes'] = electrodes

    # Extract electrode connections information.
    #
    # Create `Polygon` instances to allow fast collision detection between
    # end-points of each connection line and electrodes (if any).
    electrode_polygons = {}
    for electrode_i in info['electrodes']:
        try:
            polygon_i = Polygon(electrode_i['points'])
            # Compute electrode ["pole of accessibility"][1]; similar to
            # centroid, but _guaranteed to be within shape_.
            #
            # [1]: https://github.com/mapbox/polylabel
            pole_i = polylabel(polygon_i)
            electrode_i['pole_of_accessibility'] = {'x': pole_i.x,
                                                    'y': pole_i.y}
            electrode_polygons[electrode_i['id']] = polygon_i
        except:
            logger.exception('Error: `%s`' % electrode_i)

    def find_shape(x, y):
        point = Point(x, y)
        for id_, polygon_i in electrode_polygons.items():
            if polygon_i.contains(point):
                return id_
        else:
            return None

    connections_layer = _xpath('//svg:g[@inkscape:label="Connections"]')[0]
    _connections_xpath = ft.partial(connections_layer.xpath, namespaces=NSMAP)
    connection_elements = _connections_xpath('.//svg:line | .//svg:path')
    connections = []

    for element_i in connection_elements:
        info_i = dict(element_i.attrib)
        info_i['__tag__'] = element_i.tag
        info_i['__sourceline__'] = element_i.sourceline
        info_i['points'] = shape_points(element_i)
        for point, name in ((info_i['points'][0], 'source'),
                            (info_i['points'][-1], 'target')):
            id_ = find_shape(*point)
            if id_ is None:
                continue
            info_i[name] = {'id': id_,
                            'x': point[0],
                            'y': point[1]}
        connections.append(info_i)
    info['connections'] = connections

    # Compute electrode areas using vectorized form of [Shoelace formula][1].
    #
    # [1]: http://en.wikipedia.org/wiki/Shoelace_formula
    electrodes_by_id = {e['id']: e for e in info['electrodes']}
    df = pd.concat([pd.DataFrame(e['points'], columns=['x', 'y'])
                    for e in info['electrodes']],
                   keys=[e['id'] for e in info['electrodes']])
    df.index.names = 'id', 'vertex_i'
    # "rank" denotes the number of vertices in the respective electrode SVG shape.
    df['rank'] = df.groupby(level='id')['x'] .transform('count').astype(int)
    area_a = df.x.copy()
    area_b = df.y.copy()
    vertex = df.index.get_level_values('vertex_i')

    area_a[vertex == df['rank'] - 1] *= df.loc[vertex == 0, 'y'].values
    area_a[vertex < df['rank'] - 1] *= df.loc[vertex > 0, 'y'].values

    area_b[vertex == df['rank'] - 1] *= df.loc[vertex == 0, 'x'].values
    area_b[vertex < df['rank'] - 1] *= df.loc[vertex > 0, 'x'].values

    df['area_a'] = area_a
    df['area_b'] = area_b

    area_components = df.groupby(level='id')[['area_a', 'area_b']].sum()
    shape_areas = .5 * (area_components['area_b'] - area_components['area_a'])

    for id_i, area_i in shape_areas.iteritems():
        electrodes_by_id[id_i]['area'] = abs(area_i)
        electrodes_by_id[id_i]['direction'] = ('clockwise' if area_i >= 0
                                               else 'counter-clockwise')
    return info


def draw(chip_info, ax=None, labels=True, groupby='id', unit='pixel'):
    '''
    Draw the specified device, along with rays casted normal to the electrode
    line segments that intersect with a line segment of a neighbouring
    electrode.


    Parameters
    ----------
    svg_source : str or file-like or pandas.DataFrame
        File path, URI, or file-like object for SVG device file.

        If specified as ``pandas.DataFrame``, assume argument is in format
        returned by :func:`svg_model.svg_shapes_to_df`.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        Axis to draw on.
    labels : bool, optional
        Draw channel labels (default: ``True``).

    Returns
    -------
    `dict`
        Result `dict` includes::
        - ``axis``: axis to which the device was drawn
          (`matplotlib.axes._subplots.AxesSubplot`)
        - ``df_shapes``: table of electrode shape vertices (`pandas.DataFrame`)
        - ``electrode_channels``: mapping from channel number to electrode ID
          (`pandas.Series`).
        - ``channel_patches``: mapping from channel number to corresponding
          `matplotlib` electrode `Patch`.  May be used, e.g., to set color and
          alpha (`pandas.Series`).
    '''
    from matplotlib.patches import Polygon as Polygon_
    import matplotlib.pyplot as plt

    if unit != 'pixel':
        chip_info = to_unit(chip_info, unit)
    df_paths = pd.concat([pd.DataFrame(p['points'], columns=list('xy'))
                          for p in chip_info['electrodes']],
                         keys=[p[groupby] for p in chip_info['electrodes']])

    patches = OrderedDict(sorted([(id_,
                                   Polygon_(df_shape_i[['x', 'y']].values,
                                            closed=False, label=id_))
                                  for id_, df_shape_i in
                                  df_paths.groupby(level=0)]))

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.set_aspect(True)

    # For colors, see: https://gist.github.com/cfobel/fd939073cf13a309d7a9
    for patch_i in patches.values():
        # Light blue
        patch_i.set_facecolor('#88bde6')
        # Medium grey
        patch_i.set_edgecolor('#4d4d4d')
        ax.add_patch(patch_i)

    ax.set_xlim(df_paths.x.min(), df_paths.x.max())
    ax.set_ylim(df_paths.y.max(), df_paths.y.min())

    if labels:
        for electrode_i in chip_info['electrodes']:
            center_i = electrode_i['pole_of_accessibility']
            ax.text(center_i['x'], center_i['y'], electrode_i['id'],
                    horizontalalignment='center', verticalalignment='center',
                    color='white', fontsize=10, bbox={'facecolor': 'black',
                                                      'alpha': 0.2, 'pad': 5})

    ppi = chip_info['__metadata__']['ppi'] * ureg.pixel / ureg.inch
    width = chip_info['__metadata__']['width']
    height = chip_info['__metadata__']['height']
    if unit != 'pixel':
        width = (width * ureg.pixel / ppi).to(unit).magnitude
        height = (height * ureg.pixel / ppi).to(unit).magnitude
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)

    return {'axis': ax, 'patches': patches}


def to_unit(chip_info, unit_name):
    output = deepcopy(chip_info)
    ppi = chip_info['__metadata__']['ppi'] * ureg.pixel / ureg.inch
    unit = ureg.parse_units(unit_name)

    for e in output['electrodes']:
        e['points'] = (np.array(e['points']) * ureg.pixel /
                       ppi).to(unit).magnitude.tolist()
        e['area'] = (np.array(e['area']) * (ureg.pixel ** 2) /
                     (ppi ** 2)).to(unit ** 2).magnitude
        for dim in 'xy':
            e['pole_of_accessibility'][dim] = (e['pole_of_accessibility'][dim]
                                               * ureg.pixel /
                                               ppi).to(unit).magnitude
    for c in output['connections']:
        c['points'] = (np.array(e['points']) * ureg.pixel /
                       ppi).to(unit).magnitude.tolist()
        for k in ('source', 'target'):
            if k in c:
                for dim in 'x', 'y':
                    c[k][dim] = (c[k][dim] * ureg.pixel /
                                 ppi).to(unit).magnitude
    output['__metadata__']['unit'] = unit_name
    return output


def as_obj(chip_info):
    '''Convert nested chip info dictionary into nested object with attributes.

    This is helpful, for example, to enable tab-completion in an interactive
    shell or Jupyter notebook session.

    Parameters
    ----------
    chip_info : dict
        Nested dictionary containing chip information.

        See return type of :func:`load()`.

    Returns
    -------
    argparse.Namespace
        Nested namespace object equivalent to specified nested chip info
        `dict`.


    .. versionadded:: 0.2.0
    '''
    return json.loads(json.dumps(chip_info),
                      object_hook=lambda d: Namespace(**d))
