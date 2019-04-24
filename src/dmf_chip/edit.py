from __future__ import absolute_import
import functools as ft
import warnings

from logging_helpers import _L
from lxml.etree import QName, Element
import lxml.etree
import networkx as nx
import numpy as np
import pandas as pd

from .core import ureg
from .load import draw, load
from six.moves import zip

__all__ = ['detect_neighbours', 'draw_with_segment_rays',
           'write_connections_layer']


DEFAULT_DISTANCE_THRESHOLD = 0.175 * ureg.mm


def detect_neighbours(chip_info,
                      distance_threshold=DEFAULT_DISTANCE_THRESHOLD):
    segments = get_segment_rays(chip_info, magnitude=distance_threshold)
    return get_all_intersections(segments)


def draw_with_segment_rays(chip_info,
                           distance_threshold=DEFAULT_DISTANCE_THRESHOLD,
                           axis=None):
    import matplotlib.pyplot as plt

    if axis is None:
        fig, axis = plt.subplots(figsize=(50, 50))
    result = draw(chip_info, ax=axis)
    # result = draw(chip_info)
    axis = result['axis']

    for p in result['patches'].values():
        p.set_alpha(.3)

    light_green = '#90cd97'
    dark_green = '#059748'

    df_intersections = detect_neighbours(chip_info, distance_threshold=.175 *
                                         ureg.mm)
    for idx_i, segment_i in df_intersections.iterrows():
        axis.arrow(segment_i['x_mid'], segment_i['y_mid'],
                segment_i['x_normal'], segment_i['y_normal'],
                width=.25,
                edgecolor=dark_green, facecolor=light_green)


def get_all_intersections(df_rays):
    '''
    Parameters
    ----------
    segment_rays : pandas.DataFrame
        See return type of :func:`get_segment_rays()`.
    '''
    intersections = []

    for i, ((id_i, vertex_i), segment_i) in enumerate(df_rays.iterrows()):
        p = segment_i[['x_mid', 'y_mid']].values
        r = segment_i[['x_normal', 'y_normal']].values

        df_intersections_i = get_intersections(df_rays, p, r)

        # Do not include self electrode in consideration for neighbours.
        self_mask = df_intersections_i.index.get_level_values('id') == id_i
        df_intersections_i = df_intersections_i.loc[~self_mask]
        if df_intersections_i.shape[0]:
            intersections.append(((id_i, vertex_i), df_intersections_i))
    if not intersections:
        return pd.DataFrame()

    index, values = list(zip(*intersections))
    df_result = pd.concat(values, keys=index)
    df_result.index.names = ['id', 'vertex_i',
                             'id_neighbour', 'vertex_i_neighbour']
    return df_result


def get_intersections(df_rays, p, r):
    # See: https://stackoverflow.com/a/565282/345236
    q = df_rays[['x1', 'y1']].values
    s = df_rays[['x2', 'y2']].values - q

    r_x_s = np.cross(r, s)
    r_x_s[r_x_s == 0] = np.NaN
    t = np.cross((q - p), s) / r_x_s
    u = np.cross((q - p), r) / r_x_s

    df_tu = pd.DataFrame(np.column_stack([t, u]), columns=list('tu'),
                         index=df_rays.index)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        df_i = df_rays.join(df_tu).loc[(r_x_s != 0)
                                       & (t >= 0) & (t <= 1) & (u >= 0)
                                       & (u <= 1)]
    intersect_points = p + df_i.t.values[:, None] * r
    return df_i.join(pd.DataFrame(intersect_points, columns=['x_intersect',
                                                             'y_intersect'],
                                  index=df_i.index)).drop(['t', 'u'], axis=1)


def _electrode_segment_rays(electrode, magnitude):
    '''Compute ray cast "outwards" for each line segment of electrode shape.

    Parameters
    ----------
    electrode : dict
        See ``electrodes`` item in :func:`dmf_chip.load()`.
    magnitude : float
        Magnitude of ray vectors (in pixels).

    Returns
    -------
    pandas.DataFrame
        Each row corresponds to a ray vector cast from the respective line
        segment in the electrode shape, with the following columns::

         - ``x1``, ``y1``: start point of line segment
         - ``x2``, ``y2``: end point of line segment
         - ``x_mid``, ``y_mid``: mid point of line segment
         - ``length``: Cartesian length of line segment
         - ``x_normal``, ``y_normal``: end point of cast ray
    '''
    points = np.array(electrode['points'])
    if electrode['direction'] == 'counter-clockwise':
        points = points[::-1]
    # Vector direction/magnitude for each segment (relative to origin).
    v = .5 * (points[1:] - points[:-1])
    # Mid-point of segment.
    x_mid, y_mid = .5 * (points[1:] + points[:-1]).T
    length = np.sqrt((v ** 2)).sum(axis=1)
    v_scaled = magnitude * (v / length[:, None])
    x_normal = -v_scaled[:, 1]
    y_normal = v_scaled[:, 0]

    x1, y1 = points[:-1].T
    x2, y2 = points[1:].T
    result = pd.DataFrame(np.column_stack((x1, y1, x2, y2, x_mid, y_mid, length,
                                         x_normal, y_normal)),
                        columns=['x1', 'y1', 'x2', 'y2', 'x_mid', 'y_mid',
                                 'length', 'x_normal', 'y_normal'])
    return result


def get_segment_rays(chip_info, magnitude=DEFAULT_DISTANCE_THRESHOLD):
    magnitude_px = (magnitude * chip_info['__metadata__']['ppi'] *
                    ureg.ppi).to('pixel').magnitude
    df_rays = pd.concat([_electrode_segment_rays(e_i, magnitude_px)
                         for e_i in chip_info['electrodes']],
                        keys=[e['id'] for e in chip_info['electrodes']])
    df_rays.index.names = 'id', 'vertex_i'
    return df_rays


def write_connections_layer(chip_file,
                            distance_threshold=DEFAULT_DISTANCE_THRESHOLD):
    chip_info = load(chip_file)
    df_intersections = detect_neighbours(chip_info,
                                         distance_threshold=distance_threshold)

    doc = lxml.etree.parse(chip_file)
    root = doc.getroot()
    nsmap = {k: v for k, v in root.nsmap.items() if k}
    _xpath = ft.partial(root.xpath, namespaces=nsmap)

    device_layer = _xpath('//svg:g[@inkscape:label="Device"]')[0]
    connections_layers = _xpath('//svg:g[@inkscape:label="Connections"]')

    # Remove existing neighbouring electrode connections layer(s) (if any).
    for layer in connections_layers:
        root.remove(layer)

    # Determine and use first unused layer label number.
    layer_ids = set(_xpath('//svg:g[@inkscape:label and @inkscape:groupmode='
                           '"layer"]/@id'))
    i = 1
    while True:
        layer_id = 'layer%d' % i
        if layer_id not in layer_ids:
            break
        i += 1

    connections_layer = Element(QName(nsmap['svg'], 'g'),
                                attrib={QName(nsmap['inkscape'], 'label'):
                                        'Connections',
                                        QName(nsmap['inkscape'], 'groupmode'):
                                        'layer', 'id': layer_id})

    # Construct undirected graph from detected intersections.
    edges = df_intersections.reset_index()[['id',
                                            'id_neighbour']].values.tolist()
    graph = nx.Graph(edges)

    # Create one `<svg:path>` per electrode.
    path_elements = []

    centers = pd.Series((e['pole_of_accessibility']
                         for e in chip_info['electrodes']),
                        index=[e['id'] for e in chip_info['electrodes']])

    for a, b, in graph.edges:
        a_point, b_point = centers[[a, b]]
        path_d = 'M %.2f,%.2f L %.2f,%.2f' % (a_point['x'], a_point['y'],
                                            b_point['x'], b_point['y'])
        path_elem = Element(QName(nsmap['svg'], 'path'),
                            attrib={'id': layer_id,
                                    'style': 'stroke:#000000;stroke-width:0.1',
                                    'd': path_d})
        path_elements.append(path_elem)

    connections_layer.extend(path_elements)
    device_layer.addnext(connections_layer)
    return doc


def _get_or_create(parent, name, attrib=None):
    '''Get element specified by qualified tag name or create it.

    Parameters
    ----------
    parent : lxml.etree element
        Parent element.
    name : str
        Name in form ``"<namespace alias>:<tagname>"``, e.g.,
        ``"dmf:ChipDesign"``.  If :data:`parent` does not contain a child
        matching the specified tag name and corresponding attributes, create a
        new element.
    attrib : dict, optional
        Element attributes to match (or set, if creating new element).

    Returns
    -------
    lxml.etree.Element
        Matching child element (if available) or created element.

    Examples
    --------

    Get ``<dmf:ChipDesign>`` element or create it if it does not exist:

    >>>> from dmf_chip.edit import _get_or_create
    >>>>
    >>>> # Load xml document define `_xpath` alias...
    >>>>
    >>>> metadata = _xpath('/svg:svg/svg:metadata')[0]
    >>>> chip_design = _get_or_create(metadata, 'dmf:ChipDesign')
    '''
    docroot = parent.getroottree().getroot()
    nsmap = {k: v for k, v in docroot.nsmap.items() if k}
    ns, tagname = name.split(':')
    qname = QName(nsmap[ns], tagname)
    # Short-hand to xpath using namespaces referenced in file.
    _xpath = ft.wraps(parent.xpath)(ft.partial(parent.xpath, namespaces=nsmap))
    xquery = './%s:%s' % (ns, tagname)
    if attrib is not None:
        attrib_str = ''.join('[@%s="%s"]' % (k, v) for k, v in attrib.items())
    else:
        attrib_str = ''
    xquery += attrib_str

    if not _xpath(xquery):
        element = Element(qname, attrib=attrib)
        parent.append(element)
        _L().info('Add new element: `%s:%s%s`', ns, tagname, attrib_str)
    else:
        element = _xpath(xquery)[0]
        _L().info('found element: `%s:%s%s`', ns, tagname, attrib_str)
    return element


def write_test_route(chip_file, tour_ids, id_):
    '''Write test route to SVG metadata.

    Parameters
    ----------
    chip_file : str
        Path to chip design file.
    tour_ids : list[str]
        Ordered list of electrode ids defining tour waypoints.
    id_ : str
        Test route id.

    Returns
    -------
    lxml.etree document
        In-memory document with test route element added.
    '''
    doc = lxml.etree.parse(chip_file)
    root = doc.getroot()

    if 'dmf' not in root.nsmap:
        root.nsmap['dmf'] = \
            "https://github.com/sci-bots/dmf-chip-spec/releases/tag/v0.1"

    NSMAP = {k: v for k, v in root.nsmap.items() if k}
    # Short-hand to xpath using namespaces referenced in file.
    _xpath = ft.wraps(root.xpath)(ft.partial(root.xpath, namespaces=NSMAP))

    metadata = _xpath('/svg:svg/svg:metadata')[0]

    chip_design = _get_or_create(metadata, 'dmf:ChipDesign')
    test_routes = _get_or_create(chip_design, 'dmf:TestRoutes')

    if test_routes.xpath('./dmf:TestRoute[@id="%s"]' % id_, namespaces=NSMAP):
        raise NameError('Test route already exists with id: `%s`', id_)

    test_route = _get_or_create(test_routes, 'dmf:TestRoute',
                                attrib={'id': id_, 'version': '0.1.0'})
    for id_i in tour_ids:
        element_i = Element(QName(NSMAP['dmf'], 'Waypoint'))
        element_i.text = str(id_i)
        test_route.append(element_i)
    _L().info('Added %d waypoints.', len(tour_ids))
    return doc
