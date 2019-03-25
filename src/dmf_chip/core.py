# coding: utf-8
from __future__ import absolute_import, unicode_literals, print_function
import itertools as it
import warnings

import numpy as np
import pandas as pd
import pint
import svg_model

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

__all__ = ['chip_info', 'get_all_intersections', 'get_channel_neighbours',
           'get_intersections', 'get_segments', '__version__', 'ureg']

ureg = pint.UnitRegistry()

DEFAULT_DISTANCE_THRESHOLD = 0.3 * ureg.mm

# Only read interpret SVG paths and polygons from `Device` layer as electrodes.
ELECTRODES_XPATH = (r'//svg:g[@inkscape:label="Device"]//svg:path | '
                    r'//svg:g[@inkscape:label="Device"]//svg:polygon')


def _extract_electrode_channels(df_shapes):
    if 'data-channels' not in df_shapes:
        # No electrode has an associated actuation channel.
        return pd.Series(name='channel', index=pd.Series(name='id'))
    else:
        electrode_channels = (df_shapes.drop_duplicates(['id',
                                                         'data-channels'])
                              .set_index('id')['data-channels'].dropna().map(int))
        electrode_channels.name = 'channel'
        return electrode_channels


def _resolve_source(svg_source, **kwargs):
    '''
    Parameters
    ----------
    svg_source : str or file-like or pandas.DataFrame
        File path, URI, or file-like object for SVG device file.

        If specified as ``pandas.DataFrame``, assume argument is in format
        returned by :func:`svg_model.svg_shapes_to_df`.

    Returns
    -------
    pandas.DataFrame
        See return type of :func:`svg_model.svg_shapes_to_df()`.
    '''
    if not isinstance(svg_source, pd.DataFrame):
        return svg_model.svg_shapes_to_df(svg_source, **kwargs)
    else:
        return svg_source


def get_segments(svg_source, distance_threshold=DEFAULT_DISTANCE_THRESHOLD,
                 **kwargs):
    '''
    Parameters
    ----------
    svg_source : str or file-like or pandas.DataFrame
        File path, URI, or file-like object for SVG device file.

        If specified as ``pandas.DataFrame``, assume argument is in format
        returned by :func:`svg_model.svg_shapes_to_df`.
    distance_threshold : pint.quantity.Quantity
        Maximum gap between electrodes to still be considered neighbours.
    '''
    df_shapes = _resolve_source(svg_source, **kwargs)

    # Calculate distance in pixels assuming 96 pixels per inch (PPI).
    distance_threshold_px = (distance_threshold * 96 *
                             ureg.pixels_per_inch).to('pixels')

    df_segments = (df_shapes.groupby('id').apply(lambda x: x.iloc[:-1])
                   .reset_index(drop=True)
                   .join(df_shapes.groupby('id').apply(lambda x: x.iloc[1:])
                         .reset_index(drop=True),
                         rsuffix='2'))[['id', 'vertex_i', 'vertex_i2',
                                        'x', 'y', 'x2', 'y2']]
    v = (df_segments[['x2', 'y2']].values - df_segments[['x', 'y']]).values
    mid = .5 * v + df_segments[['x', 'y']].values
    x_mid = mid[:, 0]
    y_mid = mid[:, 1]
    length = np.sqrt((v ** 2).sum(axis=1))
    v_scaled = distance_threshold_px.magnitude * v / length[:, None]
    x_normal = -v_scaled[:, 1]
    y_normal = v_scaled[:, 0]

    # Create new data frame from scratch and join it to the `df_segments`
    # frame since it is **much** faster than adding new columns directly
    # the existing `df_segments` frame.
    df_normal = pd.DataFrame(np.column_stack([x_mid, y_mid, length, x_normal,
                                              y_normal]),
                             columns=['x_mid', 'y_mid', 'length', 'x_normal',
                                      'y_normal'])
    return df_segments.join(df_normal).set_index(['id', 'vertex_i'])


def get_intersections(df_segments, p, r):
    # See: https://stackoverflow.com/a/565282/345236
    q = df_segments[['x', 'y']].values
    s = df_segments[['x2', 'y2']].values - q

    r_x_s = np.cross(r, s)
    r_x_s[r_x_s == 0] = np.NaN
    t = np.cross((q - p), s) / r_x_s
    u = np.cross((q - p), r) / r_x_s

    df_tu = pd.DataFrame(np.column_stack(
        [t, u]), columns=list('tu'), index=df_segments.index)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        df_i = df_segments.join(df_tu).loc[(r_x_s != 0)
                                           & (t >= 0) & (t <= 1) & (u >= 0)
                                           & (u <= 1)]
    intersect_points = p + df_i.t.values[:, None] * r
    return df_i.join(pd.DataFrame(intersect_points, columns=['x_intersect',
                                                             'y_intersect'],
                                  index=df_i.index)).drop(['t', 'u'], axis=1)


def get_all_intersections(df_shapes,
                          distance_threshold=DEFAULT_DISTANCE_THRESHOLD):
    '''
    Parameters
    ----------
    svg_source : str or file-like or pandas.DataFrame
        File path, URI, or file-like object for SVG device file.

        If specified as ``pandas.DataFrame``, assume argument is in format
        returned by :func:`svg_model.svg_shapes_to_df`.
    distance_threshold : pint.quantity.Quantity
        Maximum gap between electrodes to still be considered neighbours.
    '''
    df_segments = get_segments(df_shapes,
                               distance_threshold=distance_threshold)

    intersections = []
    for i, ((id_i, vertex_i), segment_i) in enumerate(df_segments.iterrows()):
        p = segment_i[['x_mid', 'y_mid']].values
        r = segment_i[['x_normal', 'y_normal']].values

        df_intersections_i = get_intersections(df_segments, p, r)

        # Do not include self electrode in consideration for neighbours.
        self_mask = df_intersections_i.index.get_level_values('id') == id_i
        df_intersections_i = df_intersections_i.loc[~self_mask]
        if df_intersections_i.shape[0]:
            intersections.append(((id_i, vertex_i), df_intersections_i))
    index, values = zip(*intersections)
    df_result = pd.concat(values, keys=index)
    df_result.index.names = ['id', 'vertex_i',
                             'id_neighbour', 'vertex_i_neighbour']
    return df_result


def _get_neighbours(svg_source, distance_threshold=DEFAULT_DISTANCE_THRESHOLD,
                    **kwargs):
    '''
    Parameters
    ----------
    svg_source : str or file-like or pandas.DataFrame
        File path, URI, or file-like object for SVG device file.

        If specified as ``pandas.DataFrame``, assume argument is in format
        returned by :func:`svg_model.svg_shapes_to_df`.
    distance_threshold : pint.quantity.Quantity
        Maximum gap between electrodes to still be considered neighbours.

    Returns
    -------
    pandas.DataFrame
    '''
    df_shapes = _resolve_source(svg_source, **kwargs)
    df_segments = get_segments(df_shapes,
                               distance_threshold=distance_threshold)
    df_intersections = \
        get_all_intersections(df_shapes, distance_threshold=distance_threshold)

    df_neighbours = (df_intersections.reset_index([2, 3])
                     .join(df_segments, lsuffix='_neighbour'))
    df_neighbours.reset_index('id', inplace=True)
    df_neighbours.drop_duplicates(['id', 'id_neighbour'], inplace=True)
    df_neighbours.insert(0, 'direction', None)

    # Assign direction labels
    vertical = df_neighbours.x_normal.abs() < df_neighbours.y_normal.abs()
    df_neighbours.loc[vertical & (df_neighbours.y_normal < 0),
                      'direction'] = 'up'
    df_neighbours.loc[vertical & (df_neighbours.y_normal > 0),
                      'direction'] = 'down'
    df_neighbours.loc[~vertical & (df_neighbours.x_normal < 0),
                      'direction'] = 'left'
    df_neighbours.loc[~vertical & (df_neighbours.x_normal > 0),
                      'direction'] = 'right'
    df_neighbours.insert(0, 'normal_magnitude',
                         df_neighbours[['x_normal', 'y_normal']].abs()
                         .max(axis=1))

    df_neighbours.sort_values(['id', 'direction', 'normal_magnitude'],
                              inplace=True, ascending=False)
    # If multiple neighbours match a direction, only keep the first match.
    df_neighbours.drop_duplicates(['id', 'direction'], inplace=True)
    return df_neighbours


def get_id_neighbours(svg_source,
                      distance_threshold=DEFAULT_DISTANCE_THRESHOLD, **kwargs):
    # Make local copy with new index.
    df_neighbours = _get_neighbours(svg_source,
                                    distance_threshold=distance_threshold,
                                    **kwargs)
    ids = df_neighbours['id'].sort_values().drop_duplicates().values
    df = df_neighbours.set_index(['id', 'direction'])
    df.sort_index(inplace=True)

    directions = ['up', 'down', 'left', 'right']
    id_neighbours = (df.loc[[i for c in ids
                            for i in zip(it.cycle([c]), directions)],
                            'id_neighbour'])
    # XXX Work around Pandas regression where index names do not persist to
    # data frame view.
    id_neighbours.index.names = df.index.names
    id_neighbours


def get_channel_neighbours(svg_source,
                           distance_threshold=DEFAULT_DISTANCE_THRESHOLD,
                           **kwargs):
    '''
    Parameters
    ----------
    svg_source : str or file-like or pandas.DataFrame
        File path, URI, or file-like object for SVG device file.

        If specified as ``pandas.DataFrame``, assume argument is in format
        returned by :func:`svg_model.svg_shapes_to_df`.
    distance_threshold : pint.quantity.Quantity
        Maximum gap between electrodes to still be considered neighbours.

    Returns
    -------
    pandas.Series
    '''
    df_shapes = _resolve_source(svg_source, **kwargs)

    # Make local copy with new index.
    df_neighbours = _get_neighbours(df_shapes,
                                    distance_threshold=distance_threshold)

    electrode_channels = _extract_electrode_channels(df_shapes)
    df_neighbours.insert(0, 'channel',
                         electrode_channels.reindex(df_neighbours['id']).values)
    df_neighbours.insert(0, 'channel_neighbour',
                         electrode_channels
                         .reindex(df_neighbours['id_neighbour']).values)
    df_neighbours.set_index(['channel', 'direction'], inplace=True)
    df_neighbours.sort_index(inplace=True)

    directions = ['up', 'down', 'left', 'right']
    channel_neighbours = (df_neighbours.loc[[i for c in range(120)
                                             for i in zip(it.cycle([c]),
                                                          directions)],
                                            'channel_neighbour'])
    # XXX Work around Pandas regression where index names do not persist to
    # data frame view.
    channel_neighbours.index.names = df_neighbours.index.names

    return channel_neighbours


def chip_info(svg_source, **kwargs):
    '''
    Parameters
    ----------
    svg_source : `str` or file-like or `pandas.DataFrame`
        File path, URI, or file-like object for SVG device file.

        If specified as ``pandas.DataFrame``, assume argument is in format
        returned by :func:`svg_model.svg_shapes_to_df`.

    Returns
    -------
    `dict`
        Chip info with fields:

        - ``electrode_shapes``: ``x``, ``y``, ``width``, ``height``, and
          ``area`` for each electrode (`pandas.DataFrame`).
        - ``electrode_channels``: mapping from channel number to electrode ID
          (`pandas.Series`).
        - ``channel_electrodes``: mapping from electrode ID to channel number
          (`pandas.Series`).
    '''
    df_shapes = _resolve_source(svg_source, **kwargs)

    electrode_shapes = svg_model.data_frame.get_shape_infos(df_shapes, 'id')
    electrode_channels = _extract_electrode_channels(df_shapes)
    channel_electrodes = pd.Series(electrode_channels.index,
                                   index=electrode_channels.values)

    return {'electrode_shapes': electrode_shapes,
            'electrode_channels': electrode_channels,
            'channel_electrodes': channel_electrodes}


def _neighbours(df_electrode_shapes, df_connections):
    '''
    Parameters
    ----------
    df_electrode_neighbours : pandas.DataFrame
        Data frame indexed by electrode SVG shape "id", containing (at least)
        the following columns::

         - ``x``, ``y``: top-left coordinates of SVG shape (note that
           ``y``-axis is reversed in SVG, with 0 indicating top of image).
         - ``width``, ``height``: width and height of shape bounding box
    df_electrode_neighbours : pandas.DataFrame
        Each row corresponds to connection between two electrode SVG shapes
        denoted with respective column names, ``source`` and ``target``.

    Returns
    -------
    pandas.DataFrame
        Table indexed by electrode identifier (e.g., ``electrode001``) with
        columns, each indicating the identifier of the neighbouring electrode
        in the respective direction (``NaN`` if no neighbour in corresponding
        direction):
        - up
        - down
        - left
        - right

    Example
    -------

                                up          down          left         right
        id
        electrode000  electrode001           NaN  electrode043  electrode002
        electrode001  electrode088  electrode000  electrode043  electrode006
        electrode002  electrode006  electrode037  electrode000  electrode005
        electrode003  electrode063  electrode042           NaN           NaN
        electrode004  electrode066  electrode009  electrode068  electrode009
        ...

    '''

    df_centers = (df_electrode_shapes[['x', 'y']] +
                  .5 * df_electrode_shapes[['width', 'height']].values)
    df_by_source = df_connections.set_index('source')
    df_by_target = (df_connections.set_index('target')
                    .rename(columns={'source': 'target'}))
    df_neighbours = df_by_source.append(df_by_target)
    df_neighbours.index.name = 'source'

    # Add **source** x/y center coordinates
    df_neighbours = df_neighbours.join(df_centers.loc[df_neighbours.index
                                                      .drop_duplicates()])

    # Add **target** x/y center coordinates
    df_target_centers = df_centers.loc[df_neighbours.target]
    df_neighbours['target_x'] = df_target_centers.x.values
    df_neighbours['target_y'] = df_target_centers.y.values

    df_neighbours['x_delta'] = df_neighbours.target_x - df_neighbours.x
    df_neighbours['y_delta'] = df_neighbours.target_y - df_neighbours.y

    # Index by target electrode identifier.
    df_neighbours.reset_index(inplace=True)
    df_neighbours.set_index('target', inplace=True)

    # Find neighbour in each direction.
    up = df_neighbours.loc[df_neighbours.y_delta <
                           0].groupby('source')['y_delta'].idxmin()
    down = df_neighbours.loc[df_neighbours.y_delta >
                             0].groupby('source')['y_delta'].idxmax()
    left = df_neighbours.loc[df_neighbours.x_delta <
                             0].groupby('source')['x_delta'].idxmin()
    right = df_neighbours.loc[df_neighbours.x_delta >
                              0].groupby('source')['x_delta'].idxmax()

    df_electrode_neighbours = pd.DataFrame(None,
                                           index=df_electrode_shapes.index)
    df_electrode_neighbours['up'] = up
    df_electrode_neighbours['down'] = down
    df_electrode_neighbours['left'] = left
    df_electrode_neighbours['right'] = right
    return df_electrode_neighbours
